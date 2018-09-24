clc;clear;
path('../similarityMatrix',path);
fileName = 'H3N2-68-03.tab';
thresh = 10; 
norm=1;
% rank =10;
iter = 2000;
tol = 1e-8;
lam1=1e-4;
lam2=3.28e-13;
lam3=3.28e-13;

% [recMatrix,virusName] =matrixCompletion3(fileName, thresh, rank, iter, tol, norm,lam1,lam2);
% [dataHI, virusName, serumName, reference] = readTable(fileName);
%Load and normalize files
[dataNorm,typeI, typeII,virusName, serumName] = normalizeTable(fileName,thresh,norm);

[nVirus,nSerum]=size(dataNorm);

%Year of virus extraction
for i=1:nVirus
   virusName=strtrim(virusName); 
   aa{i}= virusName{i}(end-1:end);
   virusNameYear(i)=str2num(aa{i});
   if virusNameYear(i)>10
       virusNameYear(i)=virusNameYear(i)+1900;
   else virusNameYear(i)=virusNameYear(i)+2000;
   end
end
virusNameYear=virusNameYear';

%Year of serum extraction
for j=1:nSerum
   serumName=strtrim(serumName);
   bb{j}= serumName{j}(end-1:end);
   serumNameYear(j)=str2num(bb{j});
    if serumNameYear(j)>10
       serumNameYear(j)=serumNameYear(j)+1900;
   else serumNameYear(j)=serumNameYear(j)+2000;
   end   
end
%Different type values in HI
typeIcol_01=typeI(:,1);
typeIcol_02=typeI(:,2);
typeI_E=zeros(nVirus,nSerum);
typeI_E(sub2ind(size(typeI_E),typeIcol_01,typeIcol_02))=1;

typeIIcol_11=typeII(:,1);
typeIIcol_22=typeII(:,2);
typeII_E=zeros(nVirus,nSerum);
typeII_E(sub2ind(size(typeII_E),typeIIcol_11,typeIIcol_22))=1;

% Get sliding window
Size=32;
Count=37-Size;
load SeqSM_500_2_Prlic2000.mat
load SeqSMsera_500_2_Prlic2000.mat
for i=1:Count
    k=i+1967;
    index1=find(virusNameYear>=k & virusNameYear<=k+Size-1 );
    index2=find(serumNameYear>=k & serumNameYear<=k+Size-1 );
    SubdataNorms{i}=dataNorm(index1,index2);
    SubSeqSMs{i}=SeqSM(index1,index1);
    SubSeqSMseras{i}=SeqSMsera(index2,index2);
    
    SubtypeI_Es{i}=typeI_E(index1,index2);
    SubtypeII_Es{i}=typeII_E(index1,index2);
    [m,n]=find(SubtypeII_Es{i}'==1);
    SubtypeIIs{i}=[n,m];
end


%Iterate by window swipe
rank_store=9;
Mid=0;
for KK=1:numel(rank_store)
    rank=rank_store(KK)
   rr=rank-1;
   SumUnMatrix=zeros(nVirus,nSerum);
   Times=zeros(nVirus,nSerum);
for j=1:Count
    [X S Y dist] = OptSpaceII_sera(sparse(SubdataNorms{j}),SubtypeIIs{j}, SubSeqSMs{j},SubSeqSMseras{j},rr,iter,tol,lam1,lam2,lam3);
    recMatrixes{j} = X*S*Y';
    
    Tent1_1=zeros(nVirus,nSerum);
    Tent2_2=zeros(nVirus,nSerum);
    
    k=j+1967;
    index1=find(virusNameYear>=k & virusNameYear<=k+Size-1);
    index2=find(serumNameYear>=k & serumNameYear<=k+Size-1);
    
    Tent1_1(index1,index2)=recMatrixes{j};
    SumUnMatrix=SumUnMatrix+Tent1_1;
    
    [p,q]=size(recMatrixes{j});
    Tent2_2(index1,index2)=ones(p,q);
    Times=Times+Tent2_2;
end
for j=1:Count
    [X S Y dist] = OptSpaceII_sera(sparse(SubdataNorms{j}),SubtypeIIs{j}, SubSeqSMs{j},SubSeqSMseras{j},rank,iter,tol,lam1,lam2,lam3);
    recMatrixes{j} = X*S*Y';
    
    Tent1_1=zeros(nVirus,nSerum);
    Tent2_2=zeros(nVirus,nSerum);
    
    k=j+1967;
    index1=find(virusNameYear>=k & virusNameYear<=k+Size-1);
    index2=find(serumNameYear>=k & serumNameYear<=k+Size-1);
    
    Tent1_1(index1,index2)=recMatrixes{j};
    SumUnMatrix=SumUnMatrix+Tent1_1;
    
    [p,q]=size(recMatrixes{j});
    Tent2_2(index1,index2)=ones(p,q);
    Times=Times+Tent2_2;
end

[p0,q0]=find(Times~=0);
RemianDataNorm=dataNorm;
RemianDataNorm(sub2ind(size(RemianDataNorm),p0,q0))=0;

for m1=1:nVirus
    for n1=1:nSerum
        if SumUnMatrix(m1,n1)~=0
            AverUnMatrix(m1,n1)=SumUnMatrix(m1,n1)/ Times(m1,n1);
        else AverUnMatrix(m1,n1)=0;
        end
    end
end
save AverUnMatrix;

UncomMatrix=RemianDataNorm+AverUnMatrix;

%10 fold cross validation
[crossmatrixs] = getcrossMatrixs(UncomMatrix);
E={};
 for j = 1:10
     E{j} = full(spones(crossmatrixs{j}));
 end   
 for j = 1:10
        trainmatrix{j}=UncomMatrix-crossmatrixs{j};
 end  
 for j=1:10
     EE{j} = full(spones(trainmatrix{j}));  
      PossII= typeII_E.* EE{j};     
     [m0,n0]=find(PossII'==1);
     Pos=[n0,m0];
     triantypeIIs{j}=Pos;   
 end

 SUM=0;
  for k = 1:10
   [X S Y dist] = OptSpaceII_sera(sparse(trainmatrix{k}),triantypeIIs{k},SeqSM,SeqSMsera,rank,iter,tol,lam1,lam2,lam3);
    recMatrix = X*S*Y';
    PossI=E{k}.* typeI_E;
    fmatrix=(dataNorm-recMatrix).*PossI;
    RMSE=sqrt(sum(sum(fmatrix.^2))/nnz(PossI));
    SUM=SUM+RMSE;
    end 
 AverageRMSE=SUM/10;
 Mid=Mid+1;
 Total(Mid,1)=rank;
 Total(Mid,2)=AverageRMSE;
end
 save Total_500_2_Prlic2000;    









