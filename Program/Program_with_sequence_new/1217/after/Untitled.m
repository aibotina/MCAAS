%clear;clc;
a=[1 2 3 6;4 5 6 9;7 8 9 11;12 15 2 5];
Vec1=1;
Vec2=3
An=[];
aa=[1;2;3;4];
 b=a(Vec1,aa);
 c=a(Vec2,aa);
An=[An,b];
An=[An,c];


%% get the distance matrix for each cluster
pIDistAntige={};
k=numel(clusterVec);
for i=1:k
    Mid=virusCluster{i};
    m=size(Mid,1);
    Vec=Mid;
    distAn=[];
    for j=1:m
     %  An=distMatrix(Mid(j),Vec);
     An=Seq_pdis(Mid(j),Vec);
      distAn=[distAn;An];
    end
   pIDistAntige{i}= distAn;    
end
save('pIDistAntige.mat','pIDistAntige');


%% get the pairwise distance for the eleven cluters
k=numel(clusterVec);
for i=1:k
    m=size(virusCluster{i},1);
    for j=1:k
        n=size(virusCluster{j},1);
     %   Bn=distMatrix(virusCluster{i},virusCluster{j});
         Bn=Seq_pdis(virusCluster{i},virusCluster{j});
        pEDistAntige(i,j)=sum(sum(Bn))/(m*n);
    end
end 
save('pEDistAntige.mat','pEDistAntige');


% %rearrange the order
mm=numel(virusCluster);
RemdsCorr=[];
for b=1:mm
    vv=virusCluster{b};
    RemdsCorr=[RemdsCorr;mdsCorr(vv,:)];
end






