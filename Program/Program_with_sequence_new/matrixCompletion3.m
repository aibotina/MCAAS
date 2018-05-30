function [recMatrix,virusName] = matrixCompletion3(fileName, thresh, rank, iter, tol, norm,lam1,lam2)
% check fileName
% check if folder temp exits
pwd;
folderName = pwd;

if ~exist('Temp', 'dir')
    mkdir('Temp')
end

outFile = [folderName, '/Temp/Complete.txt'];
[dataNorm, typeII,virusName, serumName] = normalizeTable(fileName,thresh,norm);
[X S Y dist] = OptSpaceII(sparse(dataNorm),typeII,rank,iter,tol,lam1,lam2,rho);
recMatrix = X*S*Y'; 
end
%==========================================================================
function [dataNorm,typeOne, typeTwo,virusName, serumName] = normalizeTable(fileName, thresh, norm)
if nargin < 2
    help normalizeTable
    return
    
elseif nargin <3
    norm =1;  % specifies normalization scheme 
    
end
% check fileName
if ~ischar(fileName)
    error('FileName must be a character array!')
end
%% Read data from file fileName
[dataHI, virusName, serumName, reference] = readTable(fileName);

[nVirus, nSerum] = size(dataHI);

% tableHI stors the numerical HI table: change typeII data to threshold
%-------------------------------------------------------------------------
tableHI = zeros(nVirus, nSerum);
typeOne=[];
typeTwo = [];

for i = 1: nVirus
     for j = 1: nSerum
         HIValue = strtrim(dataHI{i,j});
         
         % missing value
         if strcmp(HIValue,'0')
             tableHI(i,j) = 0;
             
         % lower reactor like '<20', take tableHI(i,j) = 20
         elseif strcmp(HIValue(1),'<') 
             tableHI(i,j) = str2double(HIValue(2: numel(HIValue))) ;
             %tableHI(i,j) = thresh;
             
             % update typeTwo matrix
             pos = [i,j];
             typeTwo = [typeTwo; pos];
             
         else
             value = str2double(HIValue);
             
             % lower reactor like 25 < 40 take tableHI(i,j) = 25 
             if value < thresh
                 tableHI(i,j) = value;
                 %tableHI(i,j) = thresh;
                 
                  % update typeTwo matrix
                 pos = [i,j];
                 typeTwo = [typeTwo; pos];
                 
             % type I value
             else
                 tableHI(i,j) = value;
                 pot=[i,j];
                 typeOne=[typeOne;pot];
             end
         end
     end
end
dataNorm = zeros(nVirus, nSerum);  % stores the data after normalization

switch norm
    case 1
        % identify max overall and max in each column.
        maxCol = max(tableHI);
        
        logCol = log2(maxCol);  % log2(max(Hj))
        
        % normalize typeI and type II entries by
        % log2(max(Hij))-log2(max(Hj))+log2(Hij) (or log2(thresh))
        for i = 1: nVirus
            for j = 1: nSerum
                % not missing value
                if tableHI(i,j) ~=0
                    dataNorm(i,j) = logCol(j) - log2(tableHI(i,j));
                end
            end
        end
        
        % pay attention: artificially + 1 to fit zhipeng's data
        logMax = max(max(dataNorm));

        logMax = ceil(logMax)+1;
        
        for i = 1: nVirus
            for j = 1: nSerum
                % not missing value
                if tableHI(i,j) ~=0
                    dataNorm(i,j) = logMax - dataNorm(i,j);
                end
            end
        end
        
        %% add the code for other normalization schemes
    case 2
        dataNorm = tableHI;
end

pwd;
currentFolder = pwd;

if ~exist('Temp', 'dir')
    mkdir('Temp')
end

normalizedTableFile= [currentFolder, '/Temp/normalizedTabe.txt'];
writeNormalTable(normalizedTableFile, dataNorm, virusName, serumName);
end

%==========================================================================
function writeNormalTable(fileName, dataHI, virusName, serumName)
% Parameters fileName, dataHI, virusName, serumName should be specified
if nargin == 0
    help writeNormalTable
    return
elseif nargin < 4
    error('Parameters: fileName dataHI virusName serumName')   
end

if ~ischar(fileName)
    error('FileName must be a character array!')
end

% Check the validality of the input data 
[dataRow,dataCol] = size(dataHI);
nVirus = length(virusName);
nSerum = length(serumName);
nSerum = dataCol;

% Check viruses
if dataRow ~= nVirus
    error('The number of viruses is differnt in dataHI and virusName!')
end

% Check serum
% if dataCol ~= nSerum
%     error('The number of serums is differnt in dataHI and serumName!')
% end
    
%% Open file fileName for writing

% Open file fileName for writing
fid = fopen(fileName,'w+');

if fid == (-1)
    error('Can not open file %s for writing.\n Check wirte permission!',...
         fileName)
end

% Write data to fileName

% Virus number and Serum number line
fprintf(fid, '%d\t%d\n', nVirus, nSerum);

% Serum line
fprintf(fid,'ID\t');

for i = 1: nSerum-1
    fprintf(fid,'%s\t', serumName{i});
end
fprintf(fid,'%s\n', serumName{nSerum});

% Reaction table
if isnumeric(dataHI)
    
    for i = 1: nVirus
        
        % Virus name
        fprintf(fid,'%s \t', virusName{i});
        
        % table
        for j = 1: nSerum-1
            fprintf(fid,'%d \t', dataHI(i,j));
        end
        
        fprintf(fid,'%d \n', dataHI(i,nSerum));
    end
else
    for i = 1: nVirus
        % Virus name
        fprintf(fid,'%s \t', virusName{i});
        
        % table
        for j = 1: nSerum-1
            fprintf(fid,'%s \t', dataHI{i,j});
        end
        
        fprintf(fid,'%s \n', dataHI{i,nSerum});
    end
end

fclose(fid);

end

%==========================================================================
function [X S Y dist] = OptSpaceII(M_E,typeII,r,niter,tol,lam1,lam2,rho,E)
if(nargin==8 || nargin==9)
	
	M_E = sparse(M_E);
	[n m] = size(M_E);

	if(nargin==8)
		E = spones(M_E);
    end

	eps = nnz(E)/sqrt(m*n) ;

	if( length(tol) == 0 )
		tol = 1e-6;
	end

	if( length(r) == 0 )

		fprintf(1,'Rank not specified. Trying to guess ...\n');
		r = guessRank(M_E) ;
		fprintf(1,'Using Rank : %d\n',r);
    end
    
    m0 = size(M_E,1);
	rho = 0.01;	
    
	if( length(niter) == 0 )
		niter = 50 ;
	end	
else
	fprintf(1,'Improper arguments (See "help OptSpace")\n');
	fprintf(1,'Usage :\n[X S Y dist] = OptSpace(A,r,niter,tol) \n') ;
	fprintf(1,'[X S Y dist] = OptSpace(A)\n');
	return;
end	
% show frobinus norm of M_E
% norm(M_E,'fro')^2
rescal_param = sqrt( nnz(E) * r / norm(M_E,'fro')^2 ) ;

M_E = M_E * rescal_param ;

fprintf(1,'Trimming ...\n');
% Trimming

M_Et = M_E ;
d = sum(E);


d_=mean(full(d));

rand('seed', 10);%rng(12345);

for col=1:m
    if ( sum(E(:,col))>2*d_ )
        list = find( E(:,col) > 0 );
        p = randperm(length(list));
        M_Et( list( p(ceil(2*d_):end) ) , col ) = 0;
    end
end

% update E to EE
EE = spones(M_Et);

d = sum(EE');
d_= mean(full(d));


for row=1:n
    if ( sum(EE(row,:))>2*d_ )
        list = find( EE(row,:) > 0 );
        p = randperm(length(list));
        M_Et(row,list( p(ceil(2*d_):end) ) ) = 0;
    end
end

fprintf(1,'Sparse SVD ...\n');
% Sparse SVD

[X0 S0 Y0] = svds(M_Et,r) ;

clear M_Et;

% Initial Guess
X0 = X0*sqrt(n) ; Y0 = Y0*sqrt(m) ;
S0 = S0 / eps ;

fprintf(1,'Iteration\tFit Error\n');

% Gradient Descent
X = X0;Y=Y0;

S = getoptS(X,Y,M_E,E,lam1);
matrix1=X*S*Y';

dist(1) = norm( (M_E - matrix1).*E ,'fro')/sqrt(nnz(E) )  ;
fprintf(1,'0\t\t%e\n',dist(1) ) ;

% backup  E
EOrigin = E;
for i = 1:niter
    % restore the orignal M_E and E
    E = EOrigin;
    
    % if the recovered the typeII value is already less than the threshold,
    % ignore the punishment.
    recMatrix = X*S*Y';
   
    for j = 1: size(typeII,1)
        if recMatrix(typeII(j,1),typeII(j,2)) < M_E(typeII(j,1),typeII(j,2))-1
            E(typeII(j,1),typeII(j,2)) = 0;
        end
    end

% Compute the Gradient 
	[W Z] = gradF_t(X,Y,S,M_E,E,m0,lam1,lam2,rho);
    % Line search for the optimum jump length	
	t = getoptT(X,W,Y,Z,S,M_E,E,m0,lam1,lam2,rho) ;
	X = X + t*W;
    Y = Y + t*Z;
    S = getoptS(X,Y,M_E,E,lam1) ;
	
% Compute the distortion	
	dist(i+1) = norm( (recMatrix - X*S*Y').*E,'fro' )/sqrt(nnz(EOrigin));
	fprintf(1,'%d\t\t%e\n',i,dist(i+1) ) ;
	if( dist(i+1) < tol )
		break ;
	end
end

S = S /rescal_param ;
end

% Function to Guess the Rank of the input Matrix
function r = guessRank(M_E)
	[n, m] = size(M_E);
	epsilon = nnz(M_E)/sqrt(m*n);

    S0 = svds(M_E,100) ;

    S1=S0(1:end-1)-S0(2:end);
    S1_ = S1./mean(S1(end-10:end));
    r1=0;
    lam3=0.05;
    while(r1<=0)
        for idx=1:length(S1_)
            cost(idx) = lam3*max(S1_(idx:end)) + idx;
        end
        [v2,i2] = min(cost);
        r1 = max(i2-1);
        lam3=lam3+0.05;
    end

	clear cost;
    for idx=1:length(S0)-1
        cost(idx) = (S0(idx+1)+sqrt(idx*epsilon)*S0(1)/epsilon  )/S0(idx);
    end
    [v2,i2] = min(cost);
    r2 = max(i2);

	r = max([r1 r2]);
end


% Function to compute the distortion
function out = F_t(X,Y,S,M_E,E,m0,lam1,lam2,rho)
[n,r] = size(X) ;
[m,r] = size(Y) ;
% SeqSM is the HA pro sequence similary matrix
load SeqSM;
s = SeqSM;
a = X*S*Y';
[m,n] = size(s);
Sum = 0;
for i = 2:m
     for j = 1:i-1
         Sum = Sum+s(i,j)*norm((a(i,:)-a(j,:)),2).^2;
     end
end
out4 =lam2*Sum;
out1 = sum( sum( ( (X*S*Y' - M_E).*E ).^2 ) )/2 ;
out2 =  rho*G(Y,m0,r) ;
out3 =  rho*G(X,m0,r) ;
out = out1+out2+out3+out4 + lam1*m*n*norm(S, 'fro')^2;
end
function out = G(X,m0,r)
z = sum(X.^2,2)/(2*m0*r) ;
y = exp( (z-1).^2 ) - 1 ;
y( find(z < 1) ) = 0 ;
out = sum(y) ;
end

% Function to compute the gradient
function [W, Z] = gradF_t(X,Y,S,M_E,E,m0,lam1,lam2,rho)
[n, r] = size(X);
[m, r] = size(Y);
XS = X*S ;
YS = Y*S' ;
XSY = XS*Y' ;
load SeqSM;
B = GetseqM(SeqSM);
G3x = B*XSY*YS;
G3y = YS*X'*B*XS;
Qx = X'* ( (M_E - XSY).*E )*YS /n;
Qy = Y'* ( (M_E - XSY).*E )'*XS /m;
W = ( (XSY - M_E).*E )*YS + X*Qx + rho*Gp(X,m0,r)+4*lam2*G3x;
Z = ( (XSY - M_E).*E )'*XS + Y*Qy + rho*Gp(Y,m0,r)+4*lam2*G3y;
end

function out = Gp(X,m0,r)
z = sum(X.^2,2) /(2*m0*r) ;
z = 2*exp( (z-1).^2 ).*(z-1) ;
z( find(z<0) ) = 0;

out = X.*repmat(z,1,r) / (m0*r) ;
end

% Function to find Sopt given X, Y
function out = getoptS(X,Y,M_E,E,lam1)

[n, r] = size(X);
[m, r] = size(Y);
C = X' * ( M_E ) * Y ; C = C(:) ;

for i = 1:r
        for j = 1:r
                ind = (j-1)*r + i ;
                temp = X' * (  (X(:,i) * Y(:,j)').*E ) * Y ;
                A(:,ind) = temp(:) ;
				A(ind, ind) = A(ind, ind) + lam1*m*n ;
        end
end

S = A\C ;
out = reshape(S,r,r) ;
end

% Function to perform line search
function out = getoptT(X,W,Y,Z,S,M_E,E,m0,lam1,lam2,rho)
norm2WZ = norm(W,'fro')^2 + norm(Z,'fro')^2;
f(1) = F_t(X, Y,S,M_E,E,m0,lam1,lam2,rho) ;

t = -1e-1 ;
for i = 1:20
        f(i+1) = F_t(X+t*W,Y+t*Z,S,M_E,E,m0,lam1,lam2,rho) ;

        if( f(i+1) - f(1) <= .5*(t)*norm2WZ )
            out = t ;
            return;
        end
        t = t/2 ;
end
out = t ;
end
  