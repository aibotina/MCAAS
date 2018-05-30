function [X S Y dist] = OptSpaceII_sera(M_E,SubTypeII,SubSeqSM,SubSeqSMsera,r,niter,tol,lam1,lam2,lam3,E)
% An algorithm for Matrix Reconstruction from a partially revealed set. 
% See "Matrix Completion from a Few Entries"(http://arxiv.org/pdf/0901.3150) for details
% Usage :
% [X S Y dist] = OptSpace(A,SubTypeII,r,niter,tol,lam);
% [X S Y dist] = OptSpace(A,r,SubTypeII,niter,tol,lam,E);
% [X S Y dist] = OptSpace(A);
% 
% INPUT :
% A     :  The partially revealed matrix.
%          Sparse matrix with zeroes at the unrevealed indices.
% SubTypeII:  The positions of SubTypeII data 
%
% r     :  The rank to be used for reconstruction. Use [] to guess the rank.
% niter :  The max. no. of iterations. Use [] to use default (50).
% tol   :  Stop iterations if norm( (XSY' - M_E).*E , 'fro' )/sqrt(|E|) < tol, where
%        - E_{ij} = 1 if M_{ij} is revealed and zero otherwise, 
%        - |E| is the size of the revealed set.				
%        - Use [] to use the default (1e-6)
% lam	:  The coefficient of regularization (\lam*||XSY^T||_F^2)
% E     :  The matrix E such that : E_{ij} = 1 if entry (i,j) is revealed, 0 otherwise
%
% OUTPUT :
% X      : A size(A,1)xr matrix
% S      : An rxr matrix
% Y      : A size(A,2)xr matrix
% such that M_hat = X*S*Y' 
% dist   : A vector containing norm( (XSY' - M_E).*E , 'fro' )/sqrt(|E|) at each
%          successive iteration
%
% Date : 21st September, 2010
% COPYRIGHT 2009 Raghunandan H. Keshavan, Andrea Montanari, Sewoong Oh



if(nargin==1)
	
	M_E = sparse(M_E);
	[n m] = size(M_E);
	E = spones(M_E);
    SubTypeII = [];
	eps = nnz(E)/sqrt(m*n) ;

	tol = 1e-8;
	lam1 = 0.0;
    lam2 = 0.0;
	fprintf(1,'Rank not specified. Trying to guess ...\n');
	r = guessRank(M_E) ;
	fprintf(1,'Using Rank : %d\n',r);
	
	m0 = size(M_E,1);
    niter = 50;
elseif(nargin==8 || nargin==9||nargin==10||nargin==11)
	
	M_E = sparse(M_E);
	[n m] = size(M_E);

	if(nargin==8||nargin==10)
		E = spones(M_E);
    end

	eps = nnz(E)/sqrt(m*n) ;

	if( length(tol) == 0 )
		tol = 1e-5;
	end

	if( length(r) == 0 )

		fprintf(1,'Rank not specified. Trying to guess ...\n');
		r = guessRank(M_E) ;
		fprintf(1,'Using Rank : %d\n',r);
	end

	m0 = size(M_E,1);
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
%S=S0;
%S = getoptS1(X,S,Y,M_E,E,lam2,lam3,SubSeqSM,SubSeqSMsera);
S = getoptS_new(X,Y,M_E,E,lam1,lam2,lam3,SubSeqSM,SubSeqSMsera) ;
%S = getoptS(X,Y,M_E,E,lam1) ;
matrix1=X*S*Y';

dist(1) = norm( (M_E - matrix1).*E ,'fro')/sqrt(nnz(E) )  ;
fprintf(1,'0\t\t%e\n',dist(1) ) ;

% backup  E
EOrigin = E;

for i = 1:niter
    % restore the orignal M_E and E
    E = EOrigin;
    
    % if the recovered the SubTypeII value is already less than the threshold,
    % ignore the punishment.
    recMatrix = X*S*Y';
   
    for j = 1: size(SubTypeII,1)
        if recMatrix(SubTypeII(j,1),SubTypeII(j,2)) < M_E(SubTypeII(j,1),SubTypeII(j,2))-1
            E(SubTypeII(j,1),SubTypeII(j,2)) = 0;
        end
    end

% Compute the Gradient 
	%[W Z] = gradF_t(X,Y,S,M_E,E,m0,lam1,lam2,SubSeqSM);
     [W Z] = gradF_t(X,Y,S,M_E,E,m0,lam1,lam2,SubSeqSM,lam3,SubSeqSMsera);

% Line search for the optimum jump length	
	%t = getoptT(X,W,Y,Z,S,M_E,E,m0,lam1,lam2,SubSeqSM) ;
    t = getoptT(X,W,Y,Z,S,M_E,E,m0,lam1,lam2,SubSeqSM,lam3,SubSeqSMsera) ;
	X = X + t*W;
    Y = Y + t*Z;
    S = getoptS_new(X,Y,M_E,E,lam1,lam2,lam3,SubSeqSM,SubSeqSMsera) ;
   % S = getoptS(X,Y,M_E,E,lam1) ;
	
% Compute the distortion	
	dist(i+1) = norm( (recMatrix - X*S*Y').*E,'fro' )/sqrt(nnz(EOrigin));
%	fprintf(1,'%d\t\t%e\n',i,dist(i+1) ) ;
	if( dist(i+1) < tol )
		break ;
	end
end

S = S /rescal_param ;
end