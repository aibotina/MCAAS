function  [xOptimal, normVec, xConst, constNorm] = lpfast(A, b, lambda,normalization)
%% lpfast solves the regression problem with positive coefficients by lasso.
%    USAGE:
%        [xOptimal, normVec, xConst, constNorm] = lpfast(A, b, lambda,normalization)
%
%    INPUT: 
%       [Compulsary]
%          A: the prediction Matrix
%          b: the response vector
%       [Optional]
%          lambda: the lasso parameter to control the sparsity
%         normalization: 0 no normalization 1 normalize by inf-norm 2 normalize by 2-norm
%
%    OUTPUT:
%         xOptimal: the positive weights
%         normVec: the value divided in each column.
%         xConst:  the constant term in linear approximation
%         constNorm: the value divided in the constant vec.
%    Jialiang Yang: jyang@cvm.msstate.edu

%% initialization
% add the constant term to prediction matrix A.
[rA, cA] = size(A);

oneA = ones(rA,1);
A = [oneA A];

% specify model parameters
tol = 1e-6;

x0 = ones(cA+1,1);

L = 1000000;
maxit = 100000;

%% normaize the prediction matrix
if normalization == 0
    normA = A;
    
    normVec = ones(1, cA+1);
else
    [normA, normVec] = normalizePositive(A, normalization);
end

for i = 1: maxit
    GradA = normA'*(normA*x0-b);
    
    x1 = max(x0 - 1/L* GradA - lambda/L, 0);
    
    if norm(x1-x0, 2) < tol
        break;
    end
       
    % update x0 and L0
    x0 = x1;
end

if i == maxit
    disp('Warning: maximum iteration reached, not convergent, Please add the value of L in lpfast and rerun the program!')
end

xConst = x1(1);
xOptimal = x1(2:cA+1);

constNorm = normVec(1);
normVec = normVec(2: cA+1);

end