function  [xOptimal, normVec, xConst, constNorm] = lassoPositive(A, b, lambda,normalization)
%% lassoPositive solves the regression problem with positive coefficients by lasso.
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
A = [oneA, A];

% specify model parameters
tol = 1e-6;

x0 = ones(cA+1,1);

L0 = 2;
maxit = 100000;

Lmin = 1.5;

if L0 < Lmin
    error('The initial choice of L0 is less than Lmin, try a larger one!');
end

rdec = 1.2; % decreasing factor

%% normaize the prediction matrix
if normalization == 0
    normA = A;
    
    normVec = ones(1, cA+1);
else
    [normA, normVec] = normalizePositive(A, normalization);
end

for i = 1: maxit
    [x1, M0] = lineSearch(normA, b, lambda, x0, L0);
    
    L0 = max(Lmin, M0/rdec);
    
    if norm(x1-x0, 2) < tol
        break;
    end
       
    % update x0 and L0
    x0 = x1;
end

if i == maxit
    disp('Warning: maximum iteration reached, not convergent!')
end

xConst = x1(1);
xOptimal = x1(2:cA+1);

constNorm = normVec(1);
normVec = normVec(2: cA+1);

end

% calculate the parameter (length step) for line search
function [xPlus, M] = lineSearch(A, b, lambda, x, L)

yinc = 1.1;

phi = 1;
psi = 0;

k = 0;
while (phi > psi)
    L = L*yinc;
    
    % calculate phi and psi
    GradA = A'*(A*x-b);
    
    xPlus = max(x - 1/L* GradA - lambda/L, 0);
    
    phi = 1/2*(A*xPlus-b)'*(A*xPlus-b) + lambda*norm(xPlus,1);
    
    psi = 1/2*(A*x-b)'*(A*x-b) + GradA'*(xPlus-x) + 1/2*L*(xPlus-x)'*(xPlus-x) + lambda*norm(xPlus,1);

    x = xPlus;
    
    k = k+1;
    
    if k > 10000
        break;
    end
end

M = L;
end