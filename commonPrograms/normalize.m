function [X,mx,stdx] = normalize(X)
% NORMALIZE  Normalize the observations of a data matrix.
%
%    [X,mx,stdx] = NORMALIZE(X) centers and scales the observations of a 
%    data matrix such that each variable (column) has unit length.
%
% Author: Karl Skoglund, IMM, DTU, kas@imm.dtu.dk

n = length(X(:,1));
mx = mean(X);
stdx = std(X,0,1);

idx = find(abs(stdx) < sqrt(eps(class(stdx)))); 
if any(idx)
  stdx(idx) = 1;
end

MX = mx(ones(n,1),:);
STDX = stdx(ones(n,1),:);
Z = (X - MX) ./ STDX;
X=Z;
