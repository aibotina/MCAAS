function [ distCol ] = getUpTriangle(distMatrix)
%%  getUpTriangle use the up-trianlge of distMatrix to make the pairwise distances vector
%
%    USAGE:
%       getUpTriangle(distFile)
%
%    Input:
%          distFile: distance Matrix
%    Output:
%          distCol: a vector stores the pairwise distances
%             
% Author: Jialiang Yang, CVM, MSU, jyang@cvm.msstate.edu
% Revision: 3/1/2012

rAln = size(distMatrix,1);

% use the up-trianlge of distMatrix to make the pairwise distances column 
rProf = rAln*(rAln-1)/2;
distCol = zeros(1,rProf);

k = 1;
for i = 1: rAln-1
    for j = i+1 : rAln
        distCol(k) = distMatrix(i, j);
        k = k+1;
    end
end

end

