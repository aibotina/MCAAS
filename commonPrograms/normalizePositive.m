function [feature, normVec] = normalizePositive(trainFeature, normalization)
%NORMALIZE normalizes the trainFeature according to normalization
%
% USAGE:
%       [feature, normVec] = normalize(trainFeature, normalization)
% INPUT:
%       trainFeature: a feature matrix
%       normalization: 1-normalize by infinity norm; 2- by 2-norm
% OUTPUT:
%      feature: the normalized feature
%      normVec: the normalization vector in each column
% 
% By Jialiang Yang MSU

%% normalization process
[nRow, nCol] = size(trainFeature);
normVec = zeros(1,nCol);
feature = zeros(nRow, nCol);

for i = 1: nCol
    vec = trainFeature(:, i);
    
    if normalization == 1
        normVec(i) = norm(vec, Inf);
    else
        normVec(i) = norm(vec, 2);
    end
    
    feature(:,i) = vec/normVec(i);
end

end

