function  [distHamming, distEuc] = computeDist(similarMatrix, distType, virusCluster, serumCluster)
%% computeDist caluclate the pairwise distances for viruses in similarMatrix
%
%    Usage:
%          computeDist    gives this help on the computeDist function.
%          [distHamming, distEuc] = computeDist(similarMatrix, distType, virusCluster, serumCluster)
%    Input:
%          similarMatrix: similarity matrix
%          distType: 'Hamming' or 'Euclidean'
%          virusCluster: the clusters of viruses
%          serumCluster: the clusters of correspnding serum
%    Output:
%          distHamming: Hamming distance matrix
%          distEuc: Euclidean distance matrix
%
%  Revision Date : 5th, Jan 2011
%  Author: Jialiang Yang  CVM, MSU, jyang@cvm.msstate.edu

%% Input checking

if nargin < 1
    help computeDist
    return
    
elseif nargin < 2
    distType = 'Hamming';
    virusCluster = [];
    serumCluster = [];
    
elseif nargin < 3
    virusCluster = [];
    serumCluster = [];
    
elseif nargin < 4
    serumCluster = [];
    
end

nCluster = numel(virusCluster);
distEuc = [];

%% Calculate distances
rMatrix = size(similarMatrix,1);

% if serumCluster is empty, cacluate average distance
if isempty(serumCluster)
    
    % calculate the Hamming distance among viruses
    % -------------------------------------------------------------------------
    distHamming = zeros(rMatrix);
    
    % calculate hamming distances
    for i = 1: rMatrix-1
        for j = i+1: rMatrix
            distHamming(i,j) = mean(abs(similarMatrix(i,:)- similarMatrix(j,:)));
            distHamming(j,i) = distHamming(i,j);
        end
    end
   
    
    % calculate the scaled Eclidean distance
    if strcmpi(distType, 'Euclidean')
        
        distEuc = zeros(rMatrix);
        
        % calculate eucliean distances
        for i = 1: rMatrix-1
            for j = i+1: rMatrix
                distEuc(i,j) = sqrt(sum((similarMatrix(i,:)- similarMatrix(j,:)).^2));
                distEuc(j,i) = distEuc(i,j);
            end
        end
        
        ratio = sum(sum(distHamming))/sum(sum(distEuc));
        
        distEuc = distEuc*ratio;
    end

% else cacluate M distance by consider the cluster information
else
    % calculate the Hamming distance among viruses
    % ---------------------------------------------------------------------
    nSerumCluster = numel(serumCluster);
    
    if nCluster ~= nSerumCluster
        error('The number of virus cluster and serum cluster are differrnt!');
    end
    
    
    distHamming = zeros(rMatrix);
    
    % identify the cluster information
    clusterIndex = zeros(rMatrix,1);
    
    for i = 1: nCluster
        clusterIndex(virusCluster{i}) = i;
    end

    % calculate hamming distances
    for i = 1: rMatrix-1
        for j = i+1: rMatrix
            
            % find the set of serum used to calculate distance
            
            serumVec = [serumCluster{clusterIndex(i)} ;serumCluster{clusterIndex(j)}];
            serumVec = unique(serumVec);
            
            distHamming(i,j) = mean(abs(similarMatrix(i,serumVec)- similarMatrix(j,serumVec)));
            distHamming(j,i) = distHamming(i,j);
        end
    end
   
    % calculate the scaled Eclidean distance
    % ---------------------------------------------------------------------
    if strcmpi(distType, 'Euclidean')
        
        distEuc = zeros(rMatrix);
        
        % calculate eucliean distances
        for i = 1: rMatrix-1
            for j = i+1: rMatrix
                
                % find the set of serum used to calculate distance
                serumVec = [serumCluster{clusterIndex(i)} ;serumCluster{clusterIndex(j)}];
                serumVec = unique(serumVec);
                
                distEuc(i,j) = sqrt(mean((similarMatrix(i,serumVec)- similarMatrix(j,serumVec)).^2));
                distEuc(j,i) = distEuc(i,j);
            end
        end 
        
    end
    % --------------------------------------------------------------------- 
end

end
