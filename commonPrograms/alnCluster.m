function virusCluster = alnCluster(alnFile, clusterNameVec, subAlnFile)
%% alnCluster selects the sub alignment in the specified clusters, and  
%  returns the indices of the viruses in each cluster.
%
%  USAGE:
%        alnCluster   with less than 2 parameters gives this help 
%        virusCluster = alnCluster(alnFile, clusterName, subAln)
%  INPUT: 
%        alnFile: alignment file name in fasta format
%        clusterNameVec: a cell array of strings contains the cluster names
%        subAln: the sub alignment file name containing viruses in the clusters
%  Output:
%        virusCluster: the clusters of viruses
%
%  Revision Date : 6th Jan, 2012
%  Author: Jialiang Yang  CVM, MSU, jyang@cvm.msstate.edu

%% Input checking
if nargin < 2
    help alnCluster
    return
elseif nargin <3
    subAlnFile = 'subCluster.aln';
else
    % check alnFile
    if ~ischar(alnFile)
        error('alnFile must be a character array!')
    end
    
    % check subAln
    if ~ischar(subAlnFile)
        error('subAln must be a character array!')
    end
    
    % check clusterNameVec
    if ~iscellstr(clusterNameVec)
        error('clusterNameVec must be a cell array of strings containing virus cluster names!')
    end
end

constCluster;

%% construct sub alignment Find the indices of virus for each cluster
[virusSeq virusName] = readAln(alnFile);

nVirus = numel(virusName);
nCluster = numel(clusterNameVec);
clusterNum = zeros(nCluster,1);

% Check the validality of clusterNameVec
for i = 1: nCluster
    indVec = strcmpi(clusterNameVec{i}, CLUSTERNAME);
    
    if ~any(indVec)
        error('cluster ''%s'' is not defined in CLUSTERNAME', clusterNameVec{i});
    else
        clusterNum(i) = find(indVec == 1);
    end
end

% retrieve virus in the clusters 
% -------------------------------------------------------------------------
vCluster = cell(nCluster,1);

% initialize
for i = 1: nCluster
    vCluster{i} = [];
end

rowSelect = [];

for i = 1: nVirus
    for j = 1: nCluster
        if any(strcmpi(virusName{i}, CLUSTER{clusterNum(j)}))
            rowSelect = [rowSelect; i];
            
            vCluster{j} = [vCluster{j} i];
        end
    end
end


% calculate virusCluster containg the relative positions in vCluster
% -------------------------------------------------------------------------
virusCluster = cell(nCluster,1);

% find the relative position
for i = 1: nCluster
    [commonVirus irowSelect] = intersect(rowSelect, vCluster{i});
    virusCluster{i} = irowSelect;
end
% -------------------------------------------------------------------------

% write to new sub aln
% -------------------------------------------------------------------------
subSeq = virusSeq(rowSelect,:);
subName = virusName(rowSelect);

writeFasta(subSeq, subName, subAlnFile);

end