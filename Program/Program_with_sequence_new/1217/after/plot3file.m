clear;clc;
fileName = 'H3N2-68-03-HI.tab';
thresh = 10;
clusterVec = {'HK68', 'EN72', 'VI75', 'TX77','BK79','SI87','BE89','BE92','WU95','SY97','FU02'};
subTable = 'H3N2-68-03sub.tab';
[virusCluster serumCluster] = tableCluster(fileName, clusterVec, subTable);
save('virusCluster.mat','virusCluster');

load recMatrix;

% A1=pdist(recMatrix,'cosine');
% B1=squareform (A1);
% distMatrix=B1;

distance = 'Hamming';
[distHamming, distEuc] = computeDist(recMatrix, distance, virusCluster, serumCluster);
if strcmpi(distance, 'Euclidean')
    distMatrix = distEuc;
 else
    distMatrix = distHamming;
end
colorVec = ['b', 'c','g','r', 'm','k','b', 'c','g','r', 'm','k','b', 'c','g','r', 'm','k','b', 'c','g','r', 'm','k'];
shapeVec = ['o', 'd', 'x', '+', '*', 's',  'v', '^', '<', '>', 'p', 'h','o', 'd', 'x', '+', '*', 's',  'v', '^', '<', '>', 'p', 'h'];
clusterName=clusterVec;
nCluster = numel(clusterName);

distCol = getUpTriangle(distMatrix);
distHammingCol = getUpTriangle(distHamming);
rMatrix = size(distMatrix,1);
mdsCorr = mdscale(distCol,3);

 mdsCorr = mdscale(distCol,3);
     
     
     % calculate ratio
     realHam = mean(distHammingCol);
    
     distEstAll = 0;
     % estimated hamming
     for i = 1: rMatrix-1
         for j = i+1: rMatrix
             distEstAll = distEstAll + sqrt(sum((mdsCorr(i,:)- mdsCorr(j,:)).^2));
         end
     end
     distEst = distEstAll/numel(distCol);
     
     ratio = realHam/distEst;
     
     mdsCorr = mdsCorr*ratio;

% plot the 3d graph
for i = 1: nCluster
   index = virusCluster{i};
   scatter3(mdsCorr(index,1), mdsCorr(index,2), mdsCorr(index,3), colorVec(i), shapeVec(i));
   hold on
end
hold off
