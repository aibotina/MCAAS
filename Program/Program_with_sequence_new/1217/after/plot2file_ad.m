clear;clc;
fileName = 'H3N2-68-03-HI.tab';
thresh = 10;
clusterVec = {'HK68', 'EN72', 'VI75', 'TX77','BK79','SI87','BE89','BE92','WU95','SY97','FU02'};
subTable = 'H3N2-68-03sub.tab';
[virusCluster serumCluster] = tableCluster(fileName, clusterVec, subTable);
load recMatrix.mat;

%  distance = 'Euclidean';
 distance = 'Hamming';

[distHamming, distEuc] = computeDist(recMatrix, distance, virusCluster, serumCluster);
 if strcmpi(distance, 'Euclidean')
    distMatrix = distEuc;
 else
    distMatrix = distHamming;
end

colorVec = ['g', 'c','b','r', 'm','y','g', 'c','b','r', 'm','k','b', 'c','g','r', 'm','k','b', 'c','g','r', 'm','k'];
shapeVec = ['o', 'd', 'x', '+', '*', 's',  'v', '^', '<', '>', 'p', 'h','o', 'd', 'x', '+', '*', 's',  'v', '^', '<', '>', 'p', 'h'];
clusterName=clusterVec;
nCluster = numel(clusterName);

h = figure;
box on;
grid on; 
% [C,IA,IC] = unique(Seq_pdis,'rows','stable');
% [CC,CIA,CIC] = unique(C','rows','stable');
distCol = getUpTriangle(distMatrix);
distHammingCol = getUpTriangle(distHamming);
rMatrix = size(distMatrix,1);
mdsCorr = mdscale(distCol,2);

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
 
 
%  % plot the 2d graph
 hold on
 minX = floor(min(mdsCorr(:,1)));
 maxX = ceil(max(mdsCorr(:,1)))+1;
 minY = floor(min(mdsCorr(:,2)));
 maxY = ceil(max(mdsCorr(:,2)));
    
% axis([minX, maxX,minY,maxY]);
% set(gca,'Xtick',minX:maxX);
% set(gca, 'Ytick', minY:maxY);

% axis([-0.1, 0.2,-0.1,0.15]);
% set(gca,'Xtick',-0.1:0.2);
% set(gca, 'Ytick', -0.1:0.15);

 for i = 1: nCluster
            index = virusCluster{i};
            scatter(mdsCorr(index,1), mdsCorr(index,2), colorVec(i), shapeVec(i));    
            text(mean(mdsCorr(index,1)),mean(mdsCorr(index,2)),clusterVec(i));
 end
 
 hold off
 
  
 













% distCol = getUpTriangle(distMatrix);
% mdsCorr = mdscale(distCol, 2);


% load distGenetic;
% distMatrix=distGenetic;
% m=11;
% distMatrix=distMatrix(1:m,1:m);
% distCol1 = getUpTriangle(distMatrix);
% mdsCorr1 = mdscale(distCol1, 2);
