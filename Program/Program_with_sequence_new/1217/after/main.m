clear;
clc;

fileName = 'H3N2-68-03-HI.tab';

thresh = 10;

clusterVec = {'HK68', 'EN72', 'VI75', 'TX77','BK79','SI87','BE89','BE92','WU95','SY97','FU02'};

subTable = 'H3N2-68-03sub.tab';

[virusCluster serumCluster] = tableCluster(fileName, clusterVec, subTable);



%[virusCluster clusterName] = yearTable(fileName);

%[distMatrix, mdsCorr] = antigenMapClusterNoTypeIII(fileName, lowReact ,{}, {}, {}, 3, 'Euclidean',6);
[distMatrix, mdsCorr,mat] = antigenMapCluster(fileName,thresh,virusCluster,serumCluster,clusterVec,3,'Euclidean',6);
distance='Euclidean';

% [distHamming, distEuc] = computeDist(recMatrix, 'Euclidean', virusCluster, serumCluster);
% [ distCol ] = getUpTriangle(distMatrix)

























