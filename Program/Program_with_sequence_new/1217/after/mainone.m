

fileName = 'H3N2-68-03-HI.tab';

thresh = 10;
norm=1;
[recMatrix,virusName] = matrixCompletion2(fileName, thresh);
save('recMatrix.mat','recMatrix');