function flag = vecInMatrix(vec, matrix)
%% judge if the vector vec is in a row of the matrix

[iRow, iCol] = size(matrix);
flag = false;

for i = 1: iRow
    currentVec = matrix(i,:);
    compVec = (vec == currentVec);
    
    if sum(compVec) == iCol
        flag = true;
        break;
    end
end

