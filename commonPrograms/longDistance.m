function longDistanceMatrix = longDistance(shortDistanceFile, distanceFile, longDistanceFile)
%% longDistance retrieves long distances from distanceFile, and set all short distances to be -1
%
%  USAGE:
%        longDistanceMatrix = longDistance(shortDistanceFile, distanceFile, longDistanceFile)
%  INPUT: 
%        shortDistanceFile: a distance matrix file stores reliable short
%                distances and keep the not reliable long distance to be -1.
%        distanceFile: the distance matrix from temporal model.
%      longDistanceFile: a file contains the long distance Matrix
%  OUTPUT:
%        longDistanceMatrix: the long distance matrix
%
%  Revision Date : 8th July, 2012
%  Author: Jialiang Yang  CVM, MSU, jyang@cvm.msstate.edu

%% Input checking
if nargin < 3
    help longDistance
else
    % check fileName
    if ~ischar(shortDistanceFile)
        error('shortDistanceFile must be a character array!')
    end
     
    if ~ischar(distanceFile)
        error('distanceFile must be a character array!')
    end
    
    if ~ischar(longDistanceFile)
        error('longDistanceFile must be a character array!')
    end

end

%% load distance files and merge them
shortDistMatrix = load(shortDistanceFile);
distanceMatrix = load(distanceFile);

nShort = size(shortDistMatrix, 1);
nDistance = size(distanceMatrix, 1);

if nShort ~= nDistance
    error('The size of short distance matrix and long distance matrix not equal!');
end

longDistanceMatrix = zeros(nShort, nShort);

for i = 1: nShort
    for j = 1: nShort
        if shortDistMatrix(i,j) == -1
            longDistanceMatrix(i,j) = distanceMatrix(i,j);
        else
            longDistanceMatrix(i,j) = -1;
        end
    end
end

save(longDistanceFile, 'longDistanceMatrix',  '-ascii', '-double', '-tabs');

end
