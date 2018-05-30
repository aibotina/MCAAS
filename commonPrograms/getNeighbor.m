function coNeighbor = getNeighbor(folderName)
%% getNeighbor generates restriction on neighborhood from folder folderName
%  
%  neighbor = getNeighbor(folderName)
%
% Author: Jialiang Yang, CVM, MSU, jyang@cvm.msstate.edu
% revision: 25/11/2011

%% Input checking
if nargin == 0
    help getNeighbor
    return
else
    if ~ischar(folderName)
        error('FolderName must be a character array!')
    end
end

%% Identify all files in the foloder.

Files = dir(strcat(folderName,'\\*.txt'));

fileNum = numel(Files);  % number of files in the folder

neighborPair = [];

% read the sites from all files into allSite
for i = 1: fileNum
    currentFile = Files(i).name;
    pos = strfind(currentFile, '_');
    currentSite = str2double(currentFile(1:pos-1));
    % get neighbor sites
    neighborFile = [folderName, '\\', currentFile];
    siteInFile = readNeighborSite(neighborFile);
    
    for j = 1: numel(siteInFile)
        if(currentSite > siteInFile(j))
            pairSite = [siteInFile(j), currentSite];
            neighborPair = [neighborPair; pairSite];
        elseif(currentSite < siteInFile(j))
            pairSite = [currentSite,siteInFile(j)];
            neighborPair = [neighborPair; pairSite];
        end
    end
end

coNeighbor = unique(neighborPair,'rows'); 

% write into file
fileName = [folderName, '.txt'];

% Open file fileName for writing
fid = fopen(fileName,'w+');

if fid == (-1)
    error('Can not open file %s for writing.\n Check wirte permission!',fileName)
end

% write coNeighbor sites into file
for i = 1: size(coNeighbor,1)
    fprintf(fid, '%d %d\n', coNeighbor(i,1), coNeighbor(i,2));
end

fclose(fid);

end