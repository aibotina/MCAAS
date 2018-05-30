function writeXml(fileName, virusName, virusColor, virusCorr)
%% writeXml writes the data into an Xml file

%% Input checking
% Parameters fileName, dataHI, virusName, serumName should be specified
if nargin == 0
    help writeXml
    return
elseif nargin < 4
    error('Parameters: fileName, virusName, virusColor, virusCorr')   
end

if ~ischar(fileName)
    error('FileName must be a character array!')
end

%% Write the data into file fileName
fid = fopen(fileName,'w+');

if fid == (-1)
    error('Can not open file %s for writing.\n Check wirte permission!',fileName)
end

nVirus = numel(virusName);

% Virus number and Serum number line
fprintf(fid,'%d\n', nVirus);
fprintf(fid, 'Name\tGroup\tX\tY\tZ\n');

for i = 1: nVirus    
    fprintf(fid, '%s\t%s\t%d\t%d\t%d\n', virusName{i}, virusColor{i}, virusCorr(i,1), virusCorr(i,2), virusCorr(i,3));
end

fclose(fid);

end