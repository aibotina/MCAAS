function writeTable(fileName, dataHI, virusName, serumName, reference)
%% writeTable writes the HI data into filename with prescribed format 
%
%    Usage:
%          writeTable    gives this help on the writeTable function.
%          writeTable(fileName, dataHI, virusName, serumName, reference)
%    Input:
%          fileName: the output file name.
%            dataHI: the HI table
%         virusName: a vector stores all the virus names
%         serumName: a vector stores all the serum names
%         reference: reference HI value for each serum 
%
%  Revision Date : 29th Dec, 2011
%  Author: Jialiang Yang  CVM, MSU, jyang@cvm.msstate.edu

%% Input checking
% Parameters fileName, dataHI, virusName, serumName should be specified
if nargin == 0
    help writeTable
    return
elseif nargin < 5
    error('Parameters: fileName dataHI virusName serumName reference')   
end

if ~ischar(fileName)
    error('FileName must be a character array!')
end

%% Write the HI table into file fileName

% Check the validality of the input data 
[dataRow,dataCol] = size(dataHI);
nVirus = numel(virusName);
nSerum = numel(serumName);

% Check viruses
if dataRow ~= nVirus
    error('The number of viruses is differnt in dataHI and virusName!')
end

% Check serum
if dataCol ~= nSerum
    error('The number of serums is differnt in dataHI and serumName!')
end

% Check reference
if numel(reference) ~= dataCol
    error('The number of referce is differnt from the col of dataHI!')
end
    
%% Open file fileName for writing

% Check if the file fileName exists
% if exist(fileName,'file')
%     error('File %s already exists', fileName)
% end

% Open file fileName for writing
%--------------------------------------------------------------------------
fid = fopen(fileName,'w+');

if fid == (-1)
    error('Can not open file %s for writing.\n Check wirte permission!',...
         fileName)
end

% Write data to fileName
%--------------------------------------------------------------------------
% Virus number and Serum number line
fprintf(fid, '%d \t %d \n', nVirus, nSerum);

% Serum line
fprintf(fid,'ID \t');

for i = 1: nSerum-1
    fprintf(fid,'%s \t', serumName{i});
end

fprintf(fid,'%s \n', serumName{nSerum});

% Reference line
fprintf(fid,'Reference \t');

for i = 1: nSerum-1
    fprintf(fid,'%s \t', reference{i});
end

fprintf(fid,'%s \n', reference{nSerum});

% Reaction table
for i = 1: nVirus
    % Virus name
    fprintf(fid,'%s \t', virusName{i});
    
    % table
    for j = 1: nSerum-1
        fprintf(fid,'%s \t', dataHI{i,j});
    end
    
    fprintf(fid,'%s \n', dataHI{i,nSerum});
end

%--------------------------------------------------------------------------

fclose(fid);

end