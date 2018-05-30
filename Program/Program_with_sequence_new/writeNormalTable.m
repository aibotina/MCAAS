function writeNormalTable(fileName, dataHI, virusName, serumName)
%% writeNormalTable writes the HI data into filename with prescribed format 
%
%    Usage:
%          writeNormalTable    gives this help.
%          writeNormalTable(fileName, dataHI, virusName, serumName)
%    Input:
%          fileName: the output file name.
%            dataHI: the HI table
%         virusName: a vector stores all the virus names
%         serumName: a vector stores all the serum names
%
%  Revision Date : 2nd Jan, 2011
%  Author: Jialiang Yang  CVM, MSU, jyang@cvm.msstate.edu

%% Input checking
% Parameters fileName, dataHI, virusName, serumName should be specified
if nargin == 0
    help writeNormalTable
    return
elseif nargin < 4
    error('Parameters: fileName dataHI virusName serumName')   
end

if ~ischar(fileName)
    error('FileName must be a character array!')
end

% Check the validality of the input data 
[dataRow,dataCol] = size(dataHI);
nVirus = length(virusName);
nSerum = length(serumName);
nSerum = dataCol;

% Check viruses
if dataRow ~= nVirus
    error('The number of viruses is differnt in dataHI and virusName!')
end

% Check serum
% if dataCol ~= nSerum
%     error('The number of serums is differnt in dataHI and serumName!')
% end
    
%% Open file fileName for writing

% Open file fileName for writing
fid = fopen(fileName,'w+');

if fid == (-1)
    error('Can not open file %s for writing.\n Check wirte permission!',...
         fileName)
end

% Write data to fileName

% Virus number and Serum number line
fprintf(fid, '%d\t%d\n', nVirus, nSerum);

% Serum line
fprintf(fid,'ID\t');

for i = 1: nSerum-1
    fprintf(fid,'%s\t', serumName{i});
end
fprintf(fid,'%s\n', serumName{nSerum});

% Reaction table
if isnumeric(dataHI)
    
    for i = 1: nVirus
        
        % Virus name
        fprintf(fid,'%s \t', virusName{i});
        
        % table
        for j = 1: nSerum-1
            fprintf(fid,'%d \t', dataHI(i,j));
        end
        
        fprintf(fid,'%d \n', dataHI(i,nSerum));
    end
else
    for i = 1: nVirus
        % Virus name
        fprintf(fid,'%s \t', virusName{i});
        
        % table
        for j = 1: nSerum-1
            fprintf(fid,'%s \t', dataHI{i,j});
        end
        
        fprintf(fid,'%s \n', dataHI{i,nSerum});
    end
end

fclose(fid);

end
