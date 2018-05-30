function [dataHI, virusName, serumName, reference] = readTable(fileName)
%% readTable reads HI or MN assay table from file fileName 
%
%  Usage:
%        readTable  gives this help on the readAln function.
%        [dataHI, virusName, serumName, reference] =
%        readTable(fileName)
%  Input: 
%       fileName: a tab deliminated file with the following format
%                    2 \tab 3
%                    ID     \tab  serum1 \tab serum2 \tab serum
%                 reference \tab     0   \tab    0   \tab   0
%                  virus1   \tab  data11 \tab data12 \tab data13
%                  virus2   \tab  data21 \tab data22 \tab data23
%  Output:
%            dataHI: the HI table
%         virusName: a vector stores all the virus names
%         serumName: a vector stores all the serum names
%         reference: reference HI value for each serum 
%
%  Revision Date : 29th Dec, 2011
%  Author: Jialiang Yang  CVM, MSU, jyang@cvm.msstate.edu

%% Input checking

if nargin == 0
    help readTable
    return
else
    if ~ischar(fileName)
        error('FileName must be a character array!')
    end
    
    % Set default values
    del = '\s+';
end


%% Read data from file fileName

% Check if the file fileName exists
if (exist(fileName,'file') || exist(fullfile(pwd,fileName),'file'))
    % ftext reads lines from the file
    fText = textread(fileName,'%s','delimiter','\n','bufsize',1000000);
else
    error('File %s does not exist!', fileName)
end

% Remove spaces at the begining and end and remove blank lines if exist
fTrim = {};

for i = 1: length(fText)
    lineTrim = strtrim(fText{i});   % Remove spaces at the begin and end
    
    % Remove blank lines
    if ~isempty(lineTrim)
        fTrim = [fTrim; lineTrim];
    end
end    
        
nVirus = length(fTrim)-3;     % nVirus denotes the number of viruses

% reads the data
% -------------------------------------------------------------------------
% The second line is the serum line
serumName = regexp(fTrim{2},del,'split'); % Split serum line
serumName(1) = [];             % Remove 'ID'
nSerum = length(serumName);

% The third line is the reference line
reference = regexp(fTrim{3},del,'split');
reference(1) = [];

% Read virus names and reaction data
dataHI = {};
%lowReactor = 0;
virusName = {};

for i = 4: nVirus + 3
    sLine = regexp(fTrim{i},del,'split'); % Split line i
    
    virusName = [virusName; sLine{1}];    % Virus name is in the first col
    
    sLine(1) = [];              % Remove virus name
    
    try
        dataHI = [dataHI; sLine];
    catch dataException
        error('Table width is inconsistent in each row!')
    end
    
end
% -------------------------------------------------------------------------
   
end
