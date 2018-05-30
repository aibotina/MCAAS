function orderSeqByTable(seqFile, tableFile, alnSeqFile)
%% SEQVYTABLE  orders the sequences in seqFile according to the order of 
%  viruses in tableFile.
%   Usage:
%         seqByTable(seqFile, tableFile, alnSeqFile)
%   Input:
%         seqFile: the sequence file in fasta format
%        tableFile: the HI table file
%        alnSeqFile: the new sequence file in fasta format
%    
% Author: Jialiang Yang, CVM, MSU, jyang@cvm.msstate.edu
% Revision: 22/3/2012

%% Input checking
if nargin < 3
    help setByTable
    return
end

%% order process
% Read the sequence names and sequences from seqFile
[seq, seqName] = readFasta(seqFile);
[dataHI, virusName, serumName, reference] = readTable(tableFile);

% reorder DNA sequences according to order of protein sequences.
nVirus = numel(virusName);

index = zeros(nVirus,1);

for i = 1: nVirus
    index(i) = find(strcmpi(virusName{i}, seqName) == 1);
end

seq = seq(index,:);
seqName = seqName(index);

writeFasta(seq, seqName, alnSeqFile);

end

% =========================================================================
function [seq, seqName] = readFasta(fileName)
%% readFasta reads sequence file in fasta format 
%
%    READFASTA  gives this help on the READFASTA function.
%
%    [seq, seqName]= READFASTA(FILENAME) reads taxa name to seqName and 
%    alignment to seq.
%    
% Author: Jialiang Yang, CVM, MSU, jyang@cvm.msstate.edu

%% Input checking
if nargin == 0
    help readFasta
    return
else
    % check fileName
    if ~ischar(fileName)
        error('FileName must be a character array!')
    end
end

%% Read the sequence names and sequences from file fileName
% check if the file fileName exists
if (exist(fileName,'file') || exist(fullfile(pwd,fileName),'file'))
    % ftext reads lines from the file
    ftext = strtrim(textread(fileName,'%s','delimiter','\n'));
else
    error('File %s does not exist!', fileName)
end
    
    
% The number of sequences is identified by character '>'
nameLines = strncmp(ftext,'>',1);

if ~any(nameLines)
    error('File %s is not in fasta format!', fileName)
end

% numSeq denotes the number of sequences
numSeq = sum(nameLines);

% seqStarts denotes the starting lines of new sequences
seqStarts = [find(nameLines); size(ftext,1)+1]; % size(ftext,1)+1 the end

seqName = cell(numSeq,1);
seq = [];

% Add sequence names and sequences to seqName and seq
try 
    for i = 1: numSeq
       
        % Add sequence names to seqName
        seqName{i} = ftext{seqStarts(i)}(2:end);
        
        % Add seqence to seq
        tempSeq = [];
        
        for j = seqStarts(i)+1 : seqStarts(i+1)-1
            tempSeq = [tempSeq, ftext{j}]; 
        end;
        
        seq = [seq; tempSeq];
    end
catch exception
     error('Incorrect fasta data format in %s!',fileName)
end

end

%==========================================================================
function [dataHI, virusName, serumName, reference] = readTable(fileName)
%% readTable reads HI or MN assay table from file fileName 
%
%  Usage:
%        readTable  gives this help on the readAln function.
%        [lowReactor,dataHI, virusName, serumName, reference] =
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
    fText = textread(fileName,'%s','delimiter','\n','bufsize',100000);
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
end

%==========================================================================
function writeFasta(seq, seqName, outFile)
%% sriteFasta write seq and seqName to outFile in fasta format 
%
%    writeFasta(seq, seqName, outFile)
%    
% Author: Jialiang Yang, CVM, MSU, jyang@cvm.msstate.edu
% Revision: 12/5/2011


%% open file and write
fid = fopen(outFile, 'w+');

if fid == (-1)
    error('Can not open file %s for writing.\n Check wirte permission!',outFile)
end

nSeq = numel(seqName);

for i = 1: nSeq
    fprintf(fid,'>%s\n', seqName{i});
    fprintf(fid,'%s\n', seq(i,:));
end

fclose(fid);

end

%==========================================================================