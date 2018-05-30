function orderSeqBySeq(seqFile, tableFile, alnSeqFile)
%% orderSeqBySeq orders the sequences in seqFile according to the order of 
%  viruses in a fas file.
%   Usage:
%         orderSeqBySeq(seqFile, tableFile, alnSeqFile)
%   Input:
%         seqFile: the sequence file in fasta format
%        tableFile: the seq file in fas format 
%        alnSeqFile: the new sequence file in fasta format
%    
% Author: Jialiang Yang, CVM, MSU, jyang@cvm.msstate.edu
% Revision: 22/3/2012

%% Input checking
if nargin < 3
    help orderSeqBySeq
    return
end

%% order process
% Read the sequence names and sequences from seqFile
[seq, seqName] = readFasta(seqFile);
[seq2, virusName] = readFasta(tableFile);

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

function [seq, seqName] = readFasta(fileName)
%% readFasta reads sequence file in fasta format 
%
%    readFasta  gives this help on the readFasta function.
%
%    [seq, seqName]= readFasta(fileName) reads taxa name to seqName and 
%    alignment to seq.
%    
% Author: Jialiang Yang, CVM, MSU, jyang@cvm.msstate.edu
% Revision: 11/13/2011

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
%--------------------------------------------------------------------------
nameLines = strncmp(ftext,'>',1);

if ~any(nameLines)
    error('File %s is not in fasta format!', fileName)
end

% numSeq denotes the number of sequences
numSeq = sum(nameLines);

% seqStarts denotes the starting lines of new sequences
seqStarts = [find(nameLines); size(ftext,1)+1]; % size(ftext,1)+1 the end

seqName = cell(numSeq,1);
seq = {};

% Add sequence names and sequences to seqName and seq
% ------------------------------------------------------------------------
% try 
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
% catch exception
%      error('Incorrect fasta data format in %s!',fileName)
% end
% -------------------------------------------------------------------------

end

%==========================================================================
function writeFasta(seq, seqName, outFile)
%% writeFasta write seq and seqName to outFile in fasta format 
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
    fprintf(fid,'%s\n', seq{i});
end

fclose(fid);

end
%==========================================================================