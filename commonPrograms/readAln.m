function [seq, seqName] = readAln(fileName)
%% readAln reads sequence alignment file in fasta format 
%
%  USAGE:
%        readAln  gives this help on the readAln function.
%        [seq, seqName]= readAln(fileName)
%  INPUT: 
%        fileName: a sequence alignment file in fasta format 
%  OUTPUT:
%            seq: the vector or cell array of the sequences
%        seqName: the cell array of sequence names
%
%  Revision Date : 29th Dec, 2011
%  Author: Jialiang Yang  CVM, MSU, jyang@cvm.msstate.edu

%% Input checking
if nargin == 0
    help readAln
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
seq = [];
% ------------------------------------------------------------------------

% Add sequence names and sequences to seqName and seq
% ------------------------------------------------------------------------
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
% -------------------------------------------------------------------------

end