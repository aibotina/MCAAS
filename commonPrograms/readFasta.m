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