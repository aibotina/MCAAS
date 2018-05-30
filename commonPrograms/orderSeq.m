function orderSeq(protienFile, DNAFile, outFile)
%% ORDERSEQ  orders the sequences in DNAFile according to the order of 
%  sequences in proteinFile.
%
%    ORDERSEQ with less than two parameters  give this help.
%
%    ORDERSEQ(protienFile, DNAFile, outFile) orders the sequences and
%             output to outFile
%    
% Author: Jialiang Yang, CVM, MSU, jyang@cvm.msstate.edu
% Revision: 12/1/2011

%% Input checking
if nargin < 2
    help orderSeq
    return
end

%% order process
% Read the sequence names and sequences from ProtienFile and DNAFile
[seqProtein, nameProtein] = readFasta(protienFile);
[seqDNA, nameDNA] = readFasta(DNAFile);

% reorder DNA sequences according to order of protein sequences.
nSeq = numel(nameProtein);
nseqDNA = numel(nameDNA);

index = zeros(nSeq,1);   % index stores the positions of 

% find the index of nameDNA in nameProtein
for i = 1: nSeq
    for j = 1: nseqDNA
        if strcmp(nameProtein{i},nameDNA{j})
            index(i) = j;
            break;
        end
    end
end

% reorder DNA sequences seqDNA and write to a fasta file
% Open file fileName for writing
fid = fopen(outFile,'w+');

if fid == (-1)
    error('Can not open file %s for writing.\n Check wirte permission!',outFile)
end


if iscell(seqDNA)
    for i = 1: nSeq
        fprintf(fid, '>%s\n', nameProtein{i});
        fprintf(fid, '%s\n', seqDNA{index(i)});
    end
else
    for i = 1: nSeq
        fprintf(fid, '>%s\n', nameProtein{i});
        fprintf(fid, '%s\n', seqDNA(index(i),:));
    end
    
end


fclose(fid);

end