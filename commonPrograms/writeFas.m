function writeFas(seq, seqName, outFile)
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
    fprintf(fid,'%s\n', seq{i});
end

fclose(fid);

end