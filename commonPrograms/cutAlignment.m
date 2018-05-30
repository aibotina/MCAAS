function cutAlignment(fileName, posVec, outFile)
%% cutAlignment cuts the alignment (fasta format) from startPos to endPos 
%  and write the new alignment into outFile (fasta format)
%      
% Author: Jialiang Yang, CVM, MSU, jyang@cvm.msstate.edu
% Revision: 12/5/2011

%% Input checking
% check fileName
if ~ischar(fileName)
    error('FileName must be a character array!')
end

% read sequence and virus name
[seq, seqName] = readAln(fileName);

seqNew = seq(:,posVec); 

writeFasta(seqNew, seqName, outFile);
end