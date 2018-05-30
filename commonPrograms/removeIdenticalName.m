function removeIdenticalName(alnFile, outFile )
%REMOVESEQOVERLAP removes the overlaped seqs in alnFile
% 
% USAGE:
%       removeSeqOverlap(alnFile, outFile )
% INPUT:
%       alnFile: the alignment file in fasta format
%       outFile: the fileName stores noOverlapped alignment file
% OUTPUT:
%       None
%
%  Revision Date : 29th Dec, 2011
%  Author: Jialiang Yang  CVM, MSU, jyang@cvm.msstate.edu

%% read and retrieve sequences
[virusSeq, virusName] = readAln(alnFile);

[~, uniqueInd] = unique(virusName);

subAln = virusSeq(uniqueInd,:);
subName = virusName(uniqueInd);

writeFasta(subAln, subName, outFile);

end

