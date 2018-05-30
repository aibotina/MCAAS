function orderSeqByYear(alnFile, outFile)
%ORDERSEQBYYEAR orders the virus sequences in alnFile by year and output
%the ordered sequences to outFile in fasta format
% 
% USAGE:
%       orderSeqByYear(alnFile, outFile)
% INPUT:
%       alnFile: the alignment file in fasta format
%       outFile: the fileName stores the ordered alignment file
% OUTPUT:
%       None
%
%  Revision Date : 29th Dec, 2011
%  Author: Jialiang Yang  CVM, MSU, jyang@cvm.msstate.edu

%% read and reorder sequences
[virusSeq, virusName] = readAln(alnFile);

virusName
virusYear = getYear(virusName);


[~, ordIndex] = sort(virusYear);

orderSeq = virusSeq(ordIndex, :);
orderName = virusName(ordIndex);

writeFasta(orderSeq, orderName, outFile)
end

