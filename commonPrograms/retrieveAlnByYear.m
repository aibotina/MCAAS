function retrieveAlnByYear(alnFile, yearVec, subAlnFile)
%% RETRIEVEALNBYYEAR retrives sub alignment by years specified in yearVec
%
%  Usage:
%        retrieveAlnByYear(alnFile, yearVec, subAln)
%
%  Input: 
%       alnFile: full alignment file in fasta format
%        yearVec: the years of the alignment to retrieve
%       subAln: the output sub alignment file in fasta format
%
%  Revision Date : May 4, 2012
%  Author: Jialiang Yang  CVM, MSU, jyang@cvm.msstate.edu

%% Input checking

if nargin == 0
    help retrieveAlnByYear
    return
else
    if ~ischar(alnFile)
        error('AlnFile must be a character array!')
    end
    
    if ~ischar(subAlnFile)
        error('subAlnFile must be a character array!')
    end
    
    if ~isnumeric(yearVec)
        error('yearVec must be a vector!');
    end
end

%% Retrieve sub table

disp('Read full alignment start');
[virusSeq, virusName] = readAln(alnFile);
disp('Read full alignment End');

nVirus = numel(virusName);
subVirusInd = [];

% separate the viruses by year
for i = 1: nVirus
    currentVirusName = virusName{i};
    
    virusYear = getYear(currentVirusName);
    
    if any(virusYear == yearVec)
        subVirusInd = [subVirusInd i];
    end
end

% write to subAln
subAln = virusSeq(subVirusInd, :);
subName = virusName(subVirusInd);

writeFasta(subAln, subName, subAlnFile);

% end