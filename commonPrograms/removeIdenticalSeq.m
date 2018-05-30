function removeIdenticalSeq(alnFile, newAlnFile, siteVec)
%% REMOVEIdenticalSeq removes identical sequences in an alignment file
%
%  Usage:
%        removeIdenticalSeq(alnFile, newAlnFile)
%
%  Input: 
%       alnFile: full alignment file in fasta format
%      newAlnFile: the output file with no identical seqs
%       siteVec: the sites to check if the sequences of the viruses
%       are identical at the sites 
%
%  Revision Date : May 4, 2012
%  Author: Jialiang Yang  CVM, MSU, jyang@cvm.msstate.edu

%% Input checking

if nargin < 2
    help removeIdenticalSeq
    return
elseif nargin < 3
    siteVec = [];
else
    if ~ischar(alnFile)
        error('AlnFile must be a character array!')
    end
    
    if ~ischar(newAlnFile)
        error('newAlnFile must be a character array!')
    end
end

%% Retrieve sub table
disp('Read full alignment start');
[virusSeq, virusName] = readAln(alnFile);
disp('Read full alignment End');

if ~isempty(siteVec)
    virusSeqTrim = virusSeq(:, siteVec);
else
    virusSeqTrim = virusSeq;
end
    
nVirus = numel(virusName);
alnNoRepeat = [];
selectInd = [];

% Check for repeats
for i = 1: nVirus
    currentSeq = virusSeqTrim(i,:);
    
    flag = 0;
    
    for j = 1: size(alnNoRepeat,1)
        compareSeq = alnNoRepeat(j,:);
        
        if strcmpi(currentSeq, compareSeq)
            flag = 1;
            break;
        end
    end
    
    if flag == 0
        alnNoRepeat = [alnNoRepeat; currentSeq];
        
        selectInd = [selectInd, i];
    end
end

alnNoRepeat = virusSeq(selectInd,:);
nameNoRepeat = virusName(selectInd);

writeFasta( alnNoRepeat, nameNoRepeat, newAlnFile);

end