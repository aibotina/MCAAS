function [dChar freq] = dominantChar(seqCol)
%% dominantChar finds the dominant char and its frequency in seqCol
%   USAGE:
%        [dChar freq] = dominantChar(seqCol)
%   INPUT:
%        [Compulsary]
%             seqCol:  a character vector
%   OUTPUT:
%       dChar: dominant character in seqCol
%       freq: the frequency of the dominant char
%              
%% identify dominant char and frequency
[char, pos, charPos]= unique(seqCol);

nPos = numel(pos);

charNum = 1;

currentFreq = sum(charPos == 1);

for i = 2: nPos
    charFreq = sum(charPos == i);
    
    if currentFreq < charFreq
        charNum = i;
        currentFreq = charFreq;
    end
end

dChar = char(charNum);
freq = currentFreq/numel(seqCol);

end

