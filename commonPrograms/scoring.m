function seq = scoring(seq1,seq2,score)
%% scoring generates a score vector of seq1 and seq2 by different scoring
%  schemes. Note: should run const.m first
%
%    USAGE:
%      seq = scoring(seq1,seq2,score)
%
%    INPUT: 
%       seq1,seq2 : vector or string with the same length
%       score: scoring schemes
%              1: 0 same; 1 different;
%              2: use PIMA scoring matrix, gap score 1
%              3: use PAM250 scoring matrix, gap score 3
%              4: use BLOSUM62 scoring matrix, gap score 4
%    OUTPUT:
%       seq: vector or string after comparision
%
% Author: Jialiang Yang, CVM, MSU, jyang@cvm.msstate.edu

%% Input checking
% Set default values: score = 3; 

if nargin < 2
    help scoring
    return
elseif nargin == 2
    score = 3;
end

const
global AAINDEX;
global PIMA;
global PAM250;
global BLOSUM62;

%% perform different scoring scores

cSeq = numel(seq1);      % cSeq denotes the number of characters in seq1
seq = zeros(1,cSeq);     % initialize seq, the sequence score vector

switch score
    case 1  % score 1: 0 same; 1 different; gap 0
        for i = 1:cSeq
            if seq1(i) == seq2(i) || seq1(i) == '-' || seq2(i) == '-'
                seq(i) = 0;
            else
                seq(i) = 1;
            end
        end       
        
    case 2 % score 2: use PIMA scoring matrix; gap 0
        for i = 1:cSeq
            if seq1(i) == seq2(i) || seq1(i) == '-' || seq2(i) == '-' || ~any(AAINDEX == seq1(i)) || ~any(AAINDEX == seq2(i))
                seq(i) = 0;
            else
                seq(i) = PIMA(AAINDEX == seq1(i), AAINDEX == seq2(i));
            end
        end
        
    case 3 % score 3: use PAM250 scoring matrix; gap 3
       for i = 1:cSeq
            if seq1(i) == '-' || seq2(i) == '-'|| ~any(AAINDEX == seq1(i)) || ~any(AAINDEX == seq2(i))
                seq(i) = 3;
            else
                seq(i) = PAM250(AAINDEX == seq1(i), AAINDEX == seq2(i));
            end
       end
        
    case 4 % score 4: use BLOSUM62 scoring matrix; gap 4
       for i = 1:cSeq
            if seq1(i) == '-' || seq2(i) == '-'|| ~any(AAINDEX == seq1(i)) || ~any(AAINDEX == seq2(i))
                seq(i) = 4;
            else
                seq(i) = BLOSUM62(AAINDEX == seq1(i), AAINDEX == seq2(i));
            end
       end
end

end
