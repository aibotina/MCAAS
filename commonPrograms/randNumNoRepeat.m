function [randVec] = randNumNoRepeat(minValue,maxValue, num)
%% RANDNUM generates num random non repeated values between 1 and maxValue
%           
%    RANDNUM gives help on RANDNUM function.
%
%    [randVec] = RANDNUM(minValue, maxValue) generatates 1 random value between minValue 
%    and maxValue
%
%    [randVec]= RANDNUM(minValue, maxvalue, num) generatates num non repeated random 
%    values between minValue and maxValue
%
% Author: Jialiang Yang, CVM, MSU, jyang@cvm.msstate.edu

%% Input checking
switch nargin
    case 1
        help randNumNoRepeat
        return
    case 2
        num = 1;
end

% Initialization
randVec = zeros(num,1);

% generate num non repeated random numbers between [1, maxValue]
k = 1;

span = maxValue -minValue+1;

while randVec(num) == 0
    % generate a random number 
    randNum = ceil(rand(1)*span);
    
    % check if it is already in randVec
    if ~any((randVec-randNum) == 0) 
        randVec(k) = randNum;
        k = k+1;
    end
end

randVec = sort(randVec+minValue-1);

end