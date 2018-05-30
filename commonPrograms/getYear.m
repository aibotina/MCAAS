function year = getYear(virusName)
%% getYear calculates the year of the virus virusName: virusName must be in 
%       the format VI/7/87
%
% USAGE:
%       year = getYear(virusName) returns the year of that virus
% INPUT:
%      virusName: a string or cell array stores names of the viruses, must be in format VI/7/87
% OUTPUT:
%      year: a number or vector stores the year of the virus
%
%  Revision Date : 29th Dec, 2011
%  Author: Jialiang Yang  CVM, MSU, jyang@cvm.msstate.edu

%% identify the year of the viruses 

if ~iscell(virusName)
    %find the position of /
    pos = strfind(virusName, '/');
    
    % year located from the second / to the end of the virusName
    year = str2double(virusName(pos(2)+1: numel(virusName)));
    
    if year < 100
        if year > 50
            year = 1900+year;
        else
            year = 2000+year;
        end
    end
else
    nVirus = numel(virusName);
    year = zeros(nVirus,1);
    
    for i = 1: nVirus
        year(i) = getYear(virusName{i});
    end
end