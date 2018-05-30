function [dataNorm, typeOne, typeTwo,virusName, serumName] = normalizeTable(fileName, thresh, norm)
%%  normalizeTable normalize the antigenic table in fileName
%
%    USAGE:
%        NORMALIZETABLE with less than 2 parameters give this help.
%        [dataNorm, typeTwo,virusName, serumName] = NORMALIZETABLE(fileName, thresh, norm)
%             normalizes the table and returns the position of typeII data.
%
%    Input:
%             fileName: the file contains the HI table
%             thresh: the threshold for typeII data
%             norm: normalization schemes 1-5.
%
%    Output:
%             dataNorm: normalized HI table
%             typeTwo: a matrix saves the positions of typeII data.
%             virusName: a vector stores all ordered virus names.
%             serumName: a vector stores all serum names.
%             
% Author: Jialiang Yang, CVM, MSU, jyang@cvm.msstate.edu
% Revision: 12/12/2011

%% Input checking

if nargin < 2
    help normalizeTable
    return
    
elseif nargin <3
    norm =1;  % specifies normalization scheme 
    
end

% check fileName
if ~ischar(fileName)
    error('FileName must be a character array!')
end
  
%% Read data from file fileName
[dataHI, virusName, serumName, reference] = readTable(fileName);

[nVirus, nSerum] = size(dataHI);

% tableHI stors the numerical HI table: change typeII data to threshold
%-------------------------------------------------------------------------
tableHI = zeros(nVirus, nSerum);
 typeOne=[];
typeTwo = [];

for i = 1: nVirus
     for j = 1: nSerum
         HIValue = strtrim(dataHI{i,j});
         
         % missing value
         if strcmp(HIValue,'0')
             tableHI(i,j) = 0;
             
         % lower reactor like '<20', take tableHI(i,j) = 20
         elseif strcmp(HIValue(1),'<') 
             tableHI(i,j) = str2double(HIValue(2: numel(HIValue)));
             %tableHI(i,j) = thresh;
             
             % update typeTwo matrix
             pos = [i,j];
             typeTwo = [typeTwo; pos];
             
         else
             value = str2double(HIValue);
             
             % lower reactor like 25 < 40 take tableHI(i,j) = 25 
             if value < thresh
                 tableHI(i,j) = value;
                 %tableHI(i,j) = thresh;
                 
                  % update typeTwo matrix
                 pos = [i,j];
                 typeTwo = [typeTwo; pos];
                 
             % type I value
             else
                tableHI(i,j) = value;
                 pot=[i,j];
                 typeOne=[typeOne;pot];
             end
         end
     end
end

% normalize data according to different scoring schemes
%--------------------------------------------------------------------------
dataNorm = zeros(nVirus, nSerum);  % stores the data after normalization

switch norm
    case 1
        % identify max overall and max in each column.
        maxCol = max(tableHI);
        
        logCol = log2(maxCol);  % log2(max(Hj))
        
        % normalize typeI and type II entries by
        % log2(max(Hij))-log2(max(Hj))+log2(Hij) (or log2(thresh))
        for i = 1: nVirus
            for j = 1: nSerum
                % not missing value
                if tableHI(i,j) ~=0
                    dataNorm(i,j) = logCol(j) - log2(tableHI(i,j));
                end
            end
        end
        
        % pay attention: artificially + 1 to fit zhipeng's data
        logMax = max(max(dataNorm));

        logMax = ceil(logMax)+1;
        
        for i = 1: nVirus
            for j = 1: nSerum
                % not missing value
                if tableHI(i,j) ~=0
                    dataNorm(i,j) = logMax - dataNorm(i,j);
                end
            end
        end
        
        %% add the code for other normalization schemes
    case 2
        dataNorm = tableHI;
end

pwd;
currentFolder = pwd;

if ~exist('\Temp', 'dir')
    mkdir('\Temp')
end

normalizedTableFile= [currentFolder, '/Temp/normalizedTabe.txt'];

writeNormalTable(normalizedTableFile, dataNorm, virusName, serumName);
end
