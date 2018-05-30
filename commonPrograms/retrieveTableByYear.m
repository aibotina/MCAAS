function retrieveTableByYear(fullTable, yearVec, subTable)
%% RETRIEVETABLEBYYEAR retrives a sub table from full tab by years 
%        specified in yearVec
%
%  Usage:
%        retrieveTableByYear gives this help.
%        retrieveTableByYear(fullTable, yearVec, subTable)
%
%  Input: 
%       fullTable: a tab deliminated file with the following format
%                    2 \tab 3
%                    ID     \tab  serum1 \tab serum2 \tab serum
%                 reference \tab     0   \tab    0   \tab   0
%                  virus1   \tab  data11 \tab data12 \tab data13
%                  virus2   \tab  data21 \tab data22 \tab data23
%        yearVec: the years of the table to retrieve
%       subTable: the output subtable file
%
%  Revision Date : May 4, 2012
%  Author: Jialiang Yang  CVM, MSU, jyang@cvm.msstate.edu

%% Input checking

if nargin == 0
    help retrieveTableByYear
    return
else
    if ~ischar(fullTable)
        error('FullTable must be a character array!')
    end
    
    if ~ischar(subTable)
        error('subTable must be a character array!')
    end
    
    if ~isnumeric(yearVec)
        error('yearVec must be a vector!');
    end
end

%% Retrieve sub table

disp('Read full Table start');
[dataHI, virusName, serumName, reference] = readTable(fullTable);
disp('Read full Table End');

nVirus = numel(virusName);
nSerum = numel(serumName);

subVirusInd = [];
subSerumInd = [];

% separate the viruses by year
for i = 1: nVirus
    currentVirusName = virusName{i};
    
    virusYear = getYear(currentVirusName);
    
    if any(virusYear == yearVec)
        subVirusInd = [subVirusInd i];
    end
end


% separate the serum by year
for i = 1: nSerum
    currentSerum = serumName{i};
    
    serumYear = getYear(currentSerum);
    
    if any(serumYear == yearVec)
        subSerumInd = [subSerumInd i];
    end
end

% write to subtable
subData = dataHI(subVirusInd,subSerumInd);
subName = virusName(subVirusInd);
subSerum = serumName(subSerumInd);
subRef = reference(subSerumInd);

writeTable(subTable, subData, subName, subSerum, subRef);

% end