function neighborSite = readNeighborSite(fileName)
%% readSite reads all amino acid site in the file fileName
%
%    readData  gives this help on the readData function.
%
%    neighborSite = readNeighborSite(fileName) reads neighbor sites from fileName
%
% Author: Jialiang Yang, CVM, MSU, jyang@cvm.msstate.edu
% revised 10/27/2011

%% Input checking

if nargin == 0
    help readNeighborSite
    return
else
    if ~ischar(fileName)
        error('FileName must be a character array!')
    end
    
    % Set default values
    del = '\s+';
end


%% Read data from file fileName

% Check if the file fileName exists
if (exist(fileName,'file') || exist(fullfile(pwd,fileName),'file'))
    % ftext reads lines from the file
    fText = textread(fileName,'%s','delimiter','\n');
else
    error('File %s does not exist!', fileName)
end

% Remove spaces at the begining and end and remove blank lines if exist
%--------------------------------------------------------------------------
fTrim = {};   % fTrim stores table after removing spances and blank lines 

for i = 1: length(fText)
    % Remove spaces at the begin and end
    lineTrim = strtrim(fText{i});   
    
    % Remove blank lines
    if ~isempty(lineTrim)
        fTrim = [fTrim; lineTrim];
    end
end    

% read sites and sotres them in vector allSite
%-------------------------------------------------------------------------
nRow = numel(fTrim);        % number of rows in the file

allSite = [];

for i = 1: nRow
    sLine = regexp(fTrim{i},del,'split'); % Split line i
    
    % sites located at 6th words
    if strcmpi(sLine{1},'ATOM')
        allSite = [allSite; str2double(sLine{6})];
    end
end

neighborSite = unique(allSite);
end