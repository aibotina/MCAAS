function [virusName, virusColor, virusCorr] = readXml(fileName);
%READXML reads xml files

%% Input checking

if nargin == 0
    help readXml
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
    fText = strtrim(textread(fileName,'%s','delimiter','\n'));
else
    error('File %s does not exist!', fileName)
end

nLine = length(fText);
nVirus = nLine - 2; % nVirus denotes the number of viruses

virusName = cell(nVirus,1);
virusColor = cell(nVirus,1);
virusCorr = zeros(nVirus,3);

for i = 1: nVirus
    sLine = regexp(fText{i+2},del,'split'); % Split line i
    
    virusName(i) = sLine(1);    % Virus name is in the first col
    virusColor(i) = sLine(2);   % color name is in the second col
    
    virusCorr(i,1) = str2double(sLine{3});
    virusCorr(i,2) = str2double(sLine{4});
    virusCorr(i,3) = str2double(sLine{5});
end

end

