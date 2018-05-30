function [virusName, elementType, virusCorr] = readXmlJmol(fileName)
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

virusName = {};
elementType = {};
virusCorr = [];


for i = 6: 2: nLine-2
    lineStr = fText{i};
    
    sLine = regexp(lineStr, del, 'split'); % Split line i
    
    % deal with virus name
    vName = sLine{2};
    nameLen = numel(vName);
    virusName = [virusName; vName(5: nameLen-1 )];
    
    % deal with elementType
    eType = sLine{3};
    typeLen = numel(eType);
    elementType = [elementType; eType(14:typeLen-1)];
    
    % deal with coordinates
    xCorrString = sLine{4};
    yCorrString = sLine{5};
    zCorrString = sLine{6};
    xLen = numel(xCorrString);
    yLen = numel(yCorrString);
    zLen = numel(zCorrString);
    
    
    xStr = xCorrString(5:xLen-1);
    xCorr = str2double(xStr);
    
    yStr = yCorrString(5:yLen-1);
    yCorr = str2double(yStr);
    
    zStr = zCorrString(5:zLen);
    zCorr = str2double(zStr);

    virusCorr = [virusCorr; xCorr yCorr zCorr];
end
end

