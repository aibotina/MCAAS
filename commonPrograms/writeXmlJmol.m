function writeXmlJmol(fileName, virusName, elementType, virusCorr)
%% writeXml writes the data into an Xml file

%% Input checking
% Parameters fileName, dataHI, virusName, serumName should be specified
if nargin == 0
    help writeXml
    return
elseif nargin < 4
    error('Parameters: fileName, virusName, virusColor, virusCorr')   
end

if ~ischar(fileName)
    error('FileName must be a character array!')
end

%% Write the data into file fileName
fid = fopen(fileName,'w+');

if fid == (-1)
    error('Can not open file %s for writing.\n Check wirte permission!',fileName)
end

nVirus = numel(virusName);

% write head lines
line1 = '<?xml version="1.0" encoding="ISO-8859-1"?>';
fprintf(fid,'%s\n', line1);
line2 = '<molecule id="METHANOL" xmlns="http://www.xml-cml.org/schema/cml2/core"';
fprintf(fid,'%s\n', line2);
line3 = 'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"';
fprintf(fid,'  %s\n', line3);
line4 = 'xsi:schemaLocation="http://www.xml-cml.org/schema/cml2/core cmlAll.xsd">';
fprintf(fid,'  %s\n', line4);
line5 = '<atomArray>';
fprintf(fid,'  %s\n', line5);

% Virus name element number and coordinate lines
for i = 1: nVirus    
    fprintf(fid, '<atom id="%s" elementType="%s" x3="%e" y3="%e" z3="%e\n', virusName{i}, elementType{i}, virusCorr(i,1), virusCorr(i,2), virusCorr(i,3));
    fprintf(fid,'"/>\n');
end

% End lines
fprintf(fid, '  </atomArray>\n');
fprintf(fid, '</molecule>');

fclose(fid);

end