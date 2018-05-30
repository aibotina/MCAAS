function [virusCluster year] = yearAln(alnFile)

[seq seqName] = readAln(alnFile);


nVirus = numel(seqName);

tempYear = zeros(nVirus,1);

for i = 1: nVirus
    tempYear(i) = getYear(seqName{i});
end

[year yearPos yearInd]= unique(tempYear);

virusCluster = {};


for i=1: numel(year)
    currentYear = find(yearInd == i);
    virusCluster= [virusCluster; currentYear];
end
end