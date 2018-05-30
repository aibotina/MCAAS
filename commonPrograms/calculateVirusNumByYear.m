function [ yearNum ] = calculateVirusNumByYear( alnFile )
%CALCULATEVIRUSNUMBYYEAR calculate the number of viruses in each year

[~, seqName] = readAln(alnFile);

nVirus = numel(seqName);

startYear = getYear(seqName{1});
endYear = getYear(seqName{nVirus});
yearSpan = endYear - startYear +1;

numYear = zeros(yearSpan,1);

currentYear = startYear;

count = 0;
k=1;

for i= 1: nVirus
    currentVirus = seqName{i};
    virusYear = getYear(currentVirus);
    
    if virusYear == currentYear
        count = count +1;
    else
        numYear(k) = count;
        k = k+1;
        currentYear = currentYear+1;
        count = 1;
    end
end

numYear(yearSpan) = count;

yearNum = zeros(yearSpan, 2);

for i = 1: yearSpan
    yearNum(i, 1) = startYear +i-1;
    yearNum(i, 2) = numYear(i);
end

end

