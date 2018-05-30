
seqSerum={158,1};
for i=1:79
   serumN=['>',serumName{i}];
   seqSerum{2*i-1}=serumN;
   flag='  ';
   for j=1:2:505
      if strcmp(H3N219682003{j},serumN)
          flag=H3N219682003{j+1};
          break;
      end
   end
   seqSerum{2*i}=flag;
end
T = cell2table(seqSerum');
writetable(T,'seqSerum.txt')