function correlationCoefficent = cc(y,yEstimate)
%% cc calculates the correlation coefficient of y and yEstimate
%
%    USAGE:
%          correlationCoefficent =  cc(y,yEstimate)
%    INPUT:
%          y,yEstimate: two numerical vectors of the same size
%    OUTPUT:
%          correlationCoefficent: the correlation coefficent of y and yEstimate
%
%% Calculate correlationCoefficent

lenY = numel(y);

correlationCoefficent  = (lenY*sum(yEstimate.* y)-sum(y)*sum(yEstimate))/ sqrt((lenY*sum(y.^2)-sum(y)^2)*(lenY*sum(yEstimate.^2)-sum(yEstimate)^2));

end