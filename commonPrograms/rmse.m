function rootMeanSquareError = rmse(y,yEstimate)
%% rmse calculates the Root mean squared error of y and yEstimate
%
%    USAGE:
%          rootMeanSquareError =  rmse(y,yEstimate)
%    INPUT:
%          y,yEstimate: two numerical vectors of the same size
%    OUTPUT:
%          rootMeanSquareError: the root mean square error of y and yEstimate
%
%% Calculate rootMeanSquareError

rootMeanSquareError = sqrt(mean((yEstimate - y).^2));

end