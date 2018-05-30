function [H3] = hFive2Three(H5)
%% HFIVE2THREE converts the index in H5 to its corresponding index in H3.
%              must run const.m first to load data
%
%    USAGE:
%          HFIVE2THREE
%          HFIVE2THREE(H5)
%
%    HFIVE2THREE  gives help on HFIVE2THREE function.
%
%    HFIVE2THREE(H5) returns the corresponding H3N2 index for site H5.
%

%% Input checking
switch nargin
    case 0
        help hFive2Three
        return
end

const;
global CORR5H3;

%% converting process

if H5 > 324
    H3 = NaN;
else
    H3 = CORR5H3(H5);
end;

end