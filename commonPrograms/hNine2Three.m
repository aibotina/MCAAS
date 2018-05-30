function [H3] = hNine2Three(H9)
%% HFIVE2THREE converts the index in H9 to its corresponding index in H3.
%              must run const.m first to load data
%
%    USAGE:
%          HNINE2THREE
%          HNINE2THREE(H9)
%
%    HFIVE2THREE  gives help on HNINE2THREE function.
%
%    HNINE2THREE(H9) returns the corresponding H3N2 index for site H9.
%

%% Input checking
switch nargin
    case 0
        help hNine2Three
        return
end

const;
global CORR9H3;

%% converting process

if H9 > 324
    H3 = NaN;
else
    H3 = CORR9H3(H9);
end;

end