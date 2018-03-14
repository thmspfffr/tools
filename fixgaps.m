function y=fixgaps(x, method)
% FIXGAPS Linearly interpolates gaps in a time series
% YOUT=FIXGAPS(YIN) linearly interpolates over NaN
% in the input time series (may be complex), but ignores
% trailing and leading NaN.
%

% R. Pawlowicz 6/Nov/99

y=x;

bd=isnan(x);
gd=find(~isnan(x));

bd([1:(min(gd)) (max(gd)):end])=0;

y(bd)=interp1(gd,x(gd),find(bd), method);