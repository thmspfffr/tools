function [beta, fval] = tp_fitpowlaw(x,freq,starting_value) 
% t needs to be in samples, e.g. 1:800
% x needs is the signal

f = @(beta) sum((x - (beta(1).*freq.^beta(2)) ).^2);

[beta,fval] = fminsearch(f,starting_value);


