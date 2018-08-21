function [lambda, fval] = tp_fitexpdecay(x,t,starting_value) 
% t needs to be in samples, e.g. 1:800
% x needs is the signal

x = x(:);
t = t(:);

f = @(lambda) sum((x - exp(-lambda.*t)).^2);

[lambda,fval] = fminsearch(f,starting_value);


