function [lambda, fval] = tp_fitexpdecay(x,t,starting_values) 
% t needs to be in samples, e.g. 1:800
% x needs is the signal

x = x(:);
t = t(:);

f = @(pars) sum((x(:) - pars(2).*exp(pars(1).*t(:))).^2);

[lambda,fval] = fminsearch(f,starting_values);


