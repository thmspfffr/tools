function f = resp_fun(x,a)

% c = 100;
% a = 0.1;
% b = 40;
% x = 0:0.01:100;

% f = 1./(1+exp(-(a*x+b)));
f = 1./(1+exp(-x/a));

% f = 1./(1+exp(-x));

% plot(x,f)