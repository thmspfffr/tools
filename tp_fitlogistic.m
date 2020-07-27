function pars = tp_fitlogistic(x,y,pars_init)

x = x(:);
y = y(:);

f = @(pars) sum((x- (pars(1) ./ (1 + exp (-pars(2).*(y - pars(3)))))).^2)

[pars,fval] = fminsearch(f,pars_init);




