function pars = tp_fitnakarushton(x,y,pars_init)

x = x(:);
y = y(:);

f = @(pars) sum((x - (pars(1) * (y.^pars(3) ./ (y.^pars(3) + pars(2).^pars(3))) + pars(4))).^2);

[pars,fval] = fminsearch(f,pars_init);


