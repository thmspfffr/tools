function par = tp_ARfit(y,ord)
% Performs fit of 1st order autoregressive model
% x is time series
% ord is the model parameter

if ord == 1
  x     = ones(length(y),1);

  R  = [y(1:end-1)' x(1:end-1)];
  Y  = [y(2:end)]';

  par = R\Y;
end

