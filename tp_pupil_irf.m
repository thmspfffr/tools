function y = tp_pupil_irf(t,t_max,n) 
% Computes pupil impulse response function

% See Hoeks & Levelt (1993) for standard parameter:
% t_max = 0.930
% n = 10.1

y = t.^n .* exp(-(n*t./t_max));