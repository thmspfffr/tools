function y = tp_pupil_irf(t,t_max,n) 
% Computes pupil impulse response function
% y = tp_pupil_irf(t,t_max,n)

% See Hoeks & Levelt (1993) for standard parameter:
% t_max = 0.930
% n = 10.1

if isempty(t_max)
  t_max = 0.930;
end
if isempty(n)
  n = 10.1;
end

y = t.^n .* exp(-(n*t./t_max));