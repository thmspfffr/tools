
function [wavelet, outp]= tp_mkwavelet(f,oct,fsample,overlap)
% 
overlap = 1/(1-overlap);

f_min      = 2*f/(2^oct+1);
f_max      = 2*f/(2^-oct+1);
delta_freq = f_max-f_min; 
delta_time = 6/pi./delta_freq;
delta_time = round(delta_time*1000)/1000;
n_win      = round(delta_time*fsample);
A          = 1/sqrt((delta_time/6)*sqrt(pi));
tap        = A.*gausswin(n_win,3)';
iEXP       = exp(sqrt(-1) * ((1:n_win)-n_win/2-0.5) /fsample*f*2*pi);
wavelet    = (tap.*iEXP).';
 

outp.n_win      = n_win;
outp.n_shift    = round(n_win/overlap);
outp.freq_range = [f_min f_max];
