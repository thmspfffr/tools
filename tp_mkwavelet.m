
function [wavelet, freq_range, opt]= tp_mkwavelet(f,oct,fsample)

foi_min    = 2*f/(2^oct+1);
foi_max    = 2*f/(2^-oct+1);
delta_freq = foi_max-foi_min; % 2*std in freq domain
delta_time = 6/pi./delta_freq;
delta_time = round(delta_time*1000)/1000;
t_shift    = delta_time/2;
n_win      = round(delta_time*fsample);
n_shift    = round(t_shift*fsample);
TAPER      = gausswin(n_win,3)'; TAPER = TAPER/sum(TAPER);
iEXP       = exp(sqrt(-1) * ((1:n_win)-n_win/2-0.5) /fsample*f*2*pi);
wavelet    = (TAPER.*iEXP).';

freq_range = [foi_min foi_max];

opt.n_win = n_win;
opt.n_shift = n_shift;