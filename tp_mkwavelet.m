
function [wavelet, freq_range, opt]= tp_mkwavelet(f,oct,fsample)
% 

foi_min    = 2*f/(2^oct+1);
foi_max    = 2*f/(2^-oct+1);
delta_freq = foi_max-foi_min; 
delta_time = 6/pi./delta_freq;
delta_time = round(delta_time*1000)/1000;
t_shift    = delta_time/2;
n_win      = round(delta_time*fsample);
n_shift    = round(t_shift*fsample);
A          = 1/sqrt((delta_time/6)*sqrt(pi));
tap        = A.*gausswin(n_win,3)';
iEXP       = exp(sqrt(-1) * ((1:n_win)-n_win/2-0.5) /fsample*f*2*pi);
wavelet    = (tap.*iEXP).';



% % 

% 
% % --- compute Kernel parameter --- %
% foi_min     = 2*f/(2^oct+1);   % using arithmetic mean: fc = (fmax + fmin)/2; -> get upper and lower limit of smoothing band from center frequency and bandwidth
% foi_max     = 2*f/(2^-oct+1);  % bandwidth (octave) between fmax and fmin NOT from foi down and up!
% std_freq    = (foi_max-foi_min)/(2*sqrt(2*log(2)));  % std in freq domain
% std_time    = 1./(2*pi*std_freq);        % std of KERNEL in time domain
% n_win_float = 6*std_time*fsample; % (how many samples does window have in time domain - entire window, without padding)
% n_win       = floor(n_win_float/2)*2+1; % round to next odd integer (even integers are rounded up)
% 
% % --- construct Kernel --- %
% TAPER   = jh_gausswin2(n_win,6/2)';
% TAPER   = TAPER/sum(TAPER);
% ti      = ((1:n_win)-n_win/2-0.5)/fsample; % time centered at 0
% iEXP    = exp(sqrt(-1)*2*pi*f*ti); % => kernel will be centered on 0
% wavelet  = (TAPER.*iEXP).'; % taper complex sinusoid of frequency foi with the gaussian window
%  lK      = length(wavelet);
% % 
opt.n_win = n_win;
opt.n_shift =round(n_win/2);
freq_range = [foi_min foi_max];
