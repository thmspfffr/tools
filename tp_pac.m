function pac = tp_pac(x,f_lo,f_hi,fs) 


len=length(x); %% number of sample points in raw signal

minskip=fs; %% time lag must be at least this big

maxskip=len-fs; %% time lag must be smaller than this

% skip=ceil(numpoints.*rand(len*2,1));
% skip(find(skip>maxskip))=[];
% skip(find(skip<minskip))=[];
% skip=skip(1:len,1);
% surrogate_m=zeros(len,1); 

% analytic amplitude time series, uses eegfilt.m from EEGLAB toolbox
amplitude=abs(hilbert(eegfilt(x,fs,f_hi(1),f_hi(2))));
f
% gamma analytic phase time series, uses EEGLAB toolbox
phase=angle(hilbert(eegfilt(x,fs,f_lo(1),f_lo(2))));

% complex-valued composite signal
z=amplitude.*exp(i*phase);

% mean of z over time, prenormalized value
pac=mean(z);