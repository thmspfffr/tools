function y = tp_pupil_bandpass(x,bp,fsample,k) 
 k = 4;
 fnq = fsample/2;
 hil_hi = bp(1);
 hil_lo = 400;
 hil_Wn=[hil_hi/fnq hil_lo/fnq];
 [bhil, ahil] = butter(k, hil_Wn);
 
 tmp1 = filtfilt(bhil, ahil, x);

clear bhil ahil hil_Wn hil_hi hil_lo

k = 4;
hil_hi = bp(2);
hil_lo = 400;
fnq = fsample/2;
hil_Wn=[hil_hi/fnq hil_lo/fnq];
[bhil, ahil] = butter(k, hil_Wn);

tmp2 = filtfilt(bhil, ahil, x);

y = tmp1 - tmp2;