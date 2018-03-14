function w = tp_mkwavelet(para)

% para.fsample
% para.freqoi
% para.width
% para.width


dt = 1/para.fsample;
  
endnsample = round(0 * para.fsample);

sf = para.freqoi / 5.83;
  st = 1/(2*pi*sf);
  toi2 = -para.width*st:dt:para.width*st;
  A = 1/sqrt(st*sqrt(pi));
  tap = (A*exp(-toi2.^2/(2*st^2)))';
  acttapnumsmp = size(tap,1);
  taplen = acttapnumsmp;
  ins = ceil(endnsample./2) - floor(acttapnumsmp./2);
  prezer = zeros(ins,1);
  pstzer = zeros(endnsample - ((ins-1) + acttapnumsmp)-1,1);
  
  % produce angle with convention: cos must always be 1  and sin must always be centered in upgoing flank, so the centre of the wavelet (untapered) has angle = 0
  ind  = (-(acttapnumsmp-1)/2 : (acttapnumsmp-1)/2)'   .*  ((2.*pi./para.fsample) .* para.freqoi);
  
  % create wavelet and fft it
  w = complex(vertcat(prezer,tap.*cos(ind),pstzer), vertcat(prezer,tap.*sin(ind),pstzer));
%   wltspctrm = complex(zeros(1,endnsample));
%   wltspctrm = fft(wavelet,[],1)';
  taplen = size(w,1);

%%
% 
% f = 2.^[1:.25:7]
% 
% f_lo = f/(2^(1/2))
% f_hi = (2^(1/2))*f
% 
% bp = [f' f_lo' f_hi']
% 
% sf = 1/(2*pi*st)
% 
% A = st*sqrt(pi)^(-1/2);
% toi2 = -para.width*st:dt:para.width*st;
% A*exp(-t.^2
