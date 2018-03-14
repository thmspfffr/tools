function x_filt = tp_filt(x,pb,fs,pad,varargin)

if any(any(strcmp(varargin,'ord')))
  ord=cell2mat(varargin(find(strcmp(varargin,'ord'))+1));
end

if any(any(strcmp(varargin,'meth')))
  meth=varargin(find(strcmp(varargin,'meth'))+1);
end

% 10s zero padding
x = [zeros(pad*fs,size(x,2)); x; zeros(pad*fs,size(x,2))];

delt    = 1/fs; 
fnq     = 1/(2*delt); 
Wn      = [pb(1)/fnq pb(2)/fnq]; 
[b,a]   = butter(ord,Wn);

if pb(2) == 0
  % high pass
elseif pb(2) == 0
  % low pass
else 
  % bandpass
  for i = 1 : size(x,2)
    x_filt(:,i) =  filtfilt(b,a,x(:,i));
  end
end

x_filt = x_filt(pad*fs+1:end-pad*fs,:);
  
  
  