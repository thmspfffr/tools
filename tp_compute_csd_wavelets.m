function [csd,dataf,opt,conv_dataf] = tp_compute_csd_wavelets(data,para)
% tp_compute_csd_wavelets computes complex cross spectral density matrix based on morlets wavelets
% ----------
% INPUT
% ----------
% data: sensor level (Nchan x Ntimes)
% para.freq: center frequency [Hz]
% para.fsample: sampling rate
% ----------
% OUTPUT
% ----------
% csd: complex cross spectral density
% dataf: fourier coefficients
% opt: length and shift of wavelet (in samples)
% conv_dataf: as dataf, but w convolution yielding continous signal

if size(data,1)>size(data,2)
  error('Please check if data is of the correct format (Nchans x Ntimes)')
end
if ~isfield(para,'meth')
  para.meth = '';
end

% DEFINE WAVELETS
[wavelet,opt]=tp_mkwavelet(para.freq,0.5,para.fsample,para.overlap);
ss=repmat(wavelet,1,size(data,1));
nseg=floor((size(data,2)-opt.n_win)/opt.n_shift+1);

kk = 0;

for j=1:nseg
  
%   fprintf('%.3f%%\n',100*j/nseg)
  
  dloc2=data(:,(j-1)*opt.n_shift+1:(j-1)*opt.n_shift+opt.n_win);
  
  if any(isnan(dloc2(1,:)))
    warning('NaN detected')
    dataf(:,j) = nan(size(data,1),1);
    continue
  end
  
  kk=kk+1;
  
%   dataf(:,j)=dloc2*wavelet;
  dataf(:,j)=transpose(sum(dloc2'.*ss));
  
  if kk == 1
    csd = conj(dataf(:,j)*dataf(:,j)');
  else
    csd = csd + conj(dataf(:,j)*dataf(:,j)');
  end
  
end

% cs(:,:,f,j)+conj(datalocfft(f,:)'*datalocfft(f,:))

if strcmp(para.meth,'conv')
  conv_dataf = zeros(size(data,1),size(data,2));
  for isens = 1 : size(data,1)
    % --- compute complex signal --- %
    conv_dataf(isens,:) = conv(data(isens,:),wavelet','same');
  end
end
  

csd = csd/kk;
