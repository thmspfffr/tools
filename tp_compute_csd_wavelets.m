function [csd,nan_count] = tp_compute_csd_wavelets(data,para)
% compute complex cross spectral density matrix based on morlets wavelets
% data: sensor level (Nchan x Ntimes)
% para: freq = center frequency [Hz]
% para: fsample (sampling rate)

if size(data,1)>size(data,2)
  error('Please check if data is of the correct format (Nchans x Ntimes)')
end

% DEFINE WAVELETS
[KERNEL,~,opt]=tp_mkwavelet(para.freq,0.5,para.fsample);

nseg=floor((size(data,2)-opt.n_win)/opt.n_shift+1);

kk =0;
nan_count = 0;

for j=1:nseg
  
  fprintf('%.3f%%\n',100*j/nseg)
  
  dloc2=data(:,(j-1)*opt.n_shift+1:(j-1)*opt.n_shift+opt.n_win);
  
  if any(isnan(dloc2(1,:)))
    warning('NaN detected')
    nan_count = nan_count + 1;
    continue
  end
  
  kk=kk+1;
  
  dataf_4csd(:,j)=dloc2*KERNEL;
  
end

csd = dataf_4csd*dataf_4csd';
