
function resout=tp_data2orthopowcorr_wavelet_sens(data,para)
% calculates power correlation after orthogonalization between two sets of brain
% voxels. Usually, these two sets are either identical (then the rsult is
% for all pairs) or one set consists of a single voxel which is taken as
% seed. Data are epoched, with all epeochs (i.e. trials) stacked on top of each other.
% Each epoch is further divided into segments. Coefficients in the Fourier-domain
% are calculated for each segment. Thus the segment length determines the
% frequency resolution.
%
% usage:
% [resout]=tp_data2orthopowcorr_wavelet(data,f,filt)
%
% input:
% data: TxN data matrix for T time points and N channels
% f:    Frequency of interest (in Hz)
% filt: Spatial filter
%
% output
% resout: M1xM2 matrix of power correlations after orthogonalization.

scale=sqrt(nanmean(nanmean(data.^2)));
data=data./scale;
% ----------------------------------
res1=zeros(size(data,1),size(data,1),'single');
res2=zeros(size(data,1),size(data,1),'single');
res3=zeros(size(data,1),size(data,1),'single');
res4=zeros(size(data,1),size(data,1),'single');
res5=zeros(size(data,1),size(data,1),'single');
% ----------------------------------
% DEFINE WAVELETS
% ----------------------------------
octave = 0.5; % Frequency resolution
[KERNEL,f,opt]=tp_mkwavelet(para.freq,octave,para.fsample);
% ----------------------------------

nseg=floor((size(data,2)-opt.n_win)/opt.n_shift+1);

kk = 0;
for j=1:nseg
  
  fprintf('%.3f%%\n',100*j/nseg)
  
  dloc2=data(:,(j-1)*opt.n_shift+1:(j-1)*opt.n_shift+opt.n_win)';
  
  if any(isnan(dloc2(:,1)))
    warning('NaN detected')
    continue
  end
  kk=kk+1;
  
  dataf=(dloc2'*KERNEL)';
  
  for i1=1:size(dloc2,2)
    x1=dataf(i1);
    x2=imag(dataf*conj(x1)./abs(x1));
    y1=log10(abs(x1)^2);
    y2=log10(abs(x2).^2);
    if kk==1
      res1(i1,:)=y1*y2;
      res2(i1,:)=y1;
      res3(i1,:)=y2;
      res4(i1,:)=y1^2;
      res5(i1,:)=y2.^2;
    else
      res1(i1,:)=res1(i1,:)+y1*y2;
      res2(i1,:)=res2(i1,:)+y1;
      res3(i1,:)=res3(i1,:)+y2;
      res4(i1,:)=res4(i1,:)+y1^2;
      res5(i1,:)=res5(i1,:)+y2.^2;
    end
  end
  
end
% -------------------------------------
res1=res1/kk;
res2=res2/kk;
res3=res3/kk;
res4=res4/kk;
res5=res5/kk;
% -------------------------------------
% resout is asymetrical (see hipp 2012), resout needs to be averaged
resout=(res1-res2.*res3)./(sqrt((res4-res2.^2).*(res5-res3.^2))+eps);
resout=tanh(atanh(resout)./2+atanh(resout)'./2);
% -------------------------------------


