function [resout,variance] = tp_data2orthopowcorr_wavelet(data,filt,para)
% calculates power correlation after orthogonalization between two sets of brain
% voxels. Each epoch (usually the entire data set) is  divided into 
% segments. Coefficients in the Fourier-domain are calculated for each segment
% using Morelets wavelets.
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
%
% Implementation of Hipp et al. (2012) Nature Neuroscience
% Original code by Guido Nolte, UKE Hamburg
% Adapted by Thomas Pfeffer, UKE Hamburg
% ----------------------------------
% Note: this implementaiton yields numerically similar (equivalent)
% results as fieldtrip, but the computation is faster.
% (original implementation by G. Nolte does not include log-transform)
% ----------------------------------

scale=sqrt(nanmean(nanmean(data.^2)));
data=data./scale;
% ----------------------------------
res1=zeros(size(filt,2),size(filt,2),'single');
res2=zeros(size(filt,2),size(filt,2),'single');
res3=zeros(size(filt,2),size(filt,2),'single');
res4=zeros(size(filt,2),size(filt,2),'single');
res5=zeros(size(filt,2),size(filt,2),'single');
% ----------------------------------
% DEFINE WAVELETS
% ----------------------------------
octave = 0.5; % Frequency resolution
[wavelet,f,opt]=tp_mkwavelet(para.freq,octave,para.fsample);
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
  
  dataf=dloc2'*wavelet;
  datasf1=dataf'*filt;
  
  for i1=1:size(filt,2)
    x1=datasf1(i1);
    x2=imag(datasf1*conj(x1)./abs(x1));
% log-transform power as in hipp 2012 (10-03-2020)
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
  
  datvar(j,:) = abs(datasf1).^2;
  
end

variance = var(datvar);
% -------------------------------------
res1=res1/kk;
res2=res2/kk;
res3=res3/kk;
res4=res4/kk;
res5=res5/kk;
% -------------------------------------
% resout is asymetrical (see hipp 2012), resout needs to be averaged
resout=(res1-res2.*res3)./(sqrt((res4-res2.^2).*(res5-res3.^2))+eps);
resout=tanh((atanh(resout)./2 + atanh(resout)'./2));
% -------------------------------------
% REDUCE IF SIZE OF GRID IS INDICATIVE OF AAL PARCELATION
% (in this case, the number of vertices is >5000)
% Average across regions
% -------------------------------------
if size(resout,1)>5000
  load ~/pupmod/proc/aal_labels.mat
  reg_idx = 1:90;
  fprintf('Combining atlas regions ...\n')
  for ireg = reg_idx
    ireg
    idx{ireg} = find(aal_labels==ireg);
  end
  for i = 1 : size(idx,2)
    for j = 1 : size(idx,2)
      if i == j; fc(i,j) = nan; continue; end
      clear tmp
      for ii = 1 : length(idx{i})
        tmp(ii) = squeeze(tanh(mean(atanh(resout(idx{i}(ii),idx{j})),2)));
      end
      fc(i,j) = squeeze(tanh(mean(atanh(tmp))));
    end
  end
  resout=fc;
end
% -------------------------------------



