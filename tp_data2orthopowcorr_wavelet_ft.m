
function resout=tp_data2orthopowcorr_wavelet_ft(data,filt,para)
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

addpath ~/Documents/MATLAB/fieldtrip-20160919/; ft_defaults()

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
  
  dataf=dloc2'*KERNEL;
  datasf1(:,kk)=dataf'*filt;
  
end

resout = ft_connectivity_powcorr_ortho(datasf1);

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



