
function [resout, power]=tp_data2orthopowcorr_wavelet_weighted(data,f,filt)
% calculates power correlation after orthogonalization between two sets of brain
% voxels. Usually, these two sets are either identical (then the rsult is
% for all pairs) or one set consists of a single voxel which is taken as
% seed. Data are epoched, with all epeochs (i.e. trials) stacked on top of each other.
% Each epoch is further divided into segments. Coefficients in the Fourier-domain
% are calculated for each segment. Thus the segment length determines the
% frequency resolution.
%
% usage:
% [resout]=data2orthopowercorr(data,segleng,segshift,epleng,f,fsample,A1,A2)
%
% input:
% data: TxN data matrix for T time points and N channels
% segleng: length of a segment (in number of samples)
% segshift: amount of shift (in number of samples) (explanation: segments
%           are shifted within but not across epochs. E.g.
%           segshift=floor(segleng/2) leads to segments of  50% overlap
%            (or slightly less if segleng is odd))
% epleng: length of epoch (or trial) in number of trials)
% f:      frequency of interest (in Hz)
% fsample: sampling rate

% output
% resout: M1xM2 matrix of power correlations after orthogonalization.

scale=sqrt(nanmean(nanmean(data.trial.^2)));
data.trial=data.trial/scale;
[nchan ns1]=size(filt.sa.filt);
[epleng]=size(data.trial,2);
% ns1 = 90;
res1=zeros(ns1,ns1,'single');
res2=zeros(ns1,ns1,'single');
res3=zeros(ns1,ns1,'single');
res4=zeros(ns1,ns1,'single');
res5=zeros(ns1,ns1,'single');

kk=0;

for i=1:1
  
  % DEFINE WAVELETS
  % ----------------------------------
  octave = 0.5; % Frequency resolution
  % arithmetic mean
  foi_min    = 2*f/(2^octave+1);
  foi_max    = 2*f/(2^-octave+1);
  delta_freq = foi_max-foi_min; % 2*std in freq domain
  delta_time = 6/pi./delta_freq;
  delta_time = round(delta_time*1000)/1000;
  t_shift    = delta_time;
  n_win      = round(delta_time*data.fsample);
  n_shift    = round(t_shift*data.fsample);
  TAPER      = gausswin(n_win,3)'; TAPER = TAPER/sum(TAPER);
  iEXP       = exp(sqrt(-1) * ((1:n_win)-n_win/2-0.5) /data.fsample*f*2*pi);
  KERNEL     = (TAPER.*iEXP).';
  % ----------------------------------
  
  nseg=floor((size(data.trial,2)-n_win)/n_shift+1);
  
  %   reg_idx = 1:90;
  %   pos       = filt.sa.grid_aal4mm_indi;
  fprintf('Atlas distance weighting ...\n')
  %
  datasf1 = zeros(nseg,size(filt.sa.filt(:,:),2));
  nseg
  for j=1:nseg
%     j
%     100*j/nseg
%     tic
    %     fprintf('R%dS%d\n',ireg,j)
    dloc2=data.trial(:,(j-1)*n_shift+1:(j-1)*n_shift+n_win)';
    
    if any(isnan(dloc2(:,1)))
      warning('NaN detected')
      continue
    end
    
    dataf=dloc2'*KERNEL;
    datasf1(j,:)=dataf'*filt.sa.filt(:,:);
%     toc
  end
  datasf2 = datasf1; clear datasf1
  power = mean(abs(datasf2).^2,2);
  
  kk = 0;
  
  for j = 1:nseg
%     tic
    fprintf('%.3f%%\n',100*j/nseg)
    datasf1 = datasf2(j,:);
    
    kk = kk + 1;
    %     ns = 90;
    for i1=1:ns1
      %
      %       i1
      x1=datasf1(i1);
      x2=imag(datasf1*conj(x1)./abs(x1));
      y1=abs(x1)^2;
      y2=abs(x2).^2;
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
%     toc
  end
  
  res1=res1/kk;
  res2=res2/kk;
  res3=res3/kk;
  res4=res4/kk;
  res5=res5/kk;
  %
  resout=(res1-res2.*res3)./(sqrt((res4-res2.^2).*(res5-res3.^2))+eps);
  resout=(resout+resout')./2;
  
  % collapse across regions defined in 'labels'
  
  
%   for i = 1 : 90
%     for j = 1 : 90
%       if i == j; resout(i,j) = nan; continue; end
%       clear tmp
%       reg1 = find(filt.sa.aal_label==i);
%       reg2 = find(filt.sa.aal_label==j);
%       for ii = 1 : length(reg1)
%         tmp(ii,:,:,:) = squeeze(mean(atanh(resout(reg1(ii),reg2)),2));
%       end
%       fc(i,j) = squeeze(mean(tmp,1));
%     end
%   end
%   
%   resout = fc;
  
end

end


