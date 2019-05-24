
function resout=tp_data2orthopowcorr_wavelet(data,f,filt)
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
% nn=(1:segleng)'-segleng/2;

% mywin=hanning(segleng);
% s1=cos(nn*f*2*pi/fsample).*mywin;
% s2=sin(nn*f*2*pi/fsample).*mywin;
% ss=s1-sqrt(-1)*s2;
% if nargin<8;
%     A2=A1;
% end
[nchan ns1]=size(filt);
% [nchan ns2]=size(A2);

[epleng]=size(data.trial,2);
% ss=repmat(ss,1,nchan);

% nep=floor(data.trial/epleng);


res1=zeros(ns1,ns1);
res2=zeros(ns1,ns1);
res3=zeros(ns1,ns1);
res4=zeros(ns1,ns1);
res5=zeros(ns1,ns1);

kk=0;

%datarawf=[];
for i=1:1
  
%   data.dat = data.trial;
%   data.avg = data.dat; data.trial = [];
%   data.dimord = 'chan_time';
%   data.time = [];
%   data.time =  1/data.fsample:1/data.fsample:size(data.avg,2)/data.fsample;
% 
%   cfg=[];
%   cfg.method  ='wavelet';
%   cfg.output  = 'fourier';
%   cfg.channel = {'MEG'};
  cfg.foi     = f;
  cfg.width   = 5.83; % again, as per Hipp et al. (2012) Nat Neurosci
  tempSD      = 1./(2*pi*(cfg.foi./cfg.width)); % temporal SD in sec
  tempWd      = round(3*tempSD*data.fsample)/data.fsample; % set to 3 std dev, comp with 1000 Hz sampl rate
  cfg.toi     = tempWd.*(1:floor(data.time(end)./tempWd));
% 
%   cfg.pad='nextpow2';
%   tf=ft_freqanalysis(cfg,data);
% %   tf.fourierspctrm(:,:,:,end)=[];
%   
%   tf = squeeze(tf.fourierspctrm);
%   tf(:,isnan(tf(1,:)))=[];

  gwidth = 3;

  dt = 1/data.fsample;
  sf = f /  cfg.width;
  st = 1/(2*pi*sf);
  toi2 = -gwidth*st:dt:gwidth*st;
  A = 1/sqrt(st*sqrt(pi));
  tap = (A*exp(-toi2.^2/(2*st^2)))';
  acttapnumsmp = size(tap,1);
  ind  = (-(acttapnumsmp-1)/2 : (acttapnumsmp-1)/2)'   .*  ((2.*pi./data.fsample) .* f);

  s1=tap.*cos(ind);
  s2=tap.*sin(ind);
  ss = s1-sqrt(-1)*s2;
  ss=repmat(ss,1,nchan);
  
  segleng = acttapnumsmp;
  segshift = diff(cfg.toi)*data.fsample; segshift = segshift(1);
  
  nseg=(epleng-segleng)/segshift+1;

%   dataf = tf'*filt;
  kk =0;
%   resout = compute_orthopowcorr(dataf');
  for j=1:nseg
    j
    kk=kk+1;
%     dloc2=data.trial(:,(idx(j)-1)*segshift+1:(idx(j)-1)*segshift+segleng)';
        dloc2=data.trial(:,(j-1)*segshift+1:(j-1)*segshift+segleng)';

%     if any(isnan(dloc2(:,1)))
%       warning('NaN detected')     
%       continue
%     end
%     
    dataf=sum(dloc2.*ss)';
    datasf1=dataf'*filt;
    %datarawf=[datarawf;datasf1]
% end
    
%     end
%     datasf1=dataf(j,:);
    for i1=1:ns1;
%       
      x1=datasf1(i1);
      x2=imag(datasf1*conj(x1)./abs(x1));
      y1=abs(x1)^2;
      y2=abs(x2).^2;
      if kk==1;
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
% end
% 
res1=res1/kk;
res2=res2/kk;
res3=res3/kk;
res4=res4/kk;
res5=res5/kk;
% 
resout=(res1-res2.*res3)./(sqrt((res4-res2.^2).*(res5-res3.^2))+eps);
resout=(resout+resout')./2;
end
return;
% 
% function c = compute_orthopowcorr(mom,varargin)
% 
% refindx = 'all'; clear tapvec
% tapvec  = ones(1,size(mom,2));
% 
% if strcmp(refindx, 'all')
%   refindx = 1:size(mom,1);
% end
% 
% cmomnorm = conj(mom./abs(mom)); % only need to do conj() once
% 
% n        = size(mom,1);
% ntap     = 1;
% if ~all(tapvec==ntap)
%   error('unequal number of tapers per observation is not yet supported');
% end
% 
% ix = zeros(sum(tapvec),1);
% jx = ix;
% sx = ix;
% 
% for k = 1:numel(tapvec)
%   indx = (k-1)*ntap+(1:ntap);
%   ix(indx) = indx;
%   jx(indx) = k;
%   sx(indx) = 1./ntap;
% end
% 
% tra = sparse(ix,jx,sx,sum(tapvec),numel(tapvec));
% 
% powmom = (abs(mom).^2)*tra; % need only once
% powmom = standardise(log10(powmom), 2);
% 
% c = zeros(n, numel(refindx));%;*2);
% N = ones(n,1);
% %warning off;
% for k = 1:numel(refindx)
%   indx     = refindx(k);
%   ref      = mom(indx,:);
%   crefnorm = conj(ref./abs(ref));
%   
%   pow2 = (abs(imag(ref(N,:).*cmomnorm)).^2)*tra;
%   pow2 = standardise(log10(pow2), 2);
%   c1   = mean(powmom.*pow2, 2);
%   pow1 = (abs(imag(mom.*crefnorm(N,:))).^2)*tra;
%   pow1 = standardise(log10(pow1), 2);
%   
%   pow2 = (abs(ref).^2)*tra;
%   pow2 = standardise(log10(pow2), 2);
%   pow2 = repmat(pow2, [n 1]);
%   c2   = mean(pow1.*pow2, 2);
%   
%   c(:,k) = (c1+c2)./2;
%   %c(:,k+numel(refindx)) = c2;
% end

