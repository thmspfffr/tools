function [c, coh] = tp_powcorr_ortho(data,para,sa)

% function computes orthogonalized amplitude envelope correlations and
% phase coherence in source space (see hipp et al., 2012, nat. neurosci. for
% reference). input should be organized the following way:
% data: [time x channels]
% para:
%   - segleng:  length of segments for orthoginalization. a short segleng
%               is recommended in order to avoid problems arising from
%               nonstationarities in the data.
%   - epleng:   length of total epoch. if epleng == size(data,1), the
%               entire dataset will be analyzed. if smaller epleng is
%               chosen, time-variant FC can be computed.
%   - wavelet:  *CHANGE NAME* 3 options are available:
%                     - 'hanning' computes wavelets as impl. by guido
%                     - 'ft' computes morlet's as impl. by ft
%                     - 'bp_filt' uses 4th order butterworth bandpass filt.
%   - scnd_filt: is either 0 or 1. if 1, amplitude envelopes are filtered
%                a second time in the range from 0.04 to 0.2 hz

% tpfeffer (2016), thms.pfffr@gmail.com

if ~exist('para','var')
  para.segave = 1;
end

epleng    = para.epleng;
epshift   = para.epshift;
segleng   = para.segleng;
fsample   = para.fsample;
segshift  = para.segshift;
f         = para.foi;
scale     = sqrt(mean(mean(data.^2)));
data      = data/scale;
[n nchan] = size(data);

switch para.wavelet
  
  case 'hanning'
    
    nn    = (1:segleng)'-segleng/2;
    mywin = hanning(segleng);
    s1 = cos(nn*f*2*pi/fsample).*mywin;
    s2 = sin(nn*f*2*pi/fsample).*mywin;
    ss = s1-sqrt(-1)*s2;
    ss = repmat(ss,1,nchan);
    
  case 'ft'
    
    para.fsample = fsample;
    para.freqoi = f;
    para.gwidth = 3;
    para.width  = 7;
    
    w = tp_mkwavelet(para);
    ss = repmat(w,1,nchan);
    segleng = length(w);
    segshift = floor(segleng/8);
    
  case 'bp_filt'
    
    para.fsample = fsample;
    
    if any(size(f)~=1)
      para.freqoi = [f(1) f(2)];
      flp = f(1);           % lowpass frequency of filter
      fhi = f(2);
    else
      para.freqoi = [f-1 f+1];
      flp = f-1;           % lowpass frequency of filter
      fhi = f+1;
    end
    
    para.ord = 4;
    delt = 1/fsample;            % sampling interval
    k=4;                  % 2nd order butterworth filter
    fnq=1/(2*delt);       % Nyquist frequency
    Wn=[flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
    [bfilt,afilt]=butter(k,Wn);
    
end

nep  = floor((n-epleng)/epshift+1);

for iep = 1 : nep
  
  dloc1 = data((iep-1)*epshift+1:(iep-1)*epshift+epleng,:);
  
  if nep == 1; clear data; end
  
  % ----------------------------------------
  % compute cross spectrum
  % ----------------------------------------
  cs_data1    = [zeros(5000,size(dloc1,2)); dloc1(1:end-800,:); zeros(5000,size(dloc1,2))];
  cs_data2    = hilbert(zscore(filtfilt(bfilt,afilt,dloc1)));
  cs_data2    = cs_data2(5001:end-5000,:);
  
  segleng  = fsample;
  segshift = segleng / 2;
  epleng   = size(cs_data2,1);
  nseg     = floor((epleng-segleng)/segshift+1);
  
  cs = zeros(size(dloc1,2),size(dloc1,2)); 
  
  fprintf('Computing cross spectrum ...\n'); 
  
  for iseg = 1 : nseg
    
    cs_data3 = cs_data2((iseg-1)*segshift+1:(iseg-1)*segshift+segleng,:);
    cs = cs+(cs_data3'*cs_data3/size(cs_data3,2))/nseg;
    
  end
  clear cs_data3 cs_data2 cs_data1
  % ----------------------------------------
  
  % ----------------------------------------
  % COMPUTE SPATIAL FILTER
  % ----------------------------------------
  
  if strcmp(para.filt,'jh_lcmv')
    
    para.iscs = 1;
    if isfield(sa.sa,'L_coarse')
      % note that this means cortex3000
      load('~/Documents/MATLAB/aalmask_grid_cortex3000.mat')
      filt      = pconn_beamformer(cs,sa.sa.L_coarse,para);
      filt      = filt(:,find(aalgrid.mask));
      pos       = sa.sa.grid_cortex3000_indi;
    elseif  isfield(sa.sa,'L_aal')
      filt      = pconn_beamformer(cs,sa.sa.L_aal,para);
    elseif  isfield(sa.sa,'L_medium')
      load ~/pconn/matlab/aalmask_grid_medium
      filt      = pconn_beamformer(cs,sa.sa.L_medium,para);
      filt      = filt(:,find(aalgrid.mask));
      pos       = sa.sa.grid_medium_indi;
    end
     pos       = pos(find(aalgrid.mask),:);
     
  elseif strcmp(para.filt,'eloreta')
    pars      = [];
    pars.filt = 'eloreta';
    pars.cs   = cs;
    pars.foi  = f;
    pars.sa   = sa.sa;
    if isfield(sa.sa,'L_coarse')
      pars.grid = 'cortex';
      pos = sa.sa.grid_cortex3000_indi;
      load('~/Documents/MATLAB/aalmask_grid_cortex3000.mat')
    elseif  isfield(sa.sa,'L_aal')
      pars.grid = 'aal';
    elseif  isfield(sa.sa,'L_medium')
      pos = sa.sa.grid_medium_indi;
      pars.grid = 'medium';
      load ~/pconn/matlab/aalmask_grid_medium
    end
    filt      = get_spatfilt(pars);
    filt      = filt(:,find(aalgrid.mask));
    pos       = pos(find(aalgrid.mask),:);
  end
  % ----------------------------------------

  mom = filt'*dloc1'; clear dloc1
 
  % --------------------------------------
  % COMPUTE GAUSSIAN WEIGHTING (see Brookes et al., 2016)
  % --------------------------------------

  para.grid = 'medium';
  aal_pos = tp_grid2aal(sa.sa.grid_medium_indi',para)';
  aalgrid.mask=aalgrid.mask(aalgrid.mask~=0);

  aal_mom = zeros(91,size(mom,2));
  fprintf('AAL distance weighting ...\n')
  for ireg = 1 : 91
    
    idx = find(aalgrid.mask==ireg);
    dist = sqrt([aal_pos(ireg,1)-pos(idx,1)].^2 + [aal_pos(ireg,2)-pos(idx,2)].^2 + [aal_pos(ireg,3)-pos(idx,3)].^2);
    w = exp((-(10*dist).^2)./400);
    
    for imom = 1 : length(idx)
      aal_mom(ireg,:) = aal_mom(ireg,:)+(w(imom)*mom(imom,:)./length(idx));
    end
    
  end
  fprintf('AAL done ...\n')
  mom = aal_mom; clear aal_mom
  
  % COMPUTE CORRELATIONS BASED ON BAND-PASS FILTERED SIGNAL
  if strcmp(para.wavelet,'bp_filt')
    

    mom = (hilbert(zscore(filtfilt(bfilt,afilt,mom'))))'; 

    if para.scnd_filt

      flp   = 0.02;           % lowpass frequency of filter
      fhi   = 0.4;           % highpass
      res       = 2;
      delt      = 1/(fsample/res);            % sampling interval
      k         = 2;                  % 2nd order butterworth filter
      fnq       = 1/(2*delt);       % Nyquist frequency
      Wn        = [flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
      [bf2,af2] = butter(k,Wn);   % construct the filter
      
      fprintf('Filter, resample, filter again ...\n')

      if res ~= 1
        for il = 1 : size(mom,1)
          rmom(il,:) = resample(abs(mom(il,:)),1,res);
        end
      end
      
      mom = rmom; clear rmom
      
      mom    = hilbert(zscore(filtfilt(bf2,af2,mom)));
      
      if para.tau == inf
        segleng  = size(mom,2);
      else
        segleng  = para.tau*(fsample/res);
      end
      
      segshift = segleng / 2;
      epleng   = size(mom,2);
      nseg     = floor((epleng-segleng)/segshift+1);
      
    end
  
    fprintf('Computing orth. amp. correlations ...\n')
    
    % define the tau for orthogonalizing
    if nep == 1
      for iseg = 1 : nseg
        tmp_mom = double(mom(:,(iseg-1)*segshift+1:(iseg-1)*segshift+segleng));
        % this one computes orthopowcorrs based on fieldtrip code
        c(:,:,iseg) =  compute_orthopowcorr(tmp_mom);
      end
      c = nanmean(c,3); clear mom
    else
      c(:,:,iep) =  compute_orthopowcorr(mom); clear mom
    end
    
  else
    
    para.iscs = 0;
    
    for iseg = 1 : nseg
      
      fprintf('Seg %d/%d ...\n',iseg,nseg)
      
      dloc2 = dloc1((iseg-1)*segshift+1:(iseg-1)*segshift+segleng,:);
      
      dataf = sum(dloc2.*ss);
      mom(:,iseg) = dataf;
      
    end
    
    para.iscs = 0;
    filt = pconn_beamformer(mom,sa.sa.L_aal,para);
    
    mom = filt'*mom;
    
    % compute cross spectrum in source for coherence
    cs = complex(zeros(91,91));
    
    for iseg = 1 : nseg
      
      dloc3 = mom(:,(iseg-1)*segshift+1:(iseg-1)*segshift+segleng);
      dloc3 = dloc3';
      cs = cs+(dloc3'*dloc3/size(dloc3,2))/nseg;
      
    end
    
    coh = cs./sqrt(diag(cs)*diag(cs)');
    
    % this one computes orthopowcorrs based on fieldtrip code
    c(:,:,iep) =  compute_orthopowcorr(mom); clear mom
    
  end
end
% compute power corr


function c = compute_orthopowcorr(mom,varargin)

refindx = 'all'; clear tapvec
tapvec  = ones(1,size(mom,2));

if strcmp(refindx, 'all')
  refindx = 1:size(mom,1);
end

cmomnorm = conj(mom./abs(mom)); % only need to do conj() once

n        = size(mom,1);
ntap     = 1;
if ~all(tapvec==ntap)
  error('unequal number of tapers per observation is not yet supported');
end

ix = zeros(sum(tapvec),1);
jx = ix;
sx = ix;

for k = 1:numel(tapvec)
  indx = (k-1)*ntap+(1:ntap);
  ix(indx) = indx;
  jx(indx) = k;
  sx(indx) = 1./ntap;
end

tra = sparse(ix,jx,sx,sum(tapvec),numel(tapvec));

powmom = (abs(mom).^2)*tra; % need only once
powmom = standardise(log10(powmom), 2);

c = zeros(n, numel(refindx));%;*2);
N = ones(n,1);
%warning off;
for k = 1:numel(refindx)
  indx     = refindx(k);
  ref      = mom(indx,:);
  crefnorm = conj(ref./abs(ref));
  
  pow2 = (abs(imag(ref(N,:).*cmomnorm)).^2)*tra;
  pow2 = standardise(log10(pow2), 2);
  c1   = mean(powmom.*pow2, 2);
  pow1 = (abs(imag(mom.*crefnorm(N,:))).^2)*tra;
  pow1 = standardise(log10(pow1), 2);
  
  pow2 = (abs(ref).^2)*tra;
  pow2 = standardise(log10(pow2), 2);
  pow2 = repmat(pow2, [n 1]);
  c2   = mean(pow1.*pow2, 2);
  
  c(:,k) = (c1+c2)./2;
  %c(:,k+numel(refindx)) = c2;
end


