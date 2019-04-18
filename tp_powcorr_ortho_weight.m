function [c, coh] = tp_powcorr_ortho_weight(data,para,sa)

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
    
    error('Does not work!')
    nn    = (1:segleng)'-segleng/2;
    mywin = hanning(segleng);
    s1 = cos(nn*f*2*pi/fsample).*mywin;
    s2 = sin(nn*f*2*pi/fsample).*mywin;
    ss = s1-sqrt(-1)*s2;
    ss = repmat(ss,1,nchan);
    
  case 'ft'
    
    error('Does not work!')
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
    
    para.freqoi = [para.bpfreq(1) para.bpfreq(2)];
    flp = para.bpfreq(1);           % lowpass frequency of filter
    fhi = para.bpfreq(2);
    
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
  if nep ~= 1
    cs_data2    = hilbert(zscore(filtfilt(bfilt,afilt,dloc1)));
  else
    cs_data1    = [zeros(5000,size(dloc1,2)); dloc1; zeros(5000,size(dloc1,2))];
   	cs_data2    = hilbert(zscore(filtfilt(bfilt,afilt,cs_data1)));
    cs_data2    = cs_data2(5001:end-5000,:);
  end
  
  segleng  = fsample;
  segshift = segleng / 2;
  epleng   = size(cs_data2,1);
  nseg     = floor((epleng-segleng)/segshift+1);
  
  cs       = zeros(size(dloc1,2),size(dloc1,2));
  
  fprintf('Computing cross spectrum ...\n');
  
  for iseg = 1 : nseg
    
    cs_data3 = cs_data2((iseg-1)*segshift+1:(iseg-1)*segshift+segleng,:);
    cs       = cs+(cs_data3'*cs_data3/size(cs_data3,2))/nseg;
    
  end
  clear cs_data3 cs_data2 cs_data1
  % ----------------------------------------
  
  % ----------------------------------------
  % COMPUTE SPATIAL FILTER
  % ----------------------------------------
  
  if strcmp(para.filt,'jh_lcmv')
    
    para.iscs = 1;
    if strcmp(para.grid,'cortex')
      % note that this means cortex3000
      load('~/Documents/MATLAB/aalmask_grid_cortex3000.mat')
      filt      = pconn_beamformer(cs,sa.sa.L_coarse,para);
      filt      = filt(:,find(aalgrid.mask));
      pos       = sa.sa.grid_cortex3000_indi;
    elseif  strcmp(para.grid,'aal')
      filt      = pconn_beamformer(cs,sa.sa.L_aal,para);
    elseif  strcmp(para.grid,'medium')
      load ~/pconn/matlab/aalmask_grid_medium
      filt      = pconn_beamformer(cs,sa.sa.L_medium,para);
      filt      = filt(:,find(aalgrid.mask));
      pos       = sa.sa.grid_medium_indi;
    elseif strcmp(para.grid,'aal_6mm')
      pos       = sa.sa.grid_aal6mm_indi;
      filt      = pconn_beamformer(cs,sa.sa.L_aal_6mm,para);
    elseif strcmp(para.grid,'aal_4mm')
      pos       = sa.sa.grid_aal4mm_indi;
      filt      = pconn_beamformer(cs,sa.sa.L_aal_4mm,para);
    elseif strcmp(para.grid,'m758_4mm')
      pos       = sa.sa.grid_m758_4mm_indi;
      filt      = pconn_beamformer(cs,sa.sa.L_m758_4mm,para);
      m758      = tp_m758_grid();
      sa.sa.aal_label = m758.tissue_4mm(m758.tissue_4mm>0);
    elseif strcmp(para.grid,'m758_6mm')
      pos       = sa.sa.grid_m758_6mm_indi;
      filt      = pconn_beamformer(cs,sa.sa.L_m758_6mm,para);
      m758      = tp_m758_grid();
      sa.sa.aal_label = m758.tissue_4mm(m758.tissue_4mm>0);
    elseif strcmp(para.grid,'vtpm_4mm')
       pars.grid = 'L_vtpm_4mm';
       filt      = pconn_beamformer(cs,sa.sa.L_vtpm_4mm,para);
       pos = sa.sa.grid_vtpm_4mm_indi;
    end
    
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
    elseif strcmp(sa.sa.leadfield, 'L_medium')
      pos = sa.sa.grid_medium_indi;
      pars.grid = 'medium';
      load ~/pconn/matlab/aalmask_grid_medium
    elseif strcmp(sa.sa.leadfield, 'aal_6mm')
      pars.grid = 'L_aal_6mm';
      pos = sa.sa.grid_aal6mm_indi;
    elseif strcmp(sa.sa.leadfield, 'aal_4mm')
      pars.grid = 'L_aal_4mm';
      pos = sa.sa.grid_aal4mm_indi;
    
%     elseif strcmp(para.grid,'m758_4mm')
%       pos	= sa.sa.grid_m758_4mm_indi;
%       m758  = tp_m758_grid();
%       sa.sa.aal_label = m758.tissue_4mm(m758.tissue_4mm>0);
%     elseif strcmp(para.grid,'m758_6mm')
% %       pos       = sa.sa.grid_m758_6mm_indi;
% %       filt      = pconn_beamformer(cs,sa.sa.L_m758_6mm,para);
% %       m758      = tp_m758_grid();
%       sa.sa.aal_label = m758.tissue_4mm(m758.tissue_4mm>0);
%     
    end
    filt      = get_spatfilt(pars);
  end
  
  % apply bandpass filter
  dloc1 =  filtfilt(bfilt,afilt,dloc1);
  
  % --------------------------------------
  % APPLY GAUSSIAN WEIGHTING (see Brookes et al., 2016)
  % --------------------------------------
  if strcmp(para.grid,'aal_4mm') || strcmp(para.grid,'aal_6mm')
    aal_mom = zeros(91,size(dloc1,1));
  elseif strcmp(para.grid,'m758_4mm') || strcmp(para.grid,'m758_6mm')
    aal_mom = zeros(758,size(dloc1,1));
  elseif strcmp(para.grid,'vtpm_4mm') 
    aal_mom = zeros(920,size(dloc1,1));
    tmp = tp_create_grid('vtpm');
    aalgrid.mask = tmp.label_4mm;

  end
  
  fprintf('Atlas distance weighting ...\n')
  for ireg = 1 : size(aal_mom,1)
    if ~isfield(sa.sa,'aal_label')
      aalgrid.mask = aalgrid.mask(find(aalgrid.mask));
      idx = find(aalgrid.mask==ireg);
    else
      idx = find(sa.sa.aal_label==ireg);
    end
      
    if length(idx)>1
      com = mean(pos(idx,:));
    else
      com = pos(idx,:);
    end

    dist = sqrt((com(1)-pos(idx,1)).^2 + (com(2)-pos(idx,2)).^2 + (com(3)-pos(idx,3)).^2);
    
    % In cm or mm?
    xx = pos; for i=1:3; xx(:,1)=xx(:,i)-mean(xx(:,i)); end
    xrad = mean(sqrt(sum(xx.^2,2)));
    
    if xrad > 20
      w    = exp((-dist.^2)./400);
    else
      w    = exp((-10*dist.^2)./400);
    end
    
    tmp = filt(:,idx)'*dloc1';
    % flip sign to account for arbitrary polarity
    tmp = repmat(adjustsign(tmp'),[1 size(tmp,2)]).*tmp;
    
    if para.weigh == 1
      aal_mom(ireg,:) = sum(repmat(w,[1 size(tmp,2)]).*tmp)./length(idx);
    elseif para.weigh == 2
      aal_mom(ireg,:) = sum(tmp);
    end
    
  end
  fprintf('AAL done ...\n')
  mom = aal_mom; clear aal_mom dloc1 filt tmp m758
  % --- END GAUSSIAN WEIGHTING ---
  
  % ---------------------------------
  % COMPUTE POWER CORRELATIONS
  % ---------------------------------
  if strcmp(para.wavelet,'bp_filt')
    
    for ireg = 1:size(mom,1)
      mom(ireg,:) = hilbert(mom(ireg,:));
    end
    
    % APPLY MOVING AVERAGE TO FILTER ENVELOPES
    % (not necessary)
    % --------------------------
    if para.scnd_filt    
      overlap = 0.9;
      mom = abs(mom);
      clear mom2
      winsize = 100;
      t = 1:size(mom,2);
      for vox= 1 : size(mom,1)
        [tmp,t_avg] = osl_movavg(mom(vox,:),t,winsize,overlap,0);
        mom2(vox,:) = tmp(~isnan(t_avg));
      end
      mom = mom2; clear mom2
      for ireg = 1:size(mom,1)
        mom(ireg,:) = hilbert(mom(ireg,:));
      end
      
      segleng  = fsample/400;
      segshift = 1;
      epleng   = size(mom,2);
      nseg     = floor((epleng-segleng)/segshift+1);
    end
    % END MOVING AVERAGE ----------
    
    fprintf('Computing orth. amp. correlations ...\n')
    
    if nep == 1
      if ~para.scnd_filt
        c =  compute_orthopowcorr(mom);
      end
    else
      % not implemented (future: prepare FCD here)
      error('There is no reason for doing this unless you wanna compute FCD')
      for iep = 1 : nep
          iep
          % DEFINE STUFF HERE
          tmp_mom = double(mom(:,(iep-1)*segshift+1:(iep-1)*segshift+segleng));
          c(:,:,iep) =  compute_orthopowcorr(tmp_mom);
        end
    end
  else
    
    % DO STUFF HERE WITH WAVELETS???
    % Nothing implemented
    % --------
    
  end
end

% This subfunction computes orthogonalized power correlation
% Code is copied from fieldtrip's ft_connectivity_powcorr_ortho.m
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


