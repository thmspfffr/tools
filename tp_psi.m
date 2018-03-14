function [psi] = tp_psi(data,para,sa,cs_psi)

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
    elseif isfield(sa.sa,'L_medium')
      pos = sa.sa.grid_medium_indi;
      pars.grid = 'medium';
      load ~/pconn/matlab/aalmask_grid_medium
    elseif isfield(sa.sa,'L_aal_6mm')
      pars.grid = 'L_aal_6mm';
      pos = sa.sa.grid_aal6mm_indi;
    end
    filt      = get_spatfilt(pars);
    filt      = filt(:,find(aalgrid.mask));
    pos       = pos(find(aalgrid.mask),:);
  end
  % ----------------------------------------
  
  
  
  fprintf('Computing psi ...\n')
  % YOU NEED SOURCE LEVEL CROSS SPECTRA!
  df=1;
  [nchan nchan nf]=size(cs_psi);
  
  cs_src = zeros(size(filt,2),size(filt,2),nf);
  for f = 1 : nf
    cs_src(:,:,f) = filt'*cs_psi(:,:,f)*filt;
  end
  
%   pp=cs_src;
  
  for f=1:nf
    
    pp(:,:,f)=single(cs_src(:,:,f)./sqrt(diag(cs_src(:,:,f))*diag(cs_src(:,:,f))'));
  end
  clear cs_src dloc1 cs_psi cs
  
  for i = 1 : size(pp,2)
   
    
    ps(i,:)=sum(imag(conj(pp(i,:,1:end-df)).*pp(i,:,1+df:end)),3);
    
  end
  
  
  fprintf('Atlas distance weighting ...\n')
  for ireg = 1 : 91
    ireg
    for jreg = 1 : 91
    %     fprintf('Atlas distance weighting ... reg%d / %d ... \n',ireg,size(aal_mom,1))
    if ~isfield(sa.sa,'aal_label')
      aalgrid.mask = aalgrid.mask(find(aalgrid.mask));
      idx = find(aalgrid.mask==ireg);
    else
      idx1 = find(sa.sa.aal_label==ireg);
      idx2 = find(sa.sa.aal_label==jreg);
    end
    
    
    psi(ireg,jreg) = mean(mean(ps(idx1,idx2)));
    end
  end

  end
  
  
  
