function [dfa] = tp_dfa_src_weighted(data,para,sa)

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
    end
    
  elseif strcmp(para.filt,'eloreta')
    
%     pars      = [];
%     pars.filt = 'eloreta';
%     pars.cs   = cs;
%     pars.foi  = f;
%     pars.sa   = sa.sa;
%     if isfield(sa.sa,'L_coarse')
%       pars.grid = 'cortex';
%       pos = sa.sa.grid_cortex3000_indi;
%       load('~/Documents/MATLAB/aalmask_grid_cortex3000.mat')
%     elseif  isfield(sa.sa,'L_aal')
%       pars.grid = 'aal';
%     elseif isfield(sa.sa,'L_medium')
%       pos = sa.sa.grid_medium_indi;
%       pars.grid = 'medium';
%       load ~/pconn/matlab/aalmask_grid_medium
%     elseif isfield(sa.sa,'L_aal_6mm')
%       pars.grid = 'L_aal_6mm';
%       pos = sa.sa.grid_aal6mm_indi;
%     end
%     filt      = get_spatfilt(pars);
%     filt      = filt(:,find(aalgrid.mask));
%     pos       = pos(find(aalgrid.mask),:);
  end
  % ----------------------------------------
  % bandpass filter data before?
  dloc1 =  filtfilt(bfilt,afilt,dloc1);
  % --------------------------------------
  % COMPUTE GAUSSIAN WEIGHTING (see Brookes et al., 2016)
  % --------------------------------------
  if strcmp(para.grid,'aal_4mm') || strcmp(para.grid,'aal_6mm')
    aal_mom = zeros(91,size(dloc1,1));
  elseif strcmp(para.grid,'m758_4mm') || strcmp(para.grid,'m758_6mm')
    aal_mom = zeros(758,size(dloc1,1));
  end
  
  fprintf('Atlas distance weighting ...\n')
  for ireg = 1 : size(aal_mom,1)
%     fprintf('Atlas distance weighting ... reg%d / %d ... \n',ireg,size(aal_mom,1))
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
    
  % COMPUTE CORRELATIONS BASED ON BAND-PASS FILTERED SIGNAL
  if strcmp(para.wavelet,'bp_filt')
    
    % make sure hilbert transform is applied in the right position
    for ireg = 1:size(mom,1)
      mom(ireg,:) = abs(hilbert(mom(ireg,:)));
      
    end
    
  end
      
  dfa(iep,:) = tp_dfa(mom',[3 50],400,0.5,15);
 
end



    
    
    
    
    
    
    
    