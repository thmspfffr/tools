function csd = tp_wavelet_crossspec(data,freqoi,width)

data.dat = data.trial;
data.avg = data.dat; data.trial = [];
data.dimord = 'chan_time';
data.time = [];
data.time =  1/data.fsample:1/data.fsample:size(data.avg,2)/data.fsample;

%%% 1.3 loop over frequencies & source reconstruct
srate    = 400;

% load(['/home/tpfeffer/pp/proc/src/' sprintf('pp_cnt_sa_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,1)],'sa');
% clear nanidx
  
  cfg=[];
  cfg.method  ='wavelet';
  cfg.output  = 'fourier';
  cfg.channel = {'MEG'};
  cfg.foi     = freqoi;
  cfg.width   = width; % again, as per Hipp et al. (2012) Nat Neurosci
  tempSD      = 1./(2*pi*(cfg.foi./cfg.width)); % temporal SD in sec
  tempWd      = round(3*tempSD*srate)/srate; % set to 3 std dev, comp with 1000 Hz sampl rate
  cfg.toi     = tempWd.*(1:floor(data.time(end)./tempWd));
  
  cfg.pad='nextpow2';
  tf=ft_freqanalysis(cfg,data);
  tf.fourierspctrm(:,:,:,end)=[];
  
  nanidx = isnan(squeeze(tf.fourierspctrm(:,1,:,:)));
  
  %         % Compute mean pupil diameter for time segments
  %         clear pup pow
  %         for itoi = 1 : size(cfg.toi,2)-1
  %           t = cfg.toi(itoi);
  %           idx = find(round(data.time.*10000)==round(t*10000));
  %           pup(itoi,:) = mean(pupil(idx-(tempWd*srate-1)/2:idx+(tempWd*srate-1)/2,:));
  %         end
  
  % ---------------
  tf.time(end) = [];
  tf.time(nanidx) = [];
  tf.fourierspctrm(:,:,:,nanidx) = [];
  %         pup(nanidx) = [];
  
  csd=zeros(numel(tf.label)*[1 1]);
  % csd calc excludes artifact bins
  csdTime=tf.time(:);
  csdData=tf.fourierspctrm(:,:,:,1:end-1);
  for itbin=1:numel(csdTime)-1
    fspec=squeeze(csdData(:,:,:,itbin)).';
    for ichan=1:numel(tf.label)
      csd(:,ichan)=csd(:,ichan)+fspec(ichan)*conj(fspec);
    end
  end
  csd=csd./(numel(tf.time)-1); % avg cross-spectral dens matrix
  csdData=[]; csdTime=[];
  
%   para      = [];
%   v         = 1;
%   para.iscs = 1;
%   para.reg  = 0.05;
%   filt      = tp_beamformer(csd,sa.L_BNA_5mm,para);
% end