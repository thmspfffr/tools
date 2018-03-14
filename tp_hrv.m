%% tp_hrv.m
% pconn_hrv

% last update: 01-03-2015, tpfeffer
% implemented: trigger resample

clear
restoredefaultpath

% --------------------------------------------------------
% VERSION 1
% --------------------------------------------------------
v         = 1;
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
i_fit     = [3 50];
% --------------------------------------------------------

outdir = '/home/tpfeffer/pconn/proc/dfa/';
indir  = '/home/tpfeffer/pconn/proc/preproc/';

tp_addpaths


%%

for isubj = SUBJLIST
  for m = 1 : 3
    for ibl = 1 : 2
            
      if ~exist(sprintf([outdir 'pconn_hrv_s%d_m%d_b%d_v%d_processing.txt'],isubj,m,ibl,v))
        system(['touch ' outdir sprintf('pconn_hrv_s%d_m%d_b%d_v%d_processing.txt',isubj,m,ibl,v)]);
      else
        continue
      end

      disp(sprintf('Processing s%d b%d m%d',isubj,ibl,m));
      
        
        load([indir sprintf('pconn_preproc_data_s%d_m%d_b%d_v%d.mat',isubj,m,ibl,2)])

        cfg = [];
        cfg = cfgs;

        cfg.channel     = {'EEG059'};

        cfg.padding     = 10;
        cfg.continuous  = 'yes';

        cfg.bsfilter    = 'yes';
        cfg.lpfilter    = 'yes';
        cfg.hpfilter    = 'yes';
        
        cfg.hpfiltord = 2;

        cfg.bsfreq      = [49 51; 99 101; 149 151; 199 201];
        cfg.lpfreq      = 30;
        cfg.hpfreq      = 0.7;

        data.Fs         = 1200;
        data            = ft_preprocessing(cfg);

        cfg = [];
        cfg.detrend = 'yes';
        cfg.demean  = 'yes';
        cfg.resamplefs = 400;
        data = ft_resampledata(cfg,data);
        
        br = 0;
        th = 20:-0.1:0;
        cnt = 0;
        th_idx=[];
        hb_all=[];
        
        for ith = 1 : length(th)
          
          [~,l]=findpeaks(abs(zscore(data.trial{1}(3600:end-3600))),'MinPeakHeight',th(ith),'MinPeakDist',100);
          
          if isempty(l)   
            th_idx(ith)=0;
            hb_all(ith) = 0;
            continue
          else   
            hb = size(l,2)/((size(data.trial{1}(3600:end-3600),2)/400)/60);

            if hb>100 || hb<45
              hb_all(ith) = 0;
              th_idx(ith) = 0;
            else
              hb_all(ith) = hb;
              th_idx(ith) = 1;
            end     
          end
        end
        
        th_idx_cp = th_idx;
          
        th = 20:-0.1:0;
          
        par.th = mean(th(logical(th_idx)));        
          
        [~,l]=findpeaks(abs(zscore(data.trial{1}(3600:end-3600))),'MinPeakHeight',par.th,'MinPeakDist',100);
        
        par.hb = size(l,2)/((size(data.trial{1},2)/400)/60);
        
        % SAVE OUTPUT IF HB > 100 or < 45
        if par.hb>100 || par.hb<45
          save([outdir sprintf('pconn_hrv_dfa_s%d_m%d_b%d_v%d_slowfastweird.mat',isubj,m,ibl,v)],'hb','-v7.3');
          continue
        end

        % COMPUTE DFA
        dfa     = tp_dfa(diff(l)',i_fit,1,0.5,15);
        par.dfa = dfa.exp; 
        par.var = var(diff(l)); clear dfa l
       
        save([outdir sprintf('pconn_hrv_dfa_s%d_m%d_b%d_v%d.mat',isubj,m,ibl,v)],'par','-v7.3');

      clear data
    end
  end
end
% error('Done')
%%
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
v = 1
ord = pconn_randomization;
% ord = ord(SUBJLIST,:);

for isubj = SUBJLIST
  
  for m = 1 : 3
    for ibl = 1 : 2

    try
      im = find(ord(isubj,:)==m);
      load([outdir sprintf('pconn_hrv_dfa_s%d_m%d_b%d_v%d.mat',isubj,im,ibl,v)]);
      
      dfa(isubj,ibl,m) = par.dfa;
      hb_all(isubj,ibl,m) = par.hb;
% var(isubj,ibl,m) = par.var;
  
    catch me
    end
    
    end
  end
end

dfa     = dfa(SUBJLIST,:,:);
hb_all  = hb_all(SUBJLIST,:,:);
% var = var(SUBJLIST,:,:);
dfa(dfa==0)       = NaN;
hb_all(hb_all==0) = NaN;
% var(var==0)=nan;
dfa     = squeeze(nanmean(dfa,2));
hb_all  = squeeze(nanmean(hb_all,2));
% var = squeeze(nanmean(var,2));
%% PLOT HEART RATE

figure; set(gcf,'color','white'); hold on;
title('HEART RATE')

bar(1.5,nanmean(hb_all(:,1)),0.3,'FaceColor','r','EdgeColor','r'); ylabel('DFA-Exponent'); % ylim([0.65 0.85]);
bar(2,nanmean(hb_all(:,2)),0.3,'FaceColor','b','EdgeColor','b'); ylabel('DFA-Exponent');   %ylim([0.65 0.85])
bar(2.5,nanmean(hb_all(:,3)),0.3,'FaceColor','m','EdgeColor','m'); ylabel('DFA-Exponent');   %ylim([0.65 0.85])

lm = nanmean(hb_all)-[nanstd(hb_all)/sqrt(size(hb_all,1))]; 
um = nanmean(hb_all)+[nanstd(hb_all)/sqrt(size(hb_all,1))]; 

line([1.5 1.5],[lm(1) um(1)],'color',[0.6 0 0],'LineWidth',5);
line([2 2],[lm(2) um(2)],'color',[0 0 0.6],'LineWidth',5);
line([2.5 2.5],[lm(3) um(3)],'color',[0.6 0.1 0.6],'LineWidth',5);

ylabel('Heart rate (bpm)'); xlabel('Condition');
axis([1.3 2.8 50 80]);
print(gcf,'-depsc2',sprintf('~/pconn/proc/plots/pconn_hrv_hb_v%d.eps',v))

[~,p1] = ttest(hb_all(:,1),hb_all(:,2))
[~,p2] = ttest(hb_all(:,2),hb_all(:,3))
[~,p3] = ttest(hb_all(:,1),hb_all(:,3))




%% PLOT DFA

figure; set(gcf,'color','white'); hold on;
title('HRV DFA')

bar(1.5,nanmean(dfa(:,1)),0.3,'FaceColor','r','EdgeColor','r'); ylabel('DFA-Exponent'); % ylim([0.65 0.85]);
bar(2,nanmean(dfa(:,2)),0.3,'FaceColor','b','EdgeColor','b'); ylabel('DFA-Exponent');   %ylim([0.65 0.85])
bar(2.5,nanmean(dfa(:,3)),0.3,'FaceColor','m','EdgeColor','m'); ylabel('DFA-Exponent');   %ylim([0.65 0.85])

lm = nanmean(dfa)-[nanstd(dfa)/sqrt(size(dfa,1))]; 
um = nanmean(dfa)+[nanstd(dfa)/sqrt(size(dfa,1))]; 

line([1.5 1.5],[lm(1) um(1)],'color',[0.6 0 0],'LineWidth',5);
line([2 2],[lm(2) um(2)],'color',[0 0 0.6],'LineWidth',5);
line([2.5 2.5],[lm(3) um(3)],'color',[0.6 0.1 0.6],'LineWidth',5);

ylabel('DFA Exponent'); xlabel('Condition');

print(gcf,'-depsc2',sprintf('~/pconn/proc/plots/pconn_hrv_dfa_v%d.eps',v))

[~,p_hrv1] = ttest(dfa(:,1),dfa(:,2))
[~,p_hrv2] = ttest(dfa(:,2),dfa(:,3))
[~,p_hrv3] = ttest(dfa(:,1),dfa(:,3))


%%
if any(isnan(dfa(:)))
  [nan_row,~]=find(isnan(dfa));
end



ifoi = 3;
v = 2;
clear r p
load /home/tpfeffer/pconn/proc/src/pconn_sa_s4_m1_b1_v1.mat
% ord = pconn_randomization;
%   v = 2;
clear dfa_all

for m = 1 : 3
  
  for isubj = SUBJLIST
    im = find(ord(isubj,:)==m);
    
    load(sprintf([outdir 'pconn_sens_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v));
    
    dfa_all(:,:,isubj,m) = par.dfa; clear par
    
  end
end

load(sprintf('~/pconn/proc/dfa/pconn_sens_dfa_stat_clusterstat_f%d_contr%d_v%d.mat',3,1,2)); 
senssel = find(stats.mask);

dfa_all      = dfa_all(:,:,SUBJLIST,:);
dfa_all_diff = squeeze(nanmean(dfa_all(senssel,:,:,2),2)-nanmean(dfa_all(senssel,:,:,1),2));
% dfa_all      = squeeze(nanmean(nanmean(dfa_all,2),4));
if exist('nan_row','var')
  dfa(nan_row,:)=[];
  dfa_all(:,nan_row) = [];
  dfa_all_diff(:,nan_row) = [];
end


% d = squeeze(nanmean(dfa_all(senssel,:),1));
d = squeeze(nanmean(dfa_all_diff,1));
dfa_diff = dfa(:,2)-dfa(:,1);

figure; set(gcf,'color','white'); hold on

slp = pconn_regress(d,dfa_diff);

scatter(d,dfa_diff,250,'r','FaceColor','r');
l=[min(d) max(d)].*slp(2)+slp(1);
% line([min(d) max(d)],[l(1) l(2)],'color','r','LineWidth',5);

% slp = pconn_regress(d,dfa(:,2));
% 
% scatter(d,dfa(:,2),250,'b','FaceColor','b');
% l=[min(d) max(d)].*slp(2)+slp(1);
% line([min(d) max(d)],[l(1) l(2)],'color','b','LineWidth',5);
% 
% slp = pconn_regress(d,dfa(:,3));
% 
% scatter(d,dfa(:,3),250,'m','FaceColor','m');
% l=[min(d) max(d)].*slp(2)+slp(1);
% line([min(d) max(d)],[l(1) l(2)],'color','m','LineWidth',5);
% 
% axis tight; axis([min(d)-0.02 max(d)+0.02 0.45 1.25]);
% 
[r,p]=corrcoef(d,dfa_diff)
% [r{2},p{2}]=corrcoef(d,dfa(:,2));
% [r{3},p{3}]=corrcoef(d,dfa(:,3));

% xlabel('DFA Exponent (MEG)'); ylabel('DFA Exponent (HRV)');
% title(sprintf('DFA MEG-HB correlation, f%d, r(1) = %.2f, p(1) = %.2f, r(2) = %.2f, p(2) = %.2f, r(3) = %.2f, p(3) = %.2f',ifoi,r{1}(1,2),p{1}(1,2),r{2}(1,2),p{2}(1,2),r{3}(1,2),p{3}(1,2)));

set(gca,'TickDir','out')


print(gcf,'-depsc2',sprintf('~/pconn/proc/plots/pconn_sens_hrvmegcorr_f%d_v%d.eps',ifoi,v))

%% HRV DFA CORR WITH TOPO
ifoi = 3;
v = 2;
clear r p
  load /home/tpfeffer/pconn/proc/src/pconn_sa_s4_m1_b1_v1.mat

%   v = 2;
clear dfa_all

for m = 1 : 3
  
  for isubj = SUBJLIST
    im = find(ord(isubj,:)==m);
    
    load(sprintf([outdir 'pconn_sens_dfa_s%d_m%d_f%d_v%d.mat'],isubj,im,ifoi,v));
    
    dfa_all(:,:,isubj,m) = par.dfa; clear par
    
  end
end

dfa_all = dfa_all(:,:,SUBJLIST,:);
dfa_all = squeeze(nanmean(nanmean(dfa_all,2),4));

for ichan = 1 : size(dfa_all,1)
  [tmp p_tmp] = corrcoef(nanmean(d,2),dfa_all(ichan,:));
  r(ichan) = tmp(1,2);
  p(ichan) = p_tmp(1,2);
end

g=figure; set(g,'color','white');

pars            = [];
pars.cbar       = 0;
pars.scale      = [-1 1];
pars.markersize = 0;
pars.linewidth  = 9;
pars.resolution = 300;

showfield_colormap((p<0.05).*r,sa.locs_2D,pars);

g=figure; set(g,'color','white');

pars            = [];
pars.cbar       = 0;
pars.scale      = [0 0.5];
pars.markersize = 0;
pars.linewidth  = 9;
pars.resolution = 300;

showfield_colormap((p<0.05).*p,sa.locs_2D,pars);