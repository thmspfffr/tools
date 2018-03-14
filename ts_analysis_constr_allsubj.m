clear;

% ------------------------------------------------------
% VERSION 8
% ------------------------------------------------------
% Non-constrained fit, multiple starting parameters, 
% iterative adjustment of fitted parameters
% see email from rani moran (08/02/13)
% v = 8;
% ------------------------------------------------------
% VERSION 9
% ------------------------------------------------------
% Constrained fit, multiple starting parameters, 
% iterative adjustment of fitted parameters
v = 9;
% ------------------------------------------------------

addpath(genpath('/home/tpfeffer/timescale/'))

warning off

NCOH = 6;
nsubj=9;
remove_trials=100;

% iterations for bootstrap
kboot=10000; % just for overall (over coh-lvls) performance
nboot = 2; % NUMBER OF BOOT ITERATIONS: n_boot_idx * nboot
n_boot_idx = 10;

if nboot >= 1
    boot = 1;
else
    boot = 0;
end

weibull = @(beta,x,CL,pm) CL+(1-2*CL-pm)*(1-exp(-(x/beta(1)).^beta(2)));
direc = '/Users/tpfeffer/Documents/projects/timescale/data_ori/ori_curr_biol_271012/';
direc = '/home/tpfeffer/timescale/projects/timescale/data_ori/ori_curr_biol_271012/';

% 
% if ~exist([direc sprintf('out_ana/%d/',v)],'dir')
%   mkdir(sprintf([direc 'out_ana/%d/'],v))
% end

% OPTIONS FOR FITTING ROUTINE
% ------------------------------------------------------
options.MaxFunEvals = 500000000;
options.MaxIter = 50000000;
options.TolX = 0.00000001;
options.TolFun = 0.00000001;
options = optimset('Display', 'off') ;
options.Robust = 'on';
% ------------------------------------------------------
% LOOP OVER ALL SUBJECTS
% ------------------------------------------------------
for isubj=1:nsubj
  
  load([direc sprintf('data_subj%i',isubj)])
  DUR  = [90 60 30 15];

  clear coh
  % ------------------------------------------------------
  % MLE with bootstrap
  % Maximum-likelihood parameter estimation as described by Wichmann, FA &
  % Hill, NJ (2001) and Mjung, IJ (2003). The log likelihood is minimized
  % with the Nelder-Mead algorithm (fminsearch).
  % ------------------------------------------------------
  for icond = 1 : 2

    if exist(sprintf([direc 'out_ana/%d/workspace_subj%d_cond%d_constr_v%d_processing.mat'],v,isubj,icond,v));
      continue
    else
      fclose(fopen(sprintf([direc 'out_ana/%d/workspace_subj%d_cond%d_constr_v%d_processing.mat'],v,isubj,icond,v),'w'));
    end

    data{icond}(1:remove_trials,:)=[];

    psych{icond}.coh =  [0, 0.08, 0.16,0.24,0.32,0.4];
    coh = [psych{icond}.coh];

    boot_thresh{icond}=zeros(nboot,length(DUR));    % preallocation
    boot_slope{icond}=zeros(nboot,length(DUR));     % preallocation

    % correct rejects
    psych{icond}.cr = length(find(data{icond}(:,1) == 1 & data{icond}(:,4) == 1))/length(find(data{icond}(:,1) == 1));
    % false alarms: 1 - cr
    psych{icond}.fa = 1-psych{icond}.cr;
    psych{icond}.dur = DUR;

    for iboot=1:kboot
      nTrl = size(data{icond},1);
      rind = ceil(rand(1,nTrl)*nTrl);
      rand_data = data{icond}(rind,:);
      boot_cr(icond,iboot) = length(find(rand_data(:,1) == 1 & rand_data(:,4) == 1))/length(find(rand_data(:,1) == 1));
      boot_fa(icond,iboot) = 1-boot_cr(icond,iboot);
    end
    % ------------------------------------------------------

    for idur = 1:length(DUR)

      percCorAllCoh{icond}(idur)=length(data{icond}(data{icond}(:,4)==1&data{icond}(:,1)~=1&data{icond}(:,2)==DUR(idur)))/length(data{icond}(data{icond}(:,2)==DUR(idur)&data{icond}(:,1)~=1));
      % ------------------------------------------------------
      % nonparametric bootstrap to obtain CIs
      % ------------------------------------------------------
      if boot
        for iboot=1:kboot
          tmp_data = data{icond}(data{icond}(:,2)==DUR(idur) & data{icond}(:,4)~=2,:);
          nTrl = size(tmp_data,1);
          rind = ceil(rand(1,nTrl)*nTrl);
          rand_data = tmp_data(rind,:);
          boot_percCorAllCoh{icond}(iboot,idur)=length(rand_data(rand_data(:,4)==1&rand_data(:,1)~=1&rand_data(:,2)==DUR(idur)))/length(rand_data(rand_data(:,2)==DUR(idur)&rand_data(:,1)~=1));
        end
      end
      % ------------------------------------------------------
      % select data for maximum-likelihood estimation
      % ------------------------------------------------------
      mle_data = data{icond}(data{icond}(:,1)~=1 & data{icond}(:,2) == DUR(idur),:);
      nTrl = size(mle_data,1);
      tmp_coh = [min(mle_data(:,1)) : max(mle_data(:,1))];
      coh_lvls{icond}(idur,:) = coh(tmp_coh);

      for icoh = 1 : length(coh_lvls{icond}(idur,:))
          select_coh = tmp_coh(icoh);
          psych{icond}.RT(idur,icoh)=mean(mle_data(mle_data(:,3)>200&mle_data(:,3)<1000&mle_data(:,1)==select_coh,3));
          nTrials(idur,icoh) = length(mle_data(mle_data(:,1)==select_coh));
          prc_cor(idur,icoh) = length(mle_data(mle_data(:,1) == select_coh & mle_data(:,4) == 1))/length(mle_data(mle_data(:,1)==select_coh));
          binom(idur,icoh) = nchoosek(nTrials(idur,icoh),round((prc_cor(idur,icoh)*nTrials(idur,icoh))));
          psych{icond}.dprime(idur,icoh) = norminv(prc_cor(idur,icoh))-norminv(psych{icond}.fa); % z(H)-z(FA)
          psych{icond}.propcormax(idur,icoh) = normcdf(0.5*(norminv(prc_cor(idur,icoh))-norminv(psych{icond}.fa)));
      end
      psych{icond}.hr(idur,:) = prc_cor(idur,:);

      % CHANCE-LEVEL ESTIMATION ACCORDING TO EMAIL (JAN,7TH)
      CL{icond}(idur) = 0.5*psych{icond}.fa*(DUR(idur)*10+600)/5000;
      pm{icond}       = 0.4*psych{icond}.fa;
    end
    % ---------------------------------------------------------
    % MAXIMUM LIKELIHOOD PARAMETER ESTIMATION
    % as described by wichmann & hill (2001)
    % fit all durations with a familiy of weibull functions
    % (see gao et al., 2011)
    % ---------------------------------------------------------
    nits = 2000;
    
    beta=zeros(4,2,nits);
    fval=zeros(nits);
    
    for i = 1 : nits
      beta1 = [sort(linspace(.5*rand,.5*rand,4))' 5*rand(4,1)];
      while ~issorted(beta1(:,1)) || (~issorted(beta1(:,2)) && ~issorted(-beta1(:,2))) || (beta1(4,2)/beta1(1,2)>2 && beta1(1,2)/beta1(4,2)>2)
        beta1 = [rand randi(100,1,1)/10; rand randi(100,1,1)/10;rand randi(100,1,1)/10;rand randi(100,1,1)/10];    
      end
      disp(sprintf('%d-%d-%d',icond,idur,i)) 
      [beta(:,:,i) fval(i)] = fminsearchbnd(@(beta1) fun_constr(beta1,CL{icond},pm{icond},binom,prc_cor,coh_lvls{icond}(idur,:),nTrials),beta1,[0.01,0;0.01,0;0.01,0;0.01,0],[1,20;1,20;1,20;1,20]);
      while any(any(beta1 - beta(:,:,i) ~= 0))
        disp('ui')
        beta1(:,1) = beta(:,1,i)';
        beta1(:,2) = beta(:,2,i)';
        [beta(:,:,i) fval(i)] = fminsearchbnd(@(beta1) fun_constr(beta1,CL{icond},pm{icond},binom,prc_cor,coh_lvls{icond}(idur,:),nTrials),beta1,[0.01,0;0.01,0;0.01,0;0.01,0],[1,20;1,20;1,20;1,20]);
%       [beta(:,:,i) fval(i)]=simulannealbnd(@(beta1) fun_constr(beta1,CL{icond},pm{icond},binom,prc_cor,coh_lvls{icond}(idur,:),nTrials),beta1,[0.01,0;0.01,0;0.01,0;0.01,0],[1,20;1,20;1,20;1,20])
      end
%       end
    end
    
    [~,idx]=min(fval);
    
    psych{icond}.thresh = beta(:,1,idx)';
    psych{icond}.slope  = beta(:,2,idx)';
    
    clear beta fval
    save(sprintf([direc 'out_ana/%d/workspace_subj%d_cond%d_constr_v%d.mat'],v,isubj,icond,v));

  end
end
