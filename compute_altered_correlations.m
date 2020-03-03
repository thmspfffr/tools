function emp = compute_altered_correlations(cleandat,para)



for ifoi = para.nfreq
  
  fprintf('Compute empirical results, freq%d...\n',ifoi)
  
  % rest data
  s_fc(:,:,:,:,1) = cleandat(:,:,:,:,1,ifoi);
  % task data
  s_fc(:,:,:,:,2) = cleandat(:,:,:,:,2,ifoi);
  
  % *_p_* = number of increased correlations
  % *_n_* = number of decreased correlations
  
  % -------------------------
  % ATOMOXETINE
  % -------------------------
  % during rest (condition label = 1)
  [h,~,~,s]=ttest(atanh(s_fc(:,:,:,2,1)),atanh(s_fc(:,:,:,1,1)),'dim',3,'alpha',para.alpha);
  emp.n_p_atx(ifoi,1) = nansum(nansum((h.*sign(s.tstat))>0))./(size(s_fc,1)*size(s_fc,1)-size(s_fc,1));
  emp.n_n_atx(ifoi,1) = nansum(nansum((h.*sign(s.tstat))<0))./(size(s_fc,1)*size(s_fc,1)-size(s_fc,1));
  % during task (condition label = 2)
  [h,~,~,s]=ttest(atanh(s_fc(:,:,:,2,2)),atanh(s_fc(:,:,:,1,2)),'dim',3,'alpha',para.alpha);
  emp.n_p_atx(ifoi,2) = nansum(nansum((h.*sign(s.tstat))>0))./(size(s_fc,1)*size(s_fc,1)-size(s_fc,1));
  emp.n_n_atx(ifoi,2) = nansum(nansum((h.*sign(s.tstat))<0))./(size(s_fc,1)*size(s_fc,1)-size(s_fc,1));
  % context dependence
  emp.n_p_context_atx(ifoi) = emp.n_p_atx(ifoi,1)-emp.n_p_atx(ifoi,2);
  emp.n_n_context_atx(ifoi) = emp.n_n_atx(ifoi,1)-emp.n_n_atx(ifoi,2);
%   [h,~,~,s]=ttest(atanh(s_fc(:,:,:,2,1))-atanh(s_fc(:,:,:,1,1)),atanh(s_fc(:,:,:,2,2))-atanh(s_fc(:,:,:,1,2)),'dim',3,'alpha',para.alpha);
%   emp.n_p_context_atx(ifoi,1) = nansum(nansum((h.*sign(s.tstat))>0))./(size(s_fc,1)*size(s_fc,1)-size(s_fc,1));
%   emp.n_n_context_atx(ifoi,1) = nansum(nansum((h.*sign(s.tstat))<0))./(size(s_fc,1)*size(s_fc,1)-size(s_fc,1));
  
  % ATOMOXETINE: number of altered correlations (irrespective of sign)
  % --------------------------
  % during rest (condition label = 1)
  h=ttest(s_fc(:,:,:,2,1),s_fc(:,:,:,1,1),'dim',3,'alpha',para.alpha);
  emp.atx(ifoi,1) = nansum(nansum((h)))./(size(s_fc,1)*size(s_fc,1)-size(s_fc,1));
  % during task (condition label = 1)
  h=ttest(s_fc(:,:,:,2,2),s_fc(:,:,:,1,2),'dim',3,'alpha',para.alpha);
  emp.atx(ifoi,2) = nansum(nansum((h)))./(size(s_fc,1)*size(s_fc,1)-size(s_fc,1));
  % context dependence
  emp.context_allconn_emp_atx(ifoi) = emp.atx(ifoi,1)-emp.atx(ifoi,2);
  
  % ATOMOXETINE: local changes
  % --------------------------
  % during rest (condition label = 1)
  [h,~,~,s]=ttest(atanh(s_fc(:,:,:,2,1)),atanh(s_fc(:,:,:,1,1)),'dim',3,'alpha',para.alpha);
  emp.n_p_atx_pervoxel(:,ifoi,1) = nansum((h.*sign(s.tstat))>0)./(size(s_fc,1)-1);
  emp.n_n_atx_pervoxel(:,ifoi,1) = nansum((h.*sign(s.tstat))<0)./(size(s_fc,1)-1);
  % during task (condition label = 2)
  [h,~,~,s]=ttest(atanh(s_fc(:,:,:,2,2)),atanh(s_fc(:,:,:,1,2)),'dim',3,'alpha',para.alpha);
  emp.n_p_atx_pervoxel(:,ifoi,2) = nansum((h.*sign(s.tstat))>0)./(size(s_fc,1)-1);
  emp.n_n_atx_pervoxel(:,ifoi,2) = nansum((h.*sign(s.tstat))<0)./(size(s_fc,1)-1);
  % context dependence: atx-pbo(rest) vs. atx-pbo(task)
  [h,~,~,s]=ttest(atanh(s_fc(:,:,:,2,1))-atanh(s_fc(:,:,:,1,1)),atanh(s_fc(:,:,:,2,2))-atanh(s_fc(:,:,:,1,2)),'dim',3,'alpha',para.alpha);
  emp.n_p_context_atx_pervoxel(:,ifoi,1) = nansum((h.*sign(s.tstat))>0)./(size(s_fc,1)-1);
  emp.n_n_context_atx_pervoxel(:,ifoi,1) = nansum((h.*sign(s.tstat))<0)./(size(s_fc,1)-1);
  
  % --------------------------
  % DONEPEZIL
  % --------------------------
  % during rest (condition label = 1)
  [h,~,~,s]=ttest(atanh(s_fc(:,:,:,3,1)),atanh(s_fc(:,:,:,1,1)),'dim',3,'alpha',para.alpha);
  emp.n_p_dpz(ifoi,1) = nansum(nansum((h.*sign(s.tstat))>0))./(size(s_fc,1)*size(s_fc,1)-size(s_fc,1));
  emp.n_n_dpz(ifoi,1) = nansum(nansum((h.*sign(s.tstat))<0))./(size(s_fc,1)*size(s_fc,1)-size(s_fc,1));
  % during task (condition label = 2)
  [h,~,~,s]=ttest(atanh(s_fc(:,:,:,3,2)),atanh(s_fc(:,:,:,1,2)),'dim',3,'alpha',para.alpha);
  emp.n_p_dpz(ifoi,2) = nansum(nansum((h.*sign(s.tstat))>0))./(size(s_fc,1)*size(s_fc,1)-size(s_fc,1));
  emp.n_n_dpz(ifoi,2) = nansum(nansum((h.*sign(s.tstat))<0))./(size(s_fc,1)*size(s_fc,1)-size(s_fc,1));
  % context dependence: dpz-pbo(rest) vs. dpz-pbo(task)
  emp.n_p_context_dpz(ifoi) = emp.n_p_dpz(ifoi,1)-emp.n_p_dpz(ifoi,2);
  emp.n_n_context_dpz(ifoi) = emp.n_n_dpz(ifoi,1)-emp.n_n_dpz(ifoi,2);
%   [h,~,~,s]=ttest(atanh(s_fc(:,:,:,3,1))-atanh(s_fc(:,:,:,1,1)),atanh(s_fc(:,:,:,3,2))-atanh(s_fc(:,:,:,1,2)),'dim',3,'alpha',para.alpha);
%   emp.n_p_context_dpz(ifoi,1) = nansum(nansum((h.*sign(s.tstat))>0))./(size(s_fc,1)*size(s_fc,1)-size(s_fc,1));
%   emp.n_n_context_dpz(ifoi,1) = nansum(nansum((h.*sign(s.tstat))<0))./(size(s_fc,1)*size(s_fc,1)-size(s_fc,1));
  
  % DONEPEZIL: number of altered correlations (irrespective of sign)
  % --------------------------
  % during rest (condition label = 1)
  h=ttest(atanh(s_fc(:,:,:,3,1)),atanh(s_fc(:,:,:,1,1)),'dim',3,'alpha',para.alpha);
  emp.dpz(ifoi,1) = nansum(nansum((h)))./(size(s_fc,1)*size(s_fc,1)-size(s_fc,1));
  % during task (condition label = 2)
  h=ttest(atanh(s_fc(:,:,:,3,2)),atanh(s_fc(:,:,:,1,2)),'dim',3,'alpha',para.alpha);
  emp.dpz(ifoi,2) = nansum(nansum((h)))./(size(s_fc,1)*size(s_fc,1)-size(s_fc,1));
  % context dependence
  emp.context_allconn_emp_dpz(ifoi) = emp.dpz(ifoi,1)-emp.dpz(ifoi,2);
  
  % DONEPEZIL: local changes
  % --------------------------
  % during rest (condition label = 1)
  [h,~,~,s]=ttest(atanh(s_fc(:,:,:,3,1)),atanh(s_fc(:,:,:,1,1)),'dim',3,'alpha',para.alpha);
  emp.n_p_dpz_pervoxel(:,ifoi,1) = nansum((h.*sign(s.tstat))>0)./(size(s_fc,1)-1);
  emp.n_n_dpz_pervoxel(:,ifoi,1) = nansum((h.*sign(s.tstat))<0)./(size(s_fc,1)-1);
  % during task (condition label = 2)
  [h,~,~,s]=ttest(atanh(s_fc(:,:,:,3,2)),atanh(s_fc(:,:,:,1,2)),'dim',3,'alpha',para.alpha);
  emp.n_p_dpz_pervoxel(:,ifoi,2) = nansum((h.*sign(s.tstat))>0)./(size(s_fc,1)-1);
  emp.n_n_dpz_pervoxel(:,ifoi,2) = nansum((h.*sign(s.tstat))<0)./(size(s_fc,1)-1);
  % context dependence: atx-pbo(rest) vs. atx-pbo(task)
  [h,~,~,s]=ttest(atanh(s_fc(:,:,:,3,1))-atanh(s_fc(:,:,:,1,1)),atanh(s_fc(:,:,:,3,2))-atanh(s_fc(:,:,:,1,2)),'dim',3,'alpha',para.alpha);
  emp.n_p_context_dpz_pervoxel(:,ifoi,1) = nansum((h.*sign(s.tstat))>0)./(size(s_fc,1)-1);
  emp.n_n_context_dpz_pervoxel(:,ifoi,1) = nansum((h.*sign(s.tstat))<0)./(size(s_fc,1)-1);
  
  % --------------------------
  % TASK VS REST (during placebo only)
  % --------------------------
  % global changes
  [t_tvsr,~,~,s] = ttest(s_fc(:,:,:,1,2),s_fc(:,:,:,1,1),'dim',3,'alpha',para.alpha);
  t_tvsr = t_tvsr.*sign(s.tstat); clear s
  emp.taskvsrest_p(ifoi) = nansum(nansum(t_tvsr>0))./(size(s_fc,1)*size(s_fc,1)-size(s_fc,1));
  emp.taskvsrest_n(ifoi) = nansum(nansum(t_tvsr<0))./(size(s_fc,1)*size(s_fc,1)-size(s_fc,1));
  % local changes
  emp.taskvsrest_p_pervoxel(:,ifoi) = nansum(t_tvsr>0)./(size(s_fc,1)-1);
  emp.taskvsrest_n_pervoxel(:,ifoi) = nansum(t_tvsr<0)./(size(s_fc,1)-1);
  
end

% --------------------------
% DOUBLE DISSOCIATION ACROSS FREUQUENCIES
% --------------------------
% global effects
emp.doubledissociation_emp = emp.context_allconn_emp_atx-emp.context_allconn_emp_dpz;
% --------------------------
