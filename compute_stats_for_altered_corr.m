function outp = compute_stats_for_altered_corr(fc,para)
%COMPUTE_STATS_FOR_ALTERED_CORR 
%   outp = compute_stats_for_altered_corr(data,para) computes p-values for
%   fraction of altered correlations by means of a (unpaired) permutation
%   test. 
%   ----------------
%   INPUTS
%   ----------------
%   fc: Fuctional connectivity, with size NVOX,NVOX,NSUBJ,NCOND,NFREQ
%   - para.nfreq: Frequency indicies altered correlations should be computed
%   over (e.g. only frequencies with index [2 3 4])
%   - para.alpha: alpha-level for first-level ttest (default: 0.05)
%   - para.collect_result: 'yes' - output p-values, 'no' - output permutations
%   - para.nperm: number of permutations (default 10000)
%   - para.subs: only relevant for parallelization (default = para.nperm)

if ~isfield(para,'nperm')
  para.nperm = 10000;
end
if ~isfield(para,'collect_result')
  para.collect_result = 'yes';
end
if ~isfield(para,'alpha')
  para.alpha = 0.05;
end
if ~isfield(para,'subs')
  para.subs = para.nperm;
end
alp = para.alpha;
para.subs = 500;
para.allperms = para.nperm/para.subs;

nvox = size(fc,1)*size(fc,1)-size(fc,1);

%%
tmp = clock;
seed = ((tmp(1)+tmp(2)*tmp(3))/tmp(4)+tmp(5))*tmp(6);
rng(seed,'twister')

all_idx1 = randi(2,[size(fc,3),para.nperm]);

% if ~exist(sprintf([outdir 'pupmod_src_powcorr_permtest_perms_subs%d_nperm%d_v%d.mat'],par.subs,nperm,v))
%   all_idx1 = randi(2,[size(SUBJLIST,2),nperm]);
%   save(sprintf([outdir 'pupmod_src_powcorr_permtest_perms_subs%d_nperm%d_v%d.mat'],par.subs,nperm,v),'all_idx1');
% else
%   load(sprintf([outdir 'pupmod_src_powcorr_permtest_perms_subs%d_nperm%d_v%d.mat'],par.subs,nperm,v));
% end

fc = single(fc);

for iperm = 1 : para.allperms
  
  %   if ~exist(sprintf([outdir 'pupmod_src_powcorr_permtest_iperm%d_nperm%d_v%d_processing.txt'],iperm,nperm,v))
  %     system(['touch ' outdir sprintf('pupmod_src_powcorr_permtest_iperm%d_nperm%d_v%d_processing.txt',iperm,nperm,v)]);
  %   else
  %     continue
  %   end
  
  
  perm.tperm_n    = zeros(para.subs,length(para.nfreq));
  perm.tperm_p    = zeros(para.subs,length(para.nfreq));
  perm.tperm_all  = zeros(para.subs,length(para.nfreq));
  perm.tperm_pervoxel_neg = zeros(size(fc,1),para.subs,length(para.nfreq));
  perm.tperm_pervoxel_pos = zeros(size(fc,1),para.subs,length(para.nfreq));
  perm.tperm_pervoxel_all = zeros(size(fc,1),para.subs,length(para.nfreq));
    
  for kperm = 1 : para.subs 
    
    fprintf('Perm #%d\n',kperm);

    if strcmp(para.type,'unpaired')
      if kperm == 1
        dat = squeeze(cat(3,fc(:,:,:,1,:),fc(:,:,:,2,:)));
      end
      idx = randperm(size(fc,3)*2);
      tmpdat = dat(:,:,idx,:,:);
      permdat(:,:,:,1,:) = tmpdat(:,:,1:size(fc,3),:);
      permdat(:,:,:,2,:) = tmpdat(:,:,size(fc,3)+1:end,:);
      
    elseif  strcmp(para.type,'paired')
      iiperm = (iperm-1)*par.subs+kperm;
      idx1 = all_idx1(:,iiperm);
      idx2 = 3-idx1;
      for i = 1 : length(idx1)
        permdat(:,:,i,1,:) = fc(:,:,i,idx1(i),:);
        permdat(:,:,i,2,:) = fc(:,:,i,idx2(i),:);
      end
    else
      error('No test type provided. Needs to be either paied or unpaired=.')
    end
    
    % compute ttest during rest and atomoxetine
    [t_res1,~,~,s] = ttest(atanh(permdat(:,:,:,2,:)),atanh(permdat(:,:,:,1,:)),'dim',3,'alpha',alp);
    t_res1 = squeeze(t_res1.*sign(s.tstat)); clear s
    
    % -----------------------
    % FRACTION OF ALTERED CORRELATONS - global (across space)
    % -----------------------
    % count fraction of altered connections
    perm.tperm_n(kperm,:)=100*nansum(nansum(t_res1<0))./nvox;
    perm.tperm_p(kperm,:)=100*nansum(nansum(t_res1>0))./nvox;
    % number of fraction connections, irrespective of direction
    perm.tperm_all(kperm,:)=100*nansum(nansum(abs(t_res1)))./nvox;
    
    % -----------------------
    % FRACTION OF ALTERED CORRELATONS - local (per voxel)
    % -----------------------
    % count fraction of altered connections
    perm.tperm_pervoxel_neg(:,kperm,:)=100*nansum(t_res1<0)./size(fc,1);
    perm.tperm_pervoxel_pos(:,kperm,:)=100*nansum(t_res1>0)./size(fc,1);
    % fraction of altered connections, irrespective of direction
    perm.tperm_pervoxel_all(:,kperm,:)=100*nansum(abs(t_res1))./size(fc,1);
    
  end
  
  save(sprintf('tmp_permtest_iperm%d_nperm.mat',iperm,para.nperm),'perm')
  
  %   try
  %     pause(randi(3))
  %     load(sprintf([outdir 'pupmod_src_powcorr_permtest_iperm%d_nperm%d_v%d.mat'],iperm,nperm,v))
  %   catch me
  %     save(sprintf([outdir 'pupmod_src_powcorr_permtest_iperm%d_nperm%d_v%d.mat'],iperm,nperm,v),'par')
  %   end
  
end

if strcmp(para, 'collect_result')
  if ~isfield(para,'emp')
    error('Empirical results need to be provided in order to compute p-values!')
  end
  
  load(sprintf(['tmp_permtest_iperm%d_nperm.mat'],iperm,para.nperm))
  delete(sprintf(['tmp_permtest_iperm%d_nperm.mat'],iperm,para.nperm))
  
  emp = para.emp;
  
  if strcmp(para.correction_method,'single_threshold')
  
  % Single threshold permutation test: global first
  % ------------------
  idx_R_pos   = max(abs(perm.tperm_p),[],2);
  idx_R_neg   = max(abs(perm.tperm_n),[],2);
  idx_R_all   = max(abs(perm.tperm_all),[],2);
  
  for ifreq = 1:size(fc,5)
    outp.p_pos(ifreq) = 1-sum(abs(emp.n_pos(ifreq))>abs(idx_R_pos))/para.nperm;
    outp.p_neg(ifreq) = 1-sum(abs(emp.n_neg(ifreq))>abs(idx_R_neg))/para.nperm;
    outp.p_all(ifreq) = 1-sum(abs(emp.n_all(ifreq))>abs(idx_R_all))/para.nperm;
    
  end
  
  % Single threshold permutation test: local second
  % ------------------
  clear d_max_p d_max_n
  
  for ifreq = 1:size(fc,5)
    
    d_max_pos  = squeeze(max(perm.tperm_pervoxel_pos(:,:,ifreq),[],1));
    d_max_neg  = squeeze(max(perm.tperm_pervoxel_neg(:,:,ifreq),[],1));
    d_max_all  = squeeze(max(perm.tperm_pervoxel_neg(:,:,ifreq),[],1));
    for ivox = 1 : size(fc,1)
      outp.pval_p_atx(ivox,ifreq) = 1-sum(squeeze(emp.n_pos_pervoxel(ivox,ifreq,:))' > d_max_pos)./para.nperm;
      outp.pval_n_atx(ivox,ifreq) = 1-sum(squeeze(emp.n_neg_pervoxel(ivox,ifreq,:))' > d_max_neg)./para.nperm;
      outp.pval_n_atx(ivox,ifreq) = 1-sum(squeeze(emp.n_neg_pervoxel(ivox,ifreq,:))' > d_max_neg)./para.nperm;
    end
  end
  
  elseif strcmp(para.correction_method,'ranks')
  % Single threshold permutation test based on ranks
  % ------------------
  % global first
  
    for ifreq = 1:size(fc,5)
      [~,tmp] = sort(perm.tperm_n(:,ifreq),'ascend');
      R_n(tmp,ifreq) = 1:numel(perm.tperm_n(:,ifreq));     
      [~,tmp] = sort(perm.tperm_p(:,ifreq),'ascend');
      R_p(tmp,ifreq) = 1:numel(perm.tperm_p(:,ifreq));   
    end
    
    Rmax_n = max(cat(2,R_n,R_p),[],2);
    Rmax_p = max(cat(2,R_n,R_p),[],2);

    for ifreq = 1:size(fc,5)
      fprintf('Obtaining corrected p-values: freq%d ...\n',ifreq)
      for irank = 1 : para.nperm    
        Dmax_p(irank,ifreq)     = perm.tperm_p( R_p(:,ifreq) == Rmax_p(irank), ifreq, 1);
        Dmax_n(irank,ifreq)     = perm.tperm_n( R_n(:,ifreq) == Rmax_n(irank), ifreq, 1);
      end
      % ---------------------------
      % OBTAIN CORRECTED P-VALUES
      % ---------------------------
      outp.p_pos(ifreq)     = 1-sum(abs(emp.n_pos(ifreq))>abs(Dmax_p(:,ifreq)))/para.nperm;
      outp.p_neg(ifreq)     = 1-sum(abs(emp.n_neg(ifreq))>abs(Dmax_n(:,ifreq)))/para.nperm;

    end
  
  % local second
    
    
    
  end
else
  outp = perm;
end


