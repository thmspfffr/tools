function fac = compute_altered_corr(fc,para)
%COMPUTE_ALTERED_CORR Compute fraction of altered correlations
%   fac = compute_altered_corr(data,para) computes fraction of
%   significantly altered correlations. Altered correlations are determined
%   through a (uncorrected) t-test of the fisher-transformed correlation
%   matrices.
%   ----------------
%   INPUTS
%   ----------------
%   fc: Fuctional connectivity, with size NVOX,NVOX,NSUBJ,NCOND,NFREQ
%   para.nfreq: Frequency indicies altered correlations should be computed
%   over (e.g. only frequencies with index [2 3 4])
%   para.alpha: alpha-level for first-level ttest (default: 0.05)


if ~isfield(para,'alpha')
  para.alpha = 0.05;
end
if ~isfield(para,'nfreq')
  para.nfreq = 1:length(fc,5);
end
%   

for ifoi = para.nfreq
  
  fprintf('Computing fraction of altered correlations, freq%d...\n',ifoi)
  % -------------------------
  % EXP. VS. CONTROL (e.g. Patients vs. healthy)
  % -------------------------
  [h,~,~,s]=ttest(atanh(fc(:,:,:,2,ifoi)),atanh(fc(:,:,:,1,ifoi)),'dim',3,'alpha',para.alpha);
  fac.n_pos(ifoi) = 100*nansum(nansum((h.*sign(s.tstat))>0))./(size(fc,1)*size(fc,1)-size(fc,1));
  fac.n_neg(ifoi) = 100*nansum(nansum((h.*sign(s.tstat))<0))./(size(fc,1)*size(fc,1)-size(fc,1));
 
  % Fraction of altered correlations (positive AND negative)
  % --------------------------
  % during rest (condition label = 1)
  h=ttest(fc(:,:,:,2,ifoi),fc(:,:,:,1,ifoi),'dim',3,'alpha',para.alpha);
  fac.n_all(ifoi) = 100*nansum(nansum((h)))./(size(fc,1)*size(fc,1)-size(fc,1));
  
  % Fraction of altered correlations (for each vertex)
  % --------------------------
  % during rest (condition label = 1)
  [h,~,~,s]=ttest(atanh(fc(:,:,:,2,ifoi)),atanh(fc(:,:,:,1,ifoi)),'dim',3,'alpha',para.alpha);
  fac.n_pos_pervoxel(:,ifoi) = 100*nansum((h.*sign(s.tstat))>0)./(size(fc,1)-1);
  fac.n_neg_pervoxel(:,ifoi) = 100*nansum((h.*sign(s.tstat))<0)./(size(fc,1)-1);
 
end

