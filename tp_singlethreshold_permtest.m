function [h, p] = tp_singlethreshold_permtest(dat,para)

% Computes single threshold permutation test for paired samples
% (currently only two-tailed test supported)
% -----------------------------------
% dat: original data, N_observations x N_subj x N_cond
% para.nperm: number of permutations (default: 10000)
% para.alpha: alpha level (default: 0.05)
% para.thresholdalpha: first stats alpha level (default: 0.05)
% para.tail: -1/0/1 (default: 0, two-tailed)
% ---------------------
% thmspfffr@gmail.com
% --------------------

if ~isfield(para,'nperm')
  para.nperm = 10000;
end
if ~isfield(para,'alpha')
  para.alpha = 0.05;
end
if ~isfield(para,'thresholdalpha')
  para.thresholdalpha = 0.05;
if ~isfield(para,'tail')
  para.tail = 0;
end

[~,~,~,s]=ttest(dat(:,:,2),dat(:,:,1),'dim',2,'alpha',para.thresholdalpha,'tail','both');
t_emp = s.tstat;

for iperm = 1 : para.nperm
  
  fprintf('Permutation #%d\n', iperm);
      
  idx1 = randi(2,[size(dat,2),1]);
  idx2 = 3-idx1;
  
  for i = 1 : size(dat,2)
    
    x(:,i,1) = dat(:,i,idx1(i));
    x(:,i,2) = dat(:,i,idx2(i));
    
  end
  
  [~,~,~,s]=ttest(x(:,:,2),x(:,:,1),'dim',2,'alpha',para.thresholdalpha,'tail','both');
  t_max(iperm) = max(s.tstat);

end

for i = 1 : size(t_emp,1)
  p(i) = 1-sum(abs(t_emp(i)) > abs(t_max)) / para.nperm;
  h(i) = p(i)<para.alpha;
end

