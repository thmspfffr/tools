function [stats] = tp_permtest2D(x,y,para)
% COMPUTES CLUSTER-BASED PERMUTATION TEST IN 2D
% EXAMPLE: differences between two time series
% input: N_SUBJ x N_FREQ/N_TIME/N_X
% if ~exist(para)
%   para.alpha      = 0.05;
%   para.clustalpha = 0.05;
%   para.tail       = 0;
%   para.paried     = 0;
%   para.nperm      = 2000;
% end
%
    
para= [];
if ~isfield(para,'clusteralpha'); para.clusteralpha = 0.05; end
if ~isfield(para,'alpha'); para.alpha = 0.05; end
if ~isfield(para,'tail'); para.tail = 0; end
if ~isfield(para,'paried'); para.paired = 1; end
if ~isfield(para,'nperm'); para.nperm =10000; end

df = size(x,3)-1;
% --------------------------------------------
% EMPIRLCAL CLUSTERS
% --------------------------------------------
if para.tail       == 0
  [~,p,~,s]=ttest(x,y,'dim',3,'alpha',para.clusteralpha);
%   para.thresh = abs(tinv(para.clusteralpha,df));
else
  para.thresh = abs(tinv(para.clusteralpha,df));
  [~,p,~,s]=ttest(x,y,'dim',2,'alpha',para.clusteralpha);
end

sign = p<0.05;
clust = bwlabeln(sign,4);

for iclust = 1: max(clust(:))
  tmax_emp(iclust) = sum(s.tstat(clust==iclust));
end

%%
if para.paired
  
  dat(:,:,:,1) = x; clear x
  dat(:,:,:,2) = y; clear y
  
  [~,~,nsubj,ncond]  = size(dat);
  
  for iperm = 1 : para.nperm
    
    fprintf('Perm #%d ...\n',iperm);
    
    idx1 = randi(ncond,[nsubj,1]);
    idx2 = 3-idx1;
    
    for i = 1 : nsubj
      
      x(:,:,i) = dat(:,:,i,idx1(i));
      y(:,:,i) = dat(:,:,i,idx2(i));
      
    end
    
    [~,p,~,s]=ttest(y,x,'dim',3,'alpha',para.clusteralpha);
    
    sign = p<0.05;
    clust = bwlabeln(sign,4);
    
    if max(clust(:))>0
      for iclust = 1: max(clust(:))
        tmp(iclust) = sum(s.tstat(clust==iclust));
      end
    else
      tmp = 0;
    end
    
    tmax_perm(iperm) = max(abs(tmp));

  end
  
end

for iclust = 1 : length(tmax_emp)
 stats.p_clust(iclust) = 1 -  (sum(abs(tmax_emp(iclust)) > tmax_perm)) ./ para.nperm;
 stats.h_clust(iclust) = stats.p_clust(iclust)<para.alpha;
end





