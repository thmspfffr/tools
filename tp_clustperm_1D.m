function [stats] = tp_clustperm_1D(x,y,para)
% COMPUTES CLUSTER-BASED PERMUTATION TEST IN 2D
% EXAMPLE: differences between two time series
% if ~exist(para)
%   para.alpha      = 0.05;
%   para.clustalpha = 0.05;
%   para.tail       = 0;
%   para.paried     = 0;
%   para.nperm      = 10000;
% end
%
rng('shuffle','twister')

if ~isfield(para,'clusteralpha')
  fprintf('No cluster alpha set. Assuming default a = 0.05 ...\n')
  para.clusteralpha = 0.05;
end
if ~isfield(para,'alpha') 
  fprintf('No alpha set. Assuming default a = 0.05 ...\n')
  para.alpha = 0.025;
end
if ~isfield(para,'tail') 
  fprintf('Tail not set (Default: 0) ...\n')
  para.tail = 0;
end 
if ~isfield(para,'paired') 
  fprintf('Test not defined (Default: paired)...\n')
  para.paried = 1;
end 
if ~isfield(para,'nperm') 
  fprintf('Number of permutations undefined (Default: N = 10000)! \n')
  para.nperm = 10000;
end 

df = size(x,2)-1;
para.thresh = abs(tinv(para.clusteralpha,df));

% d_emp = mean(x,2)-mean(y,2);
[~,~,~,s]=ttest(y,x,'dim',2);
% --------------------------------------------
% EMPIRLCAL CLUSTERS
% --------------------------------------------
sign = find( abs(s.tstat) > para.thresh );
sign(find(diff(sign)>1))=[];

% possign = find( (s.tstat) >  para.thresh );

clear tmp_clust

if isempty(sign)
  empclust = min(s.tstat);
elseif (sum( diff(sign))-(length(sign)-1))==0
  empclust = sum(s.tstat(sign));
  stats.clustloc{1} = sign;
else
  i = find(diff(sign)>1);

  
  for iclust = 1 : length(i)+1
    if iclust == 1 
      empclust(iclust) = sum(s.tstat(sign(1:i(iclust))));
      stats.clustloc{iclust} = sign(1:i(iclust));
    elseif iclust == length(i)+1
      empclust(iclust) = sum(s.tstat(sign(i(iclust-1)+1:end)));
      stats.clustloc{iclust} = sign(i(iclust-1)+1:end);
    else
      empclust(iclust) = sum(s.tstat(sign(i(iclust-1)+1:i(iclust))));
      stats.clustloc{iclust} = sign(i(iclust-1)+1:i(iclust));
    end
  end
end

clear s i 
%%
if ~para.paired
  
  dat = [x(:,1);x(:,2)];
  
  for iperm = 1 : para.nperm
    
    idx = randperm(size(dat,1));
    
    x = dat(idx(1:size(dat,1)/2));
    y = dat(idx(size(dat,1)/2+1:end));
    
    d(iperm) = mean(x)-mean(y);
    
  end
  
elseif para.paired
  
  dat(:,:,1) = x; clear x
  dat(:,:,2) = y; clear y
  
  [nsamples, nsubj, ncond]  = size(dat);
  
  for iperm = 1 : para.nperm
    
    fprintf('Perm #%d\n',iperm);
    
    idx1 = randi(ncond,[nsubj,1]);
    idx2 = 3-idx1;
    
    for i = 1 : nsubj
      
      x(:,i) = dat(:,i,idx1(i));
      y(:,i) = dat(:,i,idx2(i));
      
    end
    
    [~,~,~,s]=ttest(y,x,'dim',2);
    
    negsign = find( abs(s.tstat) > para.thresh );
    negsign(find(diff(negsign)>1))=[];

%     possign = find( (s.tstat) >  para.thresh );
    
    clear tmp_clust
    
    if isempty(negsign)
      negmaxclust(iperm) = min(s.tstat);
    elseif (sum( diff(negsign))-(length(negsign)-1))==0
      negmaxclust(iperm) = sum(s.tstat(negsign));
    else
      i = find(diff(negsign)>1);
      for iclust = 1 : length(i)+1
        if iclust == 1
          tmp_clust(iclust) = sum(s.tstat(negsign(1:i(iclust))));
        elseif iclust == length(i)+1
          tmp_clust(iclust) = sum(s.tstat(negsign(i(iclust-1)+1:end)));
        else
          tmp_clust(iclust) = sum(s.tstat(negsign(i(iclust-1)+1:i(iclust))));
        end
      end
      negmaxclust(iperm) = max(tmp_clust);
    end
   
  end
  
end
%%

if length(empclust)==1
  switch para.tail
    case 0
      p = 1-sum(abs(empclust)>abs(negmaxclust))/para.nperm;
    case 1
      p = 1-sum(empclust>negmaxclust)/para.nperm;
    case -1
      p = 1-sum(empclust<negmaxclust)/para.nperm;
  end
  if p<para.alpha
      h = 1;
    else
      h = 0;
    end
else
  for i = 1 : length(empclust)
    switch para.tail
      case 0
        p(i) = 1-sum(abs(empclust(i))>abs(negmaxclust))/para.nperm;
      case 1
        p(i) = 1-sum(empclust(i)>negmaxclust)/para.nperm;
      case -1
        p(i) = 1-sum(empclust(i)<negmaxclust)/para.nperm;
    end
    if p<para.alpha
      h(i) = 1;
    else
      h(i) = 0;
    end
  end
end



stats.h = h;
stats.p = p;
stats.clust = empclust;






