function tp_clustperm_corr(z,para)

% COMPUTES CLUSTER-BASED PERMUTATION TEST IN SOURCE SPACE FOR CORRELATION
% VALUES (e.g. correlation between source activity and external variable)
% Takes the following inputs:
% - z:    correlation values
% - para: see tp_clusterperm



t_pos = s.stat.*(s.stat>s.critval);
  t_neg = s.stat.*(s.stat<-s.critval);
  
  s.stat_pos = compute_cluster(t_pos,alldat,n,para);
  s.stat_neg = compute_cluster(t_neg,alldat,n,para);
  
  if ~para.perm
    
    if ~isempty(s.stat_pos)
      fprintf('Positive clusters found: %d \n',size(s.stat_pos,2))
    end
    if ~isempty(s.stat_pos)
      fprintf('Negative clusters found: %d \n',size(s.stat_neg,2))
    end
  end
  
end


% IDENTIY CLUSTER HERE
function cluster = compute_cluster(t,alldat,n,para)


done = [];

idx = 1:size(alldat,1);
cnt = 0;

done = [done find(t==0)'];
idx(t==0) = NaN;

while sum(~isnan(idx))~=0
  
  if sum(isnan(idx))==0
    i = idx(1);
  else
    i = find(~isnan(idx),1,'first');
  end
  
  done = [done i];
  
  neigh = find(n(i,:));
  
  tmp_clust = [];
  
  sign_neigh = neigh(find(abs(t(neigh))>0));
  
  if ~isempty(sign_neigh) && length(sign_neigh) >= para.minneigh
    
    done = [done sign_neigh];
    
    tmp_clust = [i];
    
    while (~isempty(sign_neigh) || length(done)>= length(t))  && length(sign_neigh) >= para.minneigh
      
      remove = [];
      
      for ineigh = 1 : length(sign_neigh)
        if sum(abs(t(find(n(sign_neigh(ineigh),:))))>0)<para.minneigh
          remove = [remove ineigh];
        end
      end
      
      sign_neigh(remove)=[];
      
      if isempty(sign_neigh)
        break
      end
      
      if length(unique([sign_neigh(:); tmp_clust(:)]))==length(tmp_clust)
        break
      end
      
      tmp_clust = unique([tmp_clust sign_neigh(:)']);
      
      [~,nneigh]=find(n(sign_neigh,:));
      
      nneigh = unique(nneigh);
      
      done   = [done nneigh(:)'];
      
      done   = unique(done);
      
      sign_neigh = nneigh(find(abs(t(nneigh))>0));
      
    end
    
    
    if length(tmp_clust) >= para.minneigh
      cnt = cnt + 1;
      cluster{cnt}.chan	= tmp_clust;
      cluster{cnt}.stat 	= sum(t(tmp_clust));
    end
    
    clear tmp_clust sign_neigh nneigh
    
  end
  
  idx(done) = NaN; done = [];
  
end

if cnt == 0
  cluster = [];
end


function pos_cluster  = compute_clusterfreq(pos,nsamples)

  nfreq = size(pos,2);

  for ifreq = 1 : nfreq
    
    chan{ifreq} = zeros(nsamples, size(pos{ifreq},2));
    
    for iclust = 1 : size(pos{ifreq},2)
      chan{ifreq}(pos{ifreq}{iclust}.chan,iclust) = 1;
      stat{ifreq}(iclust)   = pos{ifreq}{iclust}.stat;
      pos{ifreq}{iclust}.freq = ifreq;
    end
  end
  
  mintempneigh = 10;
  clust   = 0;
  freq    = 0;
  allfreq = 1:nfreq;
  
  
  for freq  = 1:nfreq
    for clust = 1 : size(chan{freq},2) 

      afreq = find(allfreq~=freq);

      t = chan{freq}(:,clust);
      if sum(t) == 0
        continue
      end

      for ifreq = afreq

        for iclust = 1 : size(chan{ifreq},2)
          
          if isempty(pos{ifreq}{iclust})
            continue
          end
          
          isneigh(iclust) = sum((t+chan{ifreq}(:,iclust))==2)>mintempneigh;

          if isneigh(iclust)

            pos{freq}{clust}.chan = unique([pos{freq}{clust}.chan pos{ifreq}{iclust}.chan]);
            pos{freq}{clust}.stat = pos{freq}{clust}.stat+pos{ifreq}{iclust}.stat;
            pos{freq}{clust}.freq = unique([pos{freq}{clust}.freq pos{ifreq}{iclust}.freq]);
            chan{freq}(pos{freq}{clust}.chan,clust) = 1;
            chan{ifreq}(:,iclust) = 0;
            t = chan{freq}(:,clust);
            pos{ifreq}{iclust} = [];

          end
        end
        clear isneigh
      end
    end
  end 
    
  % COLLECT RESULTS
  cnt = 0;
  for ifreq = 1 : nfreq
    for iclust = 1 : size(pos{ifreq},2)
      if ~isempty(pos{ifreq}{iclust})
        cnt = cnt + 1;
        pos_cluster{cnt} = pos{ifreq}{iclust};
      end
    end
  end
  
 if ~exist('pos_cluster','var')
   pos_cluster = [];
 end

% SOME ALTERNATIVE FIELDTRIP CODE
% -------------------------------------
% posclusobs = pconn_findcluster(t_pos>0,n,para.minneigh);
% s.numclust = max(posclusobs);

% for iclust = 1 : s.numclust
%   s.clusterstat(iclust) = sum(t_pos(posclusobs==iclust));
%   s.chan{iclust}=posclusobs==iclust;
% end
% -------------------------------------

  