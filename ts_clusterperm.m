function stats = clust_perm(z,para)


para.perm = 0;
s = find_cluster_sub(z,para);

% PERMUTATION
% -----------------------------------------------------------------------
% FIRST COMPUTE STATS ON ALL CLUSTERS

para.perm = 1;
if ~(isempty(s.stat_neg) && isempty(s.stat_pos))
  
  % -----------------------------------------------------------------------
  % LOOP THROUGH ALL PERMUTATIONS
  % -----------------------------------------------------------------------
  % note that data is shuffled within find_cluster_sub (maybe change that)
  
  for iperm = 1 : para.nperm;
    
    fprintf('Permutation %d ...\n',iperm);
    
    perm = find_cluster_sub(z,para);
    
    if ~isempty(perm.stat_pos)
      
      for iclust = 1 : size(perm.stat_pos,2)
        tmp(iclust) = perm.stat_pos{iclust}.stat;
      end
      
      pos_perm(iperm) = max(tmp); clear tmp
      
    else
      pos_perm(iperm) = 0;
    end
    
    if ~isempty(perm.stat_neg)
      
      for iclust = 1 : size(perm.stat_neg,2)
        tmp(iclust) = perm.stat_neg{iclust}.stat;
      end
      neg_perm(iperm) = max(tmp); clear tmp
    else
      neg_perm(iperm) = 0;
    end
    
  end
  % end of permutation
  
  % -----------------------------------------------------------------------
  % DETERMINE CLUSTER P-VALUES
  % -----------------------------------------------------------------------
  
  if ~isempty(s.stat_pos)
    for iclust = 1 : size(s.stat_pos,2)
      s.stat.p_pos(iclust) = sum(pos_perm>s.stat_pos{iclust}.stat)/para.nperm;
    end
  end
  if ~isempty(s.stat_neg)
    for iclust = 1 : size(s.stat_neg,2)
      s.stat.p_neg(iclust) = sum(neg_perm<s.stat_neg{iclust}.stat)/para.nperm;
    end
  end
  
  % -----------------------------------------------------------------------
  % CREATE STATISTICAL MASK
  % -----------------------------------------------------------------------
  
  s.mask = zeros(size(z,1),1);
  
  if ~isempty(s.stat_pos)
    tmp = find(s.stat.p_pos<0.05);
    if ~isempty(tmp)
      for itmp = 1 : length(tmp)
        s.mask(s.stat_pos{itmp}.chan) = 1;
      end
    end
  end
  if ~isempty(s.stat_neg)
    tmp = find(s.stat.p_neg<0.025);
    if ~isempty(tmp)
      for itmp = 1 : length(tmp)
        s.mask(s.stat_neg{itmp}.chan) = 1;
      end
    end
  end
  
  % -----------------------------------------------------------------------
  % HOLM CORRECTION
  % -----------------------------------------------------------------------
  % this is in experimental stage. removes previously identified significant 
  % cluster from data and repeats search for clusters. continues doing so
  % until no more significant clusters are found. 
  % also described in hipp et al. (2011).
  
  % think about putting into separate function!
  if ~isfield(para,'correctm')
    para.correctm = [];
  end
  
  if strcmp(para.correctm,'holm')
    
    clear neg_perm pos_perm
    
    br = 0;
    
    if ~isempty(s.stat_pos)
      num_sign_clust(1) = sum(s.stat.p_pos<para.alpha);
    else
      num_sign_clust(1) = 0;
    end
    
    if ~isempty(s.stat_neg)
      num_sign_clust(2) = sum(s.stat.p_neg<para.alpha);
    else
      num_sign_clust(2) = 0;
    end
    
    while any(num_sign_clust>0)
      
      if num_sign_clust(2)>0
      
        idx = s.stat_neg{find(s.stat.p_neg<para.alpha,1,'first')}.chan;

        z(idx,:,:)=NaN;

        para.neigh(idx,:)=NaN;
        para.neigh(:,idx)=NaN;

        para.perm = 0;

        tmp_s = find_cluster_sub(z,para);

        fprintf('Removed one negative cluster - searching again!\n');

        pause(0.5);

        if isempty(tmp_s.stat_neg) && isempty(tmp_s.stat_pos)

          fprintf('No more clusters found!\n')

          br = 1;

          break

        else
        
          for iperm = 1 : para.nperm;
            
            para.perm = 1;

            fprintf('Permutation #%d ...\n',iperm);

            perm = find_cluster_sub(z,para);

            if ~isempty(perm.stat_neg)

              for iclust = 1 : size(perm.stat_neg,2)
                tmp(iclust) = perm.stat_neg{iclust}.stat;
              end

              neg_perm(iperm) = max(tmp); clear tmp

            else
              neg_perm(iperm) = 0;
            end
          end

          if ~isempty(tmp_s.stat_neg)
            for iclust = 1 : size(tmp_s.stat_neg,2)
              tmp_s.stat.p_neg(iclust) = sum(neg_perm<tmp_s.stat_neg{iclust}.stat)/para.nperm;
            end
          end
          
          if ~any(tmp_s.stat.p_neg<para.alpha)
            
            fprintf('No more significant clusters found!\n') 
            
            num_sign_clust(2) = num_sign_clust(2) - 1;
            
            break
            
          else
            
            fprintf('Significant cluster(s) found!\n') 
            
            % NOW GO BACK UP AND REMOVE SHIT AGAIN
            
          end
            
            
          

        end      
      end   
    end
  end
  % end of holm correction
  % -----------------------------------------------------------------------
  
else
  fprintf('No clusters found.\n');
end

stats = s; clear s
