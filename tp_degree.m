%% Compute degree across a variable range of frequencies
% ------------------------
% Degree is defined relative to a placebo condition ((or, more generally: a
% control condition). Note that this deviates from other definitions of
% degree (e.g. van den Brink et al., 2016).
% ------------------------

% para.alpha = 0.01;
% para.nfreq = 13;

function deg = tp_degree(fc_tmp,para)


  fprintf('Computing degree centrality across %d frequencies ...\n',para.nfreq)
  fprintf('Alpha-level at %.2f ...\n',para.alpha)

  % FC is of size nvox*nvox*nsubj*ncond*nfreq, 
  % where ncond = 1 (placebo) and 2 (drug)
  fc_tmp = squeeze(fc_tmp);
  siz = size(fc_tmp,1);

  x_ref   = zeros(siz*siz,28,2,para.nfreq,'single');
  x_mean  = zeros(siz*siz,28,2,para.nfreq,'single');
  x_std   = zeros(siz*siz,28,2,para.nfreq,'single');

  if para.absolute == 0
    for il = 1 : siz
      % x_ref:  connections from node il to all others
      x_ref((il-1)*siz+1:(il)*siz,:,:,1:para.nfreq)  = squeeze(fc_tmp(il,:,:,1:2,1:para.nfreq));
      % x_mean: mean connection strength from node il
      x_mean((il-1)*siz+1:(il)*siz,:,:,1:para.nfreq) = repmat(nanmean(squeeze(fc_tmp(il,:,:,1:2,1:para.nfreq))), [siz 1 1 1]);
      % x_std:  std of correlations from node il
      x_std((il-1)*siz+1:(il)*siz,:,:,1:para.nfreq)  = repmat(nanstd(squeeze(fc_tmp(il,:,:,1:2,1:para.nfreq))), [siz 1 1 1]);
    end
  elseif para.absolute == 1
    for il = 1 : siz
        % x_ref:  connections from node il to all others
        x_ref((il-1)*siz+1:(il)*siz,:,:,1:para.nfreq)  = squeeze(abs(fc_tmp(il,:,:,1:2,1:para.nfreq)));
        % x_mean: mean connection strength from node il
        x_mean((il-1)*siz+1:(il)*siz,:,:,1:para.nfreq) = repmat(nanmean(squeeze(abs(fc_tmp(il,:,:,1:2,1:para.nfreq)))), [siz 1 1 1]);
        % x_std:  std of correlations from node il
        x_std((il-1)*siz+1:(il)*siz,:,:,1:para.nfreq)  = repmat(nanstd(squeeze(abs(fc_tmp(il,:,:,1:2,1:para.nfreq)))), [siz 1 1 1]);
    end
  end
  
  % Compute "relative" degree: Reference connection from placebo / control
  if para.relative_degree == 1
    % z-score for placebo (or rest)
    z_ctrl = ( x_ref(:,:,1,1:para.nfreq) - x_mean(:,:,1,1:para.nfreq) ) ./  x_std(:,:,1,1:para.nfreq);
    % z-score for drug (or task)
    z_expr = ( x_ref(:,:,2,1:para.nfreq) - x_mean(:,:,1,1:para.nfreq) ) ./  x_std(:,:,1,1:para.nfreq);
  
  % Compute "absolute" degree: Reference connection from drug / exp
  elseif para.relative_degree == 0 
    % z-score for placebo (or rest)
    z_ctrl = ( x_ref(:,:,1,1:para.nfreq) - x_mean(:,:,1,1:para.nfreq) ) ./  x_std(:,:,1,1:para.nfreq);
    % z-score for drug (or task)
    z_expr = ( x_ref(:,:,2,1:para.nfreq) - x_mean(:,:,2,1:para.nfreq) ) ./  x_std(:,:,2,1:para.nfreq);
  
  end
  
  clear x_ref x_mean x_std

  [~,p_ctrl] = ttest(zeros(size(z_ctrl)),z_ctrl,'tail','left','dim',2);
  [~,p_expr] = ttest(zeros(size(z_expr)),z_expr,'tail','left','dim',2);

  th_ctrl = squeeze(p_ctrl < para.alpha/2);
  th_expr = squeeze(p_expr < para.alpha/2);

  deg = nan(siz,siz,para.nfreq,2);

  fprintf('Computing adjacency matrix ...\n')

  % obtain adjacency matrix: sum across both directions, i.e., connectin is
  % present if both A <-> B is larger than A <-> all others or A <-> B is
  % larger than B <-> all others. See Hipp et al. for details.
  for i = 1 : siz
    for j = 1 : siz

      if i==j
        continue
      end

      deg(i,j,1:para.nfreq,1) = th_ctrl( (i-1)*siz+j,1:para.nfreq ) + th_ctrl( (j-1)*siz+i,1:para.nfreq ) > 0;
      deg(i,j,1:para.nfreq,2) = th_expr( (i-1)*siz+j,1:para.nfreq ) + th_expr( (j-1)*siz+i,1:para.nfreq ) > 0;

    end
  end
  
end


 