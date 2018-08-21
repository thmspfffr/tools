function deg = tp_degree(fc_tmp,para)

% para.alpha = 0.01;
% para.nfreq = 13;

fprintf('Computing degree centrality across %d frequencies ...\n',para.nfreq)
fprintf('Alpha-level at %.2f ...\n',para.alpha)

fc_tmp = squeeze(fc_tmp);
siz = size(fc_tmp,1);

x_ref   = zeros(siz*siz,28,2,para.nfreq,'single');
x_mean  = zeros(siz*siz,28,2,para.nfreq,'single');
x_std   = zeros(siz*siz,28,2,para.nfreq,'single');

for il = 1 : siz
    
  x_ref((il-1)*siz+1:(il)*siz,:,:,1:para.nfreq)  = squeeze(abs(fc_tmp(il,:,:,1:2,1:para.nfreq)));
  x_mean((il-1)*siz+1:(il)*siz,:,:,1:para.nfreq) = repmat(nanmean(abs(squeeze(fc_tmp(il,:,:,1:2,1:para.nfreq)))), [siz 1 1 1]);
  x_std((il-1)*siz+1:(il)*siz,:,:,1:para.nfreq)  = repmat(nanstd(abs(squeeze(fc_tmp(il,:,:,1:2,1:para.nfreq)))), [siz 1 1 1]);
  
end

z_ctrl = ( x_ref(:,:,1,1:para.nfreq) - x_mean(:,:,1,1:para.nfreq) ) ./  x_std(:,:,1,1:para.nfreq);
z_expr = ( x_ref(:,:,2,1:para.nfreq) - x_mean(:,:,1,1:para.nfreq) ) ./  x_std(:,:,1,1:para.nfreq);

clear x_ref x_mean x_std

[~,p_ctrl] = ttest(zeros(size(z_ctrl)),z_ctrl,'tail','left','dim',2);
[~,p_expr] = ttest(zeros(size(z_expr)),z_expr,'tail','left','dim',2);

th_ctrl = squeeze(p_ctrl < para.alpha/2);
th_expr = squeeze(p_expr < para.alpha/2);

deg = nan(siz,siz,para.nfreq,2);

fprintf('Computing adjacency matrix ...\n')

for i = 1 : siz
  for j = 1 : siz
    
    if i==j
      continue
    end
    
    deg(i,j,1:para.nfreq,1) = th_ctrl( (i-1)*siz+j,1:para.nfreq ) + th_ctrl( (j-1)*siz+i,1:para.nfreq ) > 0;
    deg(i,j,1:para.nfreq,2) = th_expr( (i-1)*siz+j,1:para.nfreq ) + th_expr( (j-1)*siz+i,1:para.nfreq ) > 0;
  
  end
end


 