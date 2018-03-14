function y = scattered_barplot(x,y)
%     size x: N x 2, y: N x 2
% or  size x = N x 1, , y: N x 1

figure; set(gcf,'color','white'); hold on

if size(x,2) == 2 && size(y,2) == 2


  m_maxi = max([m_cnt_par1(istr) m_cnt_par2(istr) m_res_par1(istr) m_res_par2(istr) max(par_res(:)) max(par_cnt(:))]);
  m_mini = min([m_cnt_par1(istr) m_cnt_par2(istr) m_res_par1(istr) m_res_par2(istr) min(par_res(:)) min(par_cnt(:))]);

  r = rand(28,1)/20; r = r-mean(r);

  plot([1*ones(length(par_res),1)+r 1.2*ones(length(par_res),1)+r]',[par_res(:,1) par_res(:,2)]','color',[0.7 0.7 0.7]);
  scatter(1*ones(length(par_res),1)+r,par_res(:,1),100,'markerfacecolor','w','markeredgecolor',[0.7 0.7 0.7])
  scatter(1.2*ones(length(par_res),1)+r,par_res(:,2),100,'markerfacecolor','w','markeredgecolor',[0.7 0.7 0.7])

  plot([1 1.2],[m_res_par1(istr) m_res_par2(istr)],'o','markerfacecolor','r','markeredgecolor','w','markersize',14)
  plot([1 1.2]',[m_res_par1(istr) m_res_par2(istr)]','color','r');

  plot([1.4*ones(length(par_cnt),1)+r 1.6*ones(length(par_cnt),1)+r]',[par_cnt(:,1) par_cnt(:,2)]','color',[0.7 0.7 0.7]);
  scatter(1.4*ones(length(par_cnt),1)+r,par_cnt(:,1),100,'markerfacecolor','w','markeredgecolor',[0.7 0.7 0.7])
  scatter(1.6*ones(length(par_cnt),1)+r,par_cnt(:,2),100,'markerfacecolor','w','markeredgecolor',[0.7 0.7 0.7])

  plot([1.4 1.6],[m_cnt_par1(istr) m_cnt_par2(istr)],'o','markerfacecolor',[0 0.5 1],'markeredgecolor','w','markersize',14)
  plot([1.4 1.6]',[m_cnt_par1(istr) m_cnt_par2(istr)]','color',[0 0.5 1]);

end