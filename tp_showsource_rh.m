function tp_showsource_rh(par_interp,cmap,sa_meg_template,para)

viewdir   = [0.001 -0 1; 0 -.5 0; 1 0 0; -.5 0 0; -.5 0 0;1 0 0 ];

figure; set(gcf,'color','white'); hold on;
set(gcf,'Position',[600 50 1300 1050])

for iplot = 1 : 6
  
  h=subplot(3,2,iplot);
  
  if iplot == 1 || iplot == 3 || iplot == 5
    get(h,'pos');
    set(h,'pos',[ans(1)+0.10 ans(2) ans(3) ans(4) ]);
  else
    get(h,'pos');
    set(h,'pos',[ans(1)-0.10 ans(2) ans(3) ans(4) ]);
  end
  
  if iplot == 3 || iplot == 4
    get(h,'pos');
    set(h,'pos',[ans(1) ans(2)+0.08 ans(3) ans(4) ]);
  elseif iplot == 5 || iplot == 6
    get(h,'pos');
    set(h,'pos',[ans(1) ans(2)+0.16 ans(3) ans(4) ]);
  end
  
  para.myviewdir = viewdir(iplot,:);
%   a = sa_meg_template.cortex10K;
  
%   if iplot == 5
    a=cutsurface(sa_meg_template.cortex10K,[mean(sa_meg_template.cortex10K.vc(:,1)) 0 0],[1 0 0]);
%   elseif iplot == 6
%     a=cutsurface(sa_meg_template.cortex10K,[mean(sa_meg_template.cortex10K.vc(:,1)) 0 0],[-1 0 0]);
%   end
  
  pconn_showsurface(a,para,par_interp);
  colormap(cmap)
  camlight headlight
  
end