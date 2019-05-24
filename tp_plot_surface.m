function f = tp_plot_surface(par,para)

load /home/gnolte/meth/templates/sa_template.mat;
addpath /home/gnolte/meg_toolbox/meg/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/toolbox/

vc    = sa_template.vc;
g2    = sa_template.cortex10K.vc;

if ~isfield(para,'grid')
  error('para.grid must be defined!')
else
  g1 = para.grid;
end
if ~isfield(para,'dd')
  dd = 0.5;
else
  dd = para.dd;
end

par = spatfiltergauss(par(:),g1,dd,g2);
%%
if isfield(para,'fn')
  fig = figure('visible','off'); set(gcf,'color','w');
else
  fig = figure; set(gcf,'color','w');
end
% frontal: view = [180 0], light = [0 1 0]

positions = [0 0 1; 0 -1 0; 1 0 0; -1 0 0; 1 0 0; -1 0 0];
views = [0 90; 0 0 ; 90 0; -90 0; 90 0; -90 0];

spat_range = [min(sa_template.cortex10K.vc); max(sa_template.cortex10K.vc)];

% 01 02 xx xx 05 06
% yy yy xx xx zz zz
% yy yy kk kk zz zz
% 19 20 kk kk 23 24

% grid: 1st colum  = left/right,
%       2nd column = front/back,
%       3rd column = top/down
% fig.Renderer = 'Painters';
for i = [1 2 3 4]%size(positions,1)
  %   aa = subplot(2,2,i)
  if i == 1
    aa=subplot(4,6,[3 4 9 10]);
  elseif i == 2
    aa=subplot(4,6,[15 16 21 22]);
  elseif i == 3
    aa=subplot(4,6,[11 12 17 18]);
  elseif i == 4
    aa=subplot(4,6,[7 8 13 14]);
  end
  
  trisurf(sa_template.cortex10K.tri,sa_template.cortex10K.vc(:,1),sa_template.cortex10K.vc(:,2),sa_template.cortex10K.vc(:,3),par)
  shading interp; axis off; colormap(para.cmap);
  set(gca,'clim',[para.clim(1) para.clim(2)]);
  view(views(i,1),views(i,2)); daspect([1 1 1])
  h=light; material([0.4 0.6 0.2]); lighting GOURAUD
  set(h,'Position',positions(i,:));
  if isfield(para,'marker')
    hold on
    if mean(para.marker)>10; para.marker = para.marker./10; end
    plot3(para.marker(1),para.marker(2),10,'o')
  end
  %
  if i == 1
    set(aa,'XLim',[spat_range(1,1) spat_range(2,1)])
    set(aa,'Outerposition',[0.40 0.45 0.2447 0.4491])
    set(aa,'YLim',[spat_range(1,2)-1.8 spat_range(2,2)+1.8])
  elseif i == 2
    set(aa,'XLim',[spat_range(1,2) spat_range(2,2)])
    set(aa,'Outerposition',[0.365 0.13 0.2648 0.4341])
  elseif i == 3
    axis tight off
    set(aa,'YLim',[spat_range(1,2) spat_range(2,2)])
    set(aa,'ZLim',[spat_range(1,3) spat_range(2,3)])
    set(aa,'XLim',[spat_range(1,1) spat_range(2,1)])
    set(aa,'OuterPosition',[0.58 0.25 0.3 0.5])
  elseif i == 4
    axis tight off
    set(aa,'YLim',[spat_range(1,2) spat_range(2,2)])
    set(aa,'ZLim',[spat_range(1,3) spat_range(2,3)])
    set(aa,'XLim',[spat_range(1,1) spat_range(2,1)])
    set(aa,'OuterPosition',[0.17 0.25 0.3 0.5])
    %   elseif i == 5
    %     set(aa,'XLim',[spat_range(1,1) 0])
    %     elseif i == 6
    %     set(aa,'XLim',[spat_range(1,2) 0])
  end
  
end

% if filename is provided, crop figure
if isfield(para,'fn')
  print(gcf,'-dpng',para.fn)
  im = imread(para.fn);
  f = cell2mat(arrayfun(@(x) find(im(x,:,1)<255,1,'first'),1:size(im,1),'Uniformoutput',0));
  idx1 = min(f)-1;
  f = cell2mat(arrayfun(@(x) find(im(x,:,1)<255,1,'last'),1:size(im,1),'Uniformoutput',0));
  k = arrayfun(@(x) find(im(x,:,1)<255,1,'last'),1:size(im,1),'Uniformoutput',0);
  ff = arrayfun(@(x) ~isempty(k{x}),1:length(k));
  idx2 = max(f)+1;
  idx3 = find(ff==1,1,'first')-1;
  idx4 = find(ff==1,1,'last')+1;
  
  im = im(idx3:idx4,idx1:idx2,:);
  imwrite(im,para.fn,'png');
  
  fig = figure; set(gcf,'color','w');

  imshow(im);

end

