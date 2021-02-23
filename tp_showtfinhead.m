function showtfinhead(data,locs,pars);
% usage showtfinhead(data,locs,pars);
% makes head-in-head plots e.g. to show time-frequency (or frequency frequency) plots for all sensors 
% topographically arranged. In the following the third index is assumed to
% mean time, but that could also be any other quantity. 
%
% data is an NxNfxNt  matrix where N is the number of channels, Nf number of frequencies,
% and Nt number of time points, % locs is  an Nx3 or (usually) Nx5 matrix, 
% where the second and third column denote and y coordinates of sensors. 
% Columns beyond the third are ignored.  
%      1st column: a channel i gets a subplot only if locs(i,1)>.5 scal
%                  if locs is an Mx5 or Mx3 matrix with M<N
%                  the absolute values of the first column 
%                  are interpreted as indices of scalp-electrodes.
%                  If the values themselves are smaller than 0 then
%                  these electrodes are ignored  
%                 
%      2nd and 3rd column: x,y coordinates of electrodes in 2D
%      
% pars sets parameters; 
%      pars.scale sets color-scale; it is a 1x2 vector denoting min 
%                 and max of colorbar. defaults is  pars.scale=[-max(max(abs(data))),max(max(abs(data)))];
%      pars.resolution sets resolution. Default is pars.resolution=25
%                      Increasing it makes better pictures but is slower
%      pars.global_size  sets global size factor of small circles and distancec between circles.                   
%                        default: global_size=1
%      pars.relative_size  sets size factor for distances leaving circle-size unchanged                   
%                        default: relative_size=1
%      pars.head_up    moves the big circle (the head) up. Default: head_up=0
%      pars.head_right    moves the big circle (the head) to the right. Default: head_right=0
%      pars.freqaxis  Nfx1 vector: values of frequencies (defaults is 1:Nf)
%      pars.timeaxis Ntx1 vector: values of times (defaults is 1:Nt)
%      pars.tunit string to label the last index of the input data, which
%                  corrsponds to the x-axis in the plots (default is 's')
%      pars.funit string to label the second index of the input data, which
%                  corrsponds to the y-axis in the plots ((default is 'Hz') 
%     

mytext='';
xy=locs(:,2:3);[xdum ishowaxis]=min(xy(:,1)+xy(:,2));
if length(size(data))==2;[nchan nt]=size(data);data=reshape(data,nchan,1,nt);end
[nchan nf nt]=size(data);
xf=(1:nt)';
yf=(1:nf)';
tunit='s';
funit='Hz';
kplot=0;
      scal=[min(min(min(abs(data)))),max(max(max(abs(data))))];
 %     dr=reshape(data,[],1);drs=sort(dr);nnd=length(dr);
 %       scal=[drs(ceil(.01*nnd)),drs(ceil(.99*nnd))]; 
         %scal=[drs(1),drs(end)]; 
    resolution=25;
    global_size=1;
    rel_size=1;
    head_right=0;
    head_up=0;

    
if nargin>2
    if isfield(pars,'scale')
        scal=pars.scale;
        if length(scal)==1
            scal=[-scal scal];
        end
    end

    if isfield(pars,'global_size');
        global_size=pars.global_size;
 
    end
    if isfield(pars,'relative_size');
        rel_size=pars.relative_size;
 
    end
    if isfield(pars,'resolution')
        resolution=pars.resolution;
   
    end

    if isfield(pars,'head_right')
        head_right=pars.head_right;
   
    end
    if isfield(pars,'head_up')
        head_up=pars.head_up;

    end
    if isfield(pars,'mytext')
        mytext=pars.mytext;
 
    end
    if isfield(pars,'freqaxis')
        yf=pars.freqaxis;
    end
    if isfield(pars,'timeaxis')
        xf=pars.timeaxis;
    end
     if isfield(pars,'tunit');
        tunit=pars.tunit;
     end
      if isfield(pars,'funit');
        funit=pars.funit;
     end
    if isfield(pars,'plot');
        kplot=pars.plot;
     end
else
   
end



[n,m]=size(locs);
locs_all=locs;
[nall,ndum]=size(locs_all);
ind2chan=(1:nall)';
[nc,nc]=size(data);

if m==5
  if nc>n
    ind2chan=abs(locs(:,1));
  end
  indd=ind2chan(locs(:,1)>.5);
  locs=locs(locs(:,1)>.5,:);
end


[n,m]=size(locs);
if m==2;
    locs=[(1:n)',locs,locs];
elseif m==3
    locs=[locs,locs(:,2:3)];
elseif m==4
    locs=[(1:n)',locs];
end
%ind2chan=locs(:,1);

locs(:,5)=detrend(locs(:,5),0);
%locs(:,4)=detrend(locs(:,4),'constant');


coor=locs(:,4:5);

ymin=min(coor(:,1));
ymax=max(coor(:,1));
xmin=max(coor(:,1));
xmax=max(coor(:,1));


%minmin=mindis(coor(:,1),coor(:,2));
[minmin,minmax,meanmin]=mindis(coor(:,1),coor(:,2));

tot_scale=max(abs([xmin,xmax,ymin,ymax]))+minmin/2;
tot_scale=tot_scale/1.1;
locs(:,2:5)=locs(:,2:5)*.5/tot_scale;
meanmin=meanmin*.5/tot_scale;


no=zeros(n,1);
for i=1:n
    no(i)=norm(coor(i,:));
end
faktor=1/max(no)/2/1.3;
%wfaktor=faktor*1.3
wfaktor=meanmin/1.6;
wfaktor=wfaktor/1.3;
wfaktor=wfaktor*global_size;


meanx=mean(faktor*locs(:,4)*.85+.45);
meany=mean(faktor*locs(:,5)*.85+.45);

rad=max(sqrt( (faktor*locs(:,4)*.85+.45-meanx).^2+(faktor*locs(:,5)*.85+.45-meany).^2 ));
cirx=(rad*cos((1:1000)*2*pi/1000)+meanx)';ciry=(rad*sin((1:1000)*2*pi/1000)+meany)';
   faktor=.75;
faktor=faktor*global_size*rel_size;

nplot=n; figure_w
for chan=1:n
    zp=squeeze(data(indd(chan),:,:));
    %zp(1,1)=1000,
      %zp=z(ind2chan);
  
  
  
  
 subplot(nplot,nplot,n*chan);
 set(gca,'Position',[faktor*locs(chan,2)*.85+.45 faktor*(locs(chan,3)*1.-.05)+.45 wfaktor*0.070 wfaktor*0.088]);
%    set(gca,'Position',[faktor*locs(chan,4)*.85+.45 faktor*(locs(chan,5)*1.+.01)+.45 wfaktor wfaktor*.088/.07]);
%  axes('Position',[faktor*locs(chan,4)*.85+.45 faktor*(locs(chan,5)*1.+.01)+.45 wfaktor wfaktor*.088/.07]);
 
%   plot_elec_empty_lowres(zp,locs_all(:,2:3),scal,indd(chan),resolution);
  imagesc(xf,yf,zp,[-0.05 0.05]); axis off
  
  if kplot==1
      plot(xf,zp');
      axis([ min(xf) max(xf) min(min(zp)) max(max(zp))])
  else
      
  [X,Y] = meshgrid(xf,yf);
  surface(X,Y,zeros(size(zp)),zp,'edgecolor','none');shading interp;
  axis([ min(xf) max(xf) min(yf) max(yf)])
  end
  
  if chan == ishowaxis
      set(gca,'fontweight','bold','fontsize',8);
      xlabel(tunit,'fontweight','bold');
      ylabel(funit,'fontweight','bold');
  else
  %    if kplot==0;
   set(gca,'visible','off');
   %   end
   set(gca,'xTick',[])
   set(gca,'yTick',[])
  end
  caxis('manual');
  %scal=scal
  caxis(scal);
  if chan==1
      %set(gca,'Position',[faktor*locs(chan,4)*.85+.45 faktor*(locs(chan,5)*1.-.05)+.45 wfaktor*0.08 wfaktor*0.08]);
       set(gca,'Position',[faktor*locs(chan,4)*.85+.45 faktor*(locs(chan,5)*1.+.01)+.45 wfaktor*1.1 wfaktor*.088/.07]);

       if kplot==0;
      h1=colorbar;
      set(h1,'position',[.85 .2 .03 .6]);
      set(h1,'fontweight','bold','fontsize',20);
       end
  end
end;
caxis('manual');
caxis(scal);
subplot(nplot,nplot,1);
set(gca,'Position',[0.12 0.01 .7/.88 .99]);
%set(gca,'Position',[0.12 0.01 1 1]);

scalx=1.0;




% drawhead(meanx+wfaktor/4+head_right,.45+wfaktor/4+head_up,.43,scalx);
axis([ 0 1 0 1]);
if length(mytext)>0
    text(.0, .9,mytext,'fontsize',12,'fontweight','bold');
end
%text(.1,.85,figtitle,'HorizontalAlignment','center');
%title('Titel')
set(gca,'visible','off');
%get(h1)
%h=colorbar; 
%set(h,'yticklabel',
colormap jet
return;

function plot_elec_empty_lowres(z,loc,skal,chan,resolution); 

x=loc(:,1);
y=loc(:,2);

xlin = linspace(1.4*min(x),1.4*max(x),resolution);
ylin = linspace(1.4*min(y),1.4*max(y),resolution);
[X,Y] = meshgrid(xlin,ylin);
% Z = griddata(x,y,z,X,Y,'invdist');
Z = griddata(x,y,z,X,Y,'nearest');


  % Take data within head
  rmax=1.1*max(sqrt(x.^2+y.^2));
  mask = (sqrt(X.^2+Y.^2) <= rmax);
  ii = find(mask == 0);
  Z(ii) = NaN;
  
  
surface(X,Y,zeros(size(Z)),Z,'edgecolor','none');shading interp;
%caxis([ - max(abs(z)) max(abs(z))]);
%disp([ - max(abs(z)) max(abs(z))]);
%caxis([ -skal  skal]);
%disp([ -skal  skal]);

hold on;
%plot(x,y,'.k');
axis([-rmax rmax -rmax rmax]);
%colorbar;
plot(loc(chan,1),loc(chan,2),'.k');
%set(gcf,'color','none');

set(gca,'xTick',[])
set(gca,'yTick',[])
plot(.985*rmax*sin((0:1000)/1000*2*pi), .985*rmax*cos((0:1000)/1000*2*pi),'linewidth',2,'color','k'); 
%set(gcf,'color','none');

return; 

function drawhead(x,y,size,scalx);

cirx=(x+scalx*size*cos((1:1000)*2*pi/1000) )';ciry=(y+size*sin((1:1000)*2*pi/1000))';

plot(cirx,ciry,'k','linewidth',1);
hold on;

ndiff=20;
plot( [x  cirx(250-ndiff) ],[y+1.1*size ciry(250-ndiff)],'k','linewidth',1);
plot( [x  cirx(250+ndiff) ],[y+1.1*size ciry(250+ndiff)],'k','linewidth',1);


return;

function [minmin,minmax,meanmin]=mindis(x,y);

[n,m]=size(x);

minall=zeros(n,1);



for i=1:n
    md=1.e8;
    for j=1:n
        if j~=i
          %dis=max( abs(x(i)-x(j)),abs(y(i)-y(j)));
          dis=norm([x(i)-x(j),y(i)-y(j)]);
          if dis<md
             md=dis;
          end
        end     
    end
    minall(i)=md;
end



minmin=min(minall);
meanmin=mean(minall);
[minmax,imax]=max(minall);

return;

