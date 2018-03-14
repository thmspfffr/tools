function h=showvc_prog_plain(loc,tri,view_dir,para);

cmap=1;
voxelkont=0;
mycolormap='jet';
colorbars=1;
if nargin>3
  if isfield(para,'cmap')
     cmap=para.cmap;
  end
  if isfield(para,'voxelfield')
     voxelfield=para.voxelfield;
     voxelkont=1;
  end
 if isfield(para,'mycolormap')
     mycolormap=para.mycolormap;
 end
 if isfield(para,'colorbars')
     colorbars=para.colorbars;
 end

end




[ntri,ndum]=size(tri);


%colormap(newmap(30:150,:));
% idx = find(loc(:,1)>mean([-7.1200 7.1268]));
% loc(idx,:) = [];




locm=mean(loc);
[nloc,ndum]=size(loc);
relloc=loc-repmat(locm,nloc,1);
dis=(sqrt(sum((relloc.^2)')))';thresh=3;dis(dis<thresh)=thresh;

cortexcolor=[.6 .6 .6];
h=patch('vertices',loc,'faces',tri);
set(h,'FaceColor',cortexcolor);
view(view_dir);
set(h,'edgecolor','none');
set(h,'facelighting','phong');
%set(h,'specularexponent',50);
set(h,'specularstrength',0.5);
set(h,'ambientstrength',0.5);
set(h,'diffusestrength',.5)
set(h,'SpecularColorReflectance',0)

dis0=0*dis+0;
set(h,'facevertexalphadata',dis0)
set(h,'alphadatamapping','direct')
%set(h,'facealpha',.1)
camlight('headlight','infinite');
axis equal 

if voxelkont>0
   vmax=max(voxelfield);vmin=min(voxelfield);
     if vmax <=0 
         vv=vmin;
     else
         vv=vmax;
     end
   tri_new=[];
   nt=length(tri);
   vfx=voxelfield;
     for i=1:nt;
       %xx=0;for j=1:3; xx=xx+norm(voxelfield(tri(i,j)));end;
       xx1=min(abs(voxelfield(tri(i,1:3))));
       xx2=mean(abs(voxelfield(tri(i,1:3))));
       if xx1>0 
         tri_new=[tri_new;tri(i,:)];
       end
       if xx1 ==0 
         vfx(tri(i,1:3))=NaN;
       end 
     end
%      vf=voxelfield(abs(voxelfield)>=0);
%      map=colormap(mycolormap);
%      newmap=colormap_interpol(map,3);size(newmap);colormap(newmap);
%      [nc,nx]=size(newmap);
%      vfint=ceil((nc-1)*((voxelfield-min(vf))/(max(vf)-min(vf))));
%      vfint(voxelfield==0)=1;
%      vftruecolor=newmap(vfint);
   %h=patch('vertices',loc,'faces',tri_new,'FaceVertexCData',voxelfield,...
   % 'facecolor','interp','edgecolor','none','facelighting','phong');
%     h=patch('vertices',loc,'faces',tri_new,'FaceVertexCData',vftruecolor,...
%     'facecolor','interp','edgecolor','none','facelighting','phong');
    h=patch('vertices',loc,'faces',tri_new,'FaceVertexCData',vfx,...
    'facecolor','interp','edgecolor','none','facelighting','phong');
    caxis(para.colorlimits);
%     load colormapsmooth4
   map=colormap(jet);
   colorbars = 0;
   newmap=colormap_interpol(map,3);size(newmap);colormap(newmap);
   
   if colorbars==1  
       h1=colorbar;
       PP=get(h1,'position');
       PP=[PP(1)+.05 PP(2)+PP(4)/4 PP(3)/2 PP(4)/2];
       set(h1,'position',PP);
   end
       
end

set(h,'specularexponent',50000);

return;