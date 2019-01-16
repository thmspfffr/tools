function sa_out=tp_mk_sa_meg_withoutmri(sa,fn,para)
% Creates sa-strucures for MEG without MRI. The individual volume conductor
% is estimated from fiducials as a transformed template volume conductor.
% It uses fieldtrip functions to read locations of the fiducials and sensors.
% Forward calculation is done using expansions in spherical harmonics.
% First, to calculate an analytic surface (and normals) from surface
% points, and second to expand the lead field. In this code, the order are
% set to be 15 and 10, respectively. Change the parameters p1 and p2 inside
% this code  if you prefer different values.
%
% Usage: sa_out=mk_sa_meg_withoutmri(sa,fn,para);
%
% input:
% sa: template sa structure
% fn: name of directory containing CTF data
% para: optional structure. if para.cont_locs_2D is different from 1, then
%       the 2D sensor locations are not calculated. (This is included because
%       the sensor repositioning for head-in-head plots is rather slow.
%
%
% output:
% sa_out: final structure with complete structural information


cont_locs_2D=1;
nin=[];
if nargin>2
  if isfield(para,'cont_locs_2D');
    cont_locs_2D=para.cont_locs_2D;
  end
  if isfield(para,'nin');
    nin=para.nin;
  end
end

p1=15;
p2=10;

sa_out=sa;
[shape] = ft_read_headshape(fn);

x1=sa.naspalparini;
x1=x1(1:3,:);
if isfield(shape.fid,'pnt')
  y1=shape.fid.pnt;
elseif  isfield(shape.fid,'pos')
  y1=shape.fid.pos;
else
  error('no fiducials found in shape.fid')
end
x1m=mean(x1);
y1m=mean(y1);
x2=x1-repmat(x1m,3,1);
y2=y1-repmat(y1m,3,1);

Ayx=pinv(y2'*y2)*y2'*x2;
[u,s,v]=svd(Ayx);
scal=(s(1,1)+s(2,2))/2;
s(3,3)=scal;
Ayx=u*s*v';
ryx=x1m-y1m*Ayx;


sa_out.trafo.u_indi2template=Ayx;
sa_out.trafo.r_indi2template=ryx;
sa_out.trafo.readme='xtemplate=xindi*u+r, row-vectors';


vc=sa.vc{1}.vc;
[nvc,ndum]=size(vc);
vc_indi=(vc-repmat(ryx,nvc,1))*inv(Ayx);

sa_out.vc_indi.vc=vc_indi;


sens = ft_read_sens(fn);
n=length(sens.label);
inds=[];
%save sens sens
for i=1:n;
  x=sens.label{i};
  if x(1)=='M';
    inds=[inds,i];
  end
end
nchan=length(inds);
sens1=[sens.coilpos(inds,:),sens.coilori(inds,:)];
refs=[sens.coilpos(inds+nchan+min(inds)-1,:),-sens.coilori(inds+nchan+min(inds)-1,:)];

[vc_indi,center,radius,coeffs]=pointsonsurface(vc_indi,vc_indi,p1);

sa_out.fp_indi=meg_ini(vc_indi,center',p2,sens1,refs);

sa_out.locs_3D_indi=sens.coilpos(inds,:);
sa_out.sens_indi=sens;
sa_out.coils_indi=sens1;
sa_out.inds=inds;

if cont_locs_2D==1;
  if length(nin)==0;
    para1.rot=90;
    locs_2D=mk_sensors_plane(sa_out.locs_3D_indi,para1);
    para2.rot=90;
    para2.nin=50;
    locs_2D_sparse=mk_sensors_plane(sa_out.locs_3D_indi,para2);
    sa_out.locs_2D=locs_2D;
    sa_out.locs_2D_sparse=locs_2D_sparse;
  else
    para1.showfigs=0;
    para1.rot=90;
    para1.nin=nin;
    locs_2D=mk_sensors_plane(sa_out.locs_3D_indi,para1);
    sa_out.locs_2D=locs_2D;
  end
  
end

if isfield(sa,'grid_xcoarse')
  gg=sa.grid_xcoarse;
  [ng,ndum]=size(gg);
  sa_out.grid_xcoarse_indi=(gg-repmat(ryx,ng,1))*inv(Ayx);
end
if isfield(sa,'grid_coarse')
  gg=sa.grid_coarse;
  [ng,ndum]=size(gg);
  sa_out.grid_coarse_indi=(gg-repmat(ryx,ng,1))*inv(Ayx);
end
if isfield(sa,'grid_medium')
  gg=sa.grid_medium;
  [ng,ndum]=size(gg);
  sa_out.grid_medium_indi=(gg-repmat(ryx,ng,1))*inv(Ayx);
  sa_out.L_medium=grid2L(sa_out.grid_medium_indi,sa_out.fp_indi);
end
if isfield(sa,'grid_fine')
  gg=sa.grid_fine;
  [ng,ndum]=size(gg);
  sa_out.grid_fine_indi=(gg-repmat(ryx,ng,1))*inv(Ayx);
end
if isfield(sa,'grid_cortexhippo')
  gg=sa.grid_cortexhippo;
  [ng,ndum]=size(gg);
  sa_out.grid_cortexhippo_indi=(gg-repmat(ryx,ng,1))*inv(Ayx);
end
if isfield(sa,'grid_cortex3000')
  gg=sa.grid_cortex3000;
  [ng,ndum]=size(gg);
  sa_out.grid_cortex3000_indi=(gg-repmat(ryx,ng,1))*inv(Ayx);
end
if isfield(sa,'grid_aal4mm')
  gg=sa.grid_aal4mm;
  [ng,ndum]=size(gg);
  sa_out.grid_aal4mm_indi=(gg-repmat(ryx,ng,1))*inv(Ayx);
end
if isfield(sa,'grid_aal4mm')
  gg=sa.grid_aal6mm;
  [ng,ndum]=size(gg);
  sa_out.grid_aal6mm_indi=(gg-repmat(ryx,ng,1))*inv(Ayx);
end
if isfield(sa,'grid_m758_4mm')
  gg=sa.grid_m758_4mm;
  [ng,ndum]=size(gg);
  sa_out.grid_m758_4mm_indi=(gg-repmat(ryx,ng,1))*inv(Ayx);
end
if isfield(sa,'grid_m758_6mm')
  gg=sa.grid_m758_6mm;
  [ng,ndum]=size(gg);
  sa_out.grid_m758_6mm_indi=(gg-repmat(ryx,ng,1))*inv(Ayx);
end
if isfield(sa,'grid_vtpm_4mm')
  gg=sa.grid_vtpm_4mm;
  [ng,ndum]=size(gg);
  sa_out.grid_vtpm_4mm_indi=(gg-repmat(ryx,ng,1))*inv(Ayx);
end
if isfield(sa,'grid_vtpm_6mm')
  gg=sa.grid_vtpm_6mm;
  [ng,ndum]=size(gg);
  sa_out.grid_vtpm_6mm_indi=(gg-repmat(ryx,ng,1))*inv(Ayx);
end
if isfield(sa,'grid_cortex_brainstorm')
  gg=sa.grid_cortex_brainstorm;
  [ng,ndum]=size(gg);
  sa_out.grid_cortex_brainstorm_indi=(gg-repmat(ryx,ng,1))*inv(Ayx);
end
if isfield(sa,'grid_genemaps')
  gg=sa.grid_genemaps;
  [ng,ndum]=size(gg);
  sa_out.grid_genemaps_indi=(gg-repmat(ryx,ng,1))*inv(Ayx);
end
if isfield(sa,'grid_genemaps_aal')
  gg=sa.grid_genemaps_aal;
  [ng,ndum]=size(gg);
  sa_out.grid_genemaps_aal_indi=(gg-repmat(ryx,ng,1))*inv(Ayx);
end
if isfield(sa,'grid_MM')
 % MULTIMODAL PARCELLATION (MM): Glasser et al.
  gg=sa.grid_MM;
  [ng,ndum]=size(gg);
  sa_out.grid_MM_indi=(gg-repmat(ryx,ng,1))*inv(Ayx);
end

[sa.grid_cortex_lowres] = select_chans(sa_out.grid_cortex3000,400);
gg=sa.grid_cortex_lowres;
[ng,ndum]=size(gg);
sa_out.grid_cortex_lowres_indi=(gg-repmat(ryx,ng,1))*inv(Ayx);

[sa.grid_cortex800] = select_chans(sa_out.grid_cortex3000,800);
gg=sa.grid_cortex800;
[ng,ndum]=size(gg);
sa_out.grid_cortex800_indi=(gg-repmat(ryx,ng,1))*inv(Ayx);



return;



