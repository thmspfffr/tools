function [s_all,vmax_all,imax_all,dips_mom_all,dips_loc_all]=sc_music(patt,V,ns,grid,para);
% calculates ns distributions from a given pattern( e.g. data subspace)
% using the SC-MUSIC approach (Shahbazi et al., Self-Consistent MUSIC: an approach to the
% localization of true brain interactions from EEG/MEG data, Neuroimage, 2015),
% i.e. each distribution is calcutaled as a MUSIC-scan with
% best dipole from other distributions projected out. The code includes MUSIC and RAP-MUSIC
% (Mosher et al. 1999) as a special case.
%
%
% input:
%        patt  nxp matrix for n channels; This matrix contains the p
%              dimensional subspace of the data.
%              each column in patt represents a spatial pattern;
%             (only the span(patt) matters; mixing of the patterns has no effect)
%        V     nxmx3 matrix of forward models for m grid_points; V(:,i,j) is the
%              potential of a unit dipole at point i in direction j;
%        ns    number of sources to be localized
%        grid  (optional) mx3 matrix denoting locations of grid points.
%              If you omit that the output dip_loc will be empty.
%        para  (optional structure)
%              para.nite (integer number) Determines the number of iterations in the algorithm.
%                 (Default: nite=100).
%                 para.nite=0 performs RAP-MUSIC scans, and the first of these scans is the
%                 the MUSIC scan. And any value larger than zero performs SC-MUSIC.
%              para.ranstart: (inter number) (Default: ranstart=0). If para.ranstart=0 the first
%                 scan corresponds to a 'classical' MUSIC scan. If para.ranstart=1 the first
%                 topography to be projected out in the following
%                 corresponds to a dipole of random location and orientation.
%              para.v0: (Nx1 vector for N sensors) v0 is a specified topography to be
%                 projected out in the first scan.
%
%  output:
%          s_all mxns matrix; s_all(i,j) indicates fit-quality (from 0 (worst) to 1 (best)) at grid-point
%                        i of the j.th distribution. More precisely,
%                        s_all(i,j)=cos(phi) where phi is the principal
%                        angle between the signal space,patt, and the Nx3
%                        subspace defined by 3 dipoles at the i.th location.
%          vmax_all  nxns matrix; vmax_all(:,i) is the field of the best dipole of the i.th distribution
%          imax_all  ns numbers imax_all(i) denotes grid-index of best dipole of i.th  distribution
%          dips_mom_all  nsx3 matrix dips_mom_all(i,:) is the moment of the  best
%                              dipole of the i.th distribution
%          dips_loc_all nsx3 matrix dips_loc_all(i,:) is the location of the  best
%                              dipole of the i.th distribution
%


% License
%   Copyright (C) 2015 
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see http://www.gnu.org/licenses/.



if nargin<5
    para=[];
end

nite=100;
if isfield(para,'nite')
    nite=para.nite;
end


ranstart=0;
if isfield(para,'ranstart')
    ranstart=para.ranstart;
end
if isfield(para,'v0')
    ranstart=2;
end



data=orth(patt);
[nchan,ng,ndum]=size(V);
spacecorr=zeros(ng,1);
s_all=zeros(ng,ns);
vmax_all=zeros(nchan,ns);
imax_all=zeros(ns,1);

if ranstart==0;
    [s_all(:,1),vmax_all(:,1),imax_all(1)]=one_ite(V,data);
    for i=2:ns
        proj_pat=vmax_all(:,1:i-1);
        [s_all(:,i),vmax_all(:,i),imax_all(i)]=one_ite(V,data,proj_pat);
    end
    
elseif ranstart==1
    [nchan ng ndipdir]=size(V);
    for i=1:ns
        iran=ceil((ng-1)*rand(1,1));
        imax_all(i)=iran;
        mom=randn(ndipdir,1);
        v0=squeeze(V(:,iran,:))*mom;
        vmax_all(:,i)=v0/norm(v0);
    end
    
    for k=1:100
        [vmax_new,reldiff]=vmax2vmax_new(V(:,imax_all,:),vmax_all,patt);
        vmax_all=vmax_new;
    end
    %disp([k,100*reldiff]);
elseif ranstart==2
    v0=para.v0;
    [nchan,nsstart]=size(v0);
    if nsstart==1
        v0=v0/norm(v0);
        vmax_all(:,1)=v0;
        proj_pat=vmax_all(:,1);
        [s_all(:,2),vmax_all(:,2),imax_all(2)]=one_ite(V,data,proj_pat);
    elseif nsstart==2
        v0=[v0(:,1)/norm(v0(:,1)),v0(:,2)/norm(v0(:,2))];
        vmax_all=v0;
    end
end




vmax_new=vmax_all;

if nite>0
    loccount=0;
    cont=0;i=0;
    while cont==0
        i=i+1;
        s_old=s_all;
        imax_all_old=imax_all;
        for ii=1:ns
            proj_pat=vmax_new(:,[1:ii-1,ii+1:ns]);
            [s_all(:,ii),vmax_new(:,ii),imax_all(ii)]=one_ite(V,data,proj_pat);
        end
        
        
        if norm(imax_all_old-imax_all)<.001
            loccount=loccount+1;
        else
            loccount=0;
        end
        
        vmax_all=vmax_new;
        reldiff=1;
        nndir=50+round(10*rand-5);
        for k=1:49
            %lambda=0;
            %imax_all
            [vmax_new,reldiff]=vmax2vmax_new(V(:,imax_all,:),vmax_all,patt);
            %disp([k,100*reldiff]);
            vmax_all=vmax_new;
        end
        %rel=norm(1./(1-s_old)-1./(1-s_all),'fro')/norm(1./(1-s_all),'fro');
        
        rel=norm(s_old-s_all,'fro')/norm(s_old+sqrt(eps));
        %rel=norm(1./(1-s_old)-1./(1-s_all),'fro')/norm(1./(1-s_all),'fro');
        %disp([i,rel*100,reldiff*100]);
        %disp(imax_all');
        if i>nite-1 | rel<1e-4 | loccount>2
            cont=1;
        end
        
    end
    
end

dips_mom_all=vmax2dipmom(V,imax_all,vmax_all);

if nargin>3
    if length(grid)>0;
        dips_loc_all=grid(imax_all,:);
    else
        dips_loc_all=[];
    end
else
    dips_loc_all=[];
end

return;

function [s,vmax,imax]=one_ite(V,data,proj_pat)
if nargin==2;
    [nchan,ng,ndum]=size(V);
    s=zeros(ng,1);
    for i=1:ng;
        Vortholoc=orth(squeeze(V(:,i,:)));
        s(i)=calc_spacecorr(Vortholoc,data);
    end
    [smax,imax]=max(s(:));
    Vortholoc=orth(squeeze(V(:,imax,:)));
    vmax=calc_bestdir(Vortholoc,data);
else
    [nchan,ng,ndum]=size(V);
    s=zeros(ng,1);
    proj_pat=orth(proj_pat);
    data_proj=orth(data-proj_pat*(proj_pat'*data));
    for i=1:ng;
        Vortholoc=orth(squeeze(V(:,i,:)));
        s(i)=calc_spacecorr(Vortholoc,data_proj,proj_pat);
    end
    [smax,imax]=max(s);
    Vortholoc=orth(squeeze(V(:,imax,:)));
    vmax=calc_bestdir(Vortholoc,data_proj,proj_pat);
end
return


function s=calc_spacecorr(Vloc,data_pats,proj_pats)
if nargin==2
    A=data_pats'*Vloc;
    [u s v]=svd(A);
    s=s(1,1);
else
    V_proj=orth(Vloc-proj_pats*(proj_pats'*Vloc));
    data_proj=orth(data_pats-proj_pats*(proj_pats'*data_pats));
    A=data_proj'*V_proj;
    [u s v]=svd(A);
    s=s(1,1);
end
return;

function [vmax,s]=calc_bestdir(Vloc,data_pats,proj_pats)
if nargin==2
    A=data_pats'*Vloc;
    [u s v]=svd(A);
    vmax=Vloc*v(:,1);
    vmax=vmax/norm(vmax);
    s=s(1,1);
else
    [n m]=size(Vloc);
    V_proj=orth(Vloc-proj_pats*(proj_pats'*Vloc));
    A=data_pats'*V_proj;
    [u s v]=svd(A);
    BB=(Vloc'*proj_pats);
    Q=inv(eye(m)-BB*BB'+sqrt(eps));
    vmax=Vloc*(Q*Vloc'*(V_proj*v(:,1)));
    vmax=vmax/norm(vmax);
    s=s(1,1);
end
return;

function u=myortho(V);
[n m]=size(V);
[u,s,v]=svd(V);
u=u(:,1:m);
return


function   dips_mom_all=vmax2dipmom(V,imax_all,vmax_all);
ns=length(imax_all);
dips_mom_all=zeros(ns,3);
for i=1:ns
    Vloc=squeeze(V(:,imax_all(i),:));
    v=vmax_all(:,i);
    dip=inv(Vloc'*Vloc)*Vloc'*v;
    dips_mom_all(i,:)=dip'/norm(dip);
end

return;

function [vmax_new,reldiff]=vmax2vmax_new(Vs,vmax,data);

[nchan,ns,ndipdir]=size(Vs);
data_pats=orth(data);
vmax_new=vmax;
for i=1:ns
    vmax_loc=vmax(:,[1:i-1,i+1:end]);
    proj_pats=myorth(vmax_loc);
    Vloc=orth(squeeze(Vs(:,i,:)));
    V_proj=orth(Vloc-proj_pats*(proj_pats'*Vloc));
    data_pats=orth(data-proj_pats*(proj_pats'*data));
    A=data_pats'*V_proj;
    [u s v]=svd(A);
    BB=(Vloc'*proj_pats);
    Q=inv(eye(ndipdir)-BB*BB'+sqrt(eps));
    vmax_x=Vloc*(Q*Vloc'*(V_proj*v(:,1)));
    vmax_x=vmax_x/norm(vmax_x);
    %vmax_new(:,i)=vmax_x;
    vmax_new(:,i)=.9*vmax_x+.1*vmax(:,i);
end

reldiff=norm(vmax-vmax_new,'fro')/norm(vmax,'fro');

%save vmax vmax vmax_new

return;

function v0=find0(V,data,nrun);
[nchan,ng,ndipdir]=size(V);
s=zeros(nrun,1);
iran=ceil((ng-1)*rand(nrun,2));
data_orth=orth(data);
for i=1:nrun;
    Vloc=[squeeze(V(:,iran(i,1),:)),squeeze(V(:,iran(i,2),:))];
    Vloc_orth=orth(Vloc);
    A=Vloc_orth'*data_orth;
    [u,ss,v]=svd(A);
    s(i)=ss(1,1);
end
[smax,imax]=max(s);
Vloc=[squeeze(V(:,iran(imax,1),:)),squeeze(V(:,iran(imax,2),:))];
Vloc_orth=orth(Vloc);
A=data_orth'*Vloc_orth;
[u, s ,v]=svd(A);
vmax=Vloc_orth*v(:,1);
dip=inv(Vloc'*Vloc)*(Vloc'*vmax);
v1=Vloc*[dip(1:3);zeros(3,1)];
v2=Vloc*[zeros(3,1);dip(4:6)];
if norm(v1)>norm(v2)
    v0=[v1,v2];
else
    v0=[v2,v1];
end

return


function aout=myorth(ain);
[n,m]=size(ain);
[u,s,v]=svd(ain);
aout=u(:,1:m);
return




