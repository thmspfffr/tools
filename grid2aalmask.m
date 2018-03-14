function [mask,labels, ac] = tp_aalgrid()
%function mask = grid2aalmask(atlas, grid)
%
% Looks up for each coordinate in grid which AAL regions a voxel is located
%
% IN   altas     - anatomical atlas, e.g. fom fieldtrip:
%                         aal =  
%                             dim: [91 109 91]
%                             hdr: [1x1 struct]
%                       transform: [4x4 double]
%                            unit: 'mm'
%                          tissue: [91x109x91 double]
%                     tissuelabel: {1x116 cell}
%                        coordsys: 'mni'
%      grid      - [number of voxels x 3] (xzy coordinates in cm!)
%
% OUT  mask
%      labels
%
%
% by Arne Ewald, 13.1.2015

fprintf('Calculating atlas region mask. This might take a while....\n');

aal = ft_read_atlas('~/Documents/MATLAB/fieldtrip-20160919/template/atlas/aal/ROI_MNI_V4.nii');
% get coordinates from voxel indices
ac=zeros(aal.dim(1)*aal.dim(2)*aal.dim(3),3);
%ac2=zeros(atlas.dim(1)*atlas.dim(2)*atlas.dim(3),3);

% ix=18;iy=20;iz=78;
ii=0;
for ix=1:aal.dim(1)
    for iy=1:aal.dim(2)
        for iz=1:aal.dim(3)
            ii=ii+1;
            r=[ix,iy,iz,1]';
            rc=aal.transform*r;
            ac(ii,:)=rc(1:3)';  
            % same thig using FieldTrip
            %ac2(ii,:)=ft_warp_apply(atlas.transform, [ix iy iz]);
        end
    end
end

% subsample
ii=0;
for ix=1:3:aal.dim(1)
    for iy=1:3:aal.dim(2)
        for iz=1:3:aal.dim(3)
            ii=ii+1;
            r=[ix,iy,iz,1]';
            rc=aal.transform*r;
            bc(ii,:)=rc(1:3)'; 
            aal.tissue_4mm(ii) = aal.tissue(ix,iy,iz);
            % same thig using FieldTrip
            %ac2(ii,:)=ft_warp_apply(atlas.transform, [ix iy iz]);
        end
    end
end


