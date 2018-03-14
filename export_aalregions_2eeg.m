%% Init


%% TO DO: adapt paths

sh=filesep;

if isunix
    prefix='/home/aewald/';
else
    prefix='C:\';
end

addpath(([prefix 'Toolboxes' sh 'HHToolbox' sh]));
startup_hhtb;

% adding Fieldtrip
addpath(([prefix 'Toolboxes' sh 'Fieldtrip' sh 'fieldtrip-20150113' sh]));
addpath(([prefix 'Toolboxes' sh 'Fieldtrip' sh 'fieldtrip-20150113' sh 'fileio']));
addpath(([prefix 'Toolboxes' sh 'Fieldtrip' sh 'fieldtrip-20150113' sh 'utilities']));

ft_defaults

% where to save the results
pth=['C:\Code\mycode\Jan_Adip\'];

%% load MNI template

load('sa_eeg_template');
sa=sa_eeg_template;
clear sa_eeg_template

%% load AAL regions

aal = ft_read_atlas('ROI_MNI_V4.nii');
%aal.transform(1,:)=-aal.transform(1,:);
aal.tissue=permute(aal.tissue, [3,2,1]);
% tissue2=zeros(size(aal.tissue));
% for i1=1:aal.dim(1)
%     for i3=1:aal.dim(3)
%         tissue2(i1,:,i3)=-aal.tissue(i1,:,i3);
%     end
% end
% aal.tissue=tissue2;

fname=['aal_tissuelabel_eeg2'];
tl=aal.tissuelabel;
save([pth fname], 'tl');
fprintf('Successfully saved: %s.\n', [pth fname]);


%% define grids

allgrids=[];
allgrids{1}=sa.grid_xcoarse;
allgrids{2}=sa.grid_coarse;
allgrids{3}=sa.grid_medium;
allgrids{4}=sa.grid_fine;
allgrids{5}=sa.grid_cortex3000;
allgrids{6}=sa.cortex.vc;
allgrids{7}=sa.cortex10K.vc;

gridnames=[];
gridnames{1}='grid_xcoarse';
gridnames{2}='grid_coarse';
gridnames{3}='grid_medium';
gridnames{4}='grid_fine';
gridnames{5}='grid_cortex3000';
gridnames{6}='sa.cortex';
gridnames{7}='cortex10K';


%% generate AAL region mask for a specific grid

% for gidx=1:size(gridnames,2)

for gidx=2:2

mygrid=allgrids{gidx};

fprintf('** %d - Processing grid: %s.\n', gidx, gridnames{gidx});

tic
[mask, labels, ac]=grid2aalmask2_eeg(aal, mygrid);
toc

%% Check

% v=aal.tissue(:);
% para=[];
% para.colormaps={'cool'};
% para.orientation='sagittal';
% figure;
% showmri_transp_v31(sa.mri, para, [ac(v~=0,:)/10] );
% 
% 
% reg=19;
% aal.tissuelabel(reg)
% sum(mask==reg)
% 
% para=[];
% para.colormaps={'cool'};
% para.orientation='all';
% para.mricenter=mean(ac(aal.tissue(:)==reg,:))/10;
% figure;
% showmri_transp_v31(sa.mri, para, ac(aal.tissue(:)==reg,:)/10);
% 
% para=[];
% para.colormaps={'cool'};
% para.orientation='all';
% para.mricenter=mean(mygrid(mask==reg,:));
% figure;
% showmri_transp_v31(sa.mri, para, [mygrid(mask==reg, :)]);


%% Save
aalgrid=[];
aalgrid.mask=mask;
aalgrid.labels=labels;
fname=['aalmask_eeg2_' gridnames{gidx}];
save([pth fname], 'aalgrid');
fprintf('Successfully saved: %s.\n', [pth fname]);

end

fprintf('Done.\n');


