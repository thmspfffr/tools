%% CLUSTER STATS CORTEX GRID
load cortex_brainstorm

n_subjects = 28;
clear cond1 cond2
cond1.dimord = 'pos_subj_freq_time';
% cond1.dim = [3000 1];
cond1.pos = cortex_brainstorm.pos; % das sind die mni coordinates aller voxel im cube (also n_voxel_cube X 3)
% cond1.inside = false(n_voxel_cube,1);
% cond1.inside(vox) = true; % vox sind die indices der voxel, die wirklich im grid existieren. alle anderen voxel sind 0 und damit als outside definiert
cond1.time = [1];
cond1.freq = [1];

% cond1.pow = nan(n_voxel_cube,n_subjects,1,1);
cond2 = cond1;
cfg                  = [];
% cfg.dim              = [3000 1];
% cfg.
cfg.latency          = 'all';
cfg.frequency        = 'all';
cfg.parameter        = 'pow';
cfg.method           = 'montecarlo'; %analytic montecarlo
cfg.statistic        = 'depsamplesT';
cfg.computeprob      = 'yes';
cfg.correctm         = 'cluster'; %cluster fdr no
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.clustertail      = 0; %1 = right
cfg.tail             = 0; %1 = right
cfg.alpha            = 0.025;
cfg.numrandomization = 10000;
cfg.connectivity = n;
cfg.minnbchan = 0;

cfg.avgovertime      = 'yes';
cfg.avgoverfreq      = 'yes';

design = zeros(2,2*n_subjects);
design(1,:) = repmat(1:n_subjects,1,2);
design(2,:) = mod(floor([0:(2*n_subjects-1)]/(n_subjects/1)),2)+1;
cfg.design   = design;
cfg.uvar     = 1;
cfg.ivar     = 2;

cond1.pow = squeeze(dfa_all_cnt(:,1,:));
cond2.pow = squeeze(dfa_all_res(:,1,:));
[stats_tvr1] = ft_sourcestatistics(cfg, cond1, cond2);