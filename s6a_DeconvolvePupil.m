% after s1_PupilAnalysisMain, this script does deconvolution on the resulting file

clear all; close all; clc;
addpath(genpath('~/Dropbox/code/Commitment/'));
cd('~/Data/Commitment');

subjects = [0:5 7:15];
%subjects = 0;

for sj = unique(subjects),
    clc; disp(sj);
    
    clearvars -except sj subjects;
    load(sprintf('~/Data/Commitment/Pupil/P%02d/P%02d_alleye_filtered.mat', sj, sj));
    pupilchan = find(strcmp(data.label, 'EyePupil')==1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % estimate IRF for each event
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    cfg                         = [];
    cfg.regressors              = {'blink', 'int1', 'int2', 'estimation', 'feedback', 'choice1', 'choice0', 'choice-1'};
    cfg.impulse                 = 6; % in seconds;
    tic; designM                = makeDeconvDesignM(cfg, data); toc;
    tic; deconvolvedPupil       = pinv(designM) * (data.trial{1}(pupilchan, :))'; toc; % pinv more robust than inv?
    
    % reshape into the separate events that were estimated
    IRF.mean = reshape(deconvolvedPupil, [cfg.impulse*data.fsample+1 length(cfg.regressors)]);
    IRF.name = cfg.regressors;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % predict nuisance events and take out
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    cfg                         = [];
    cfg.regressors              = {'blink', 'int1', 'int2', 'estimation', 'feedback'};
    cfg.impulse                 = 6; % in seconds
    designM                     = makeDeconvDesignM(cfg, data);
    tic; deconvolvedPupil       = pinv(designM) * (data.trial{1}(pupilchan, :))'; toc; % pinv more robust than inv?
    
    % reshape into the separate events that were estimated
    IRF.predic = reshape(deconvolvedPupil, [cfg.impulse*data.fsample+1 length(cfg.regressors)]);
    
    % compute predicted response
    predic      = designM * deconvolvedPupil;
    
    % subtract the mean of both to not get weird offsets
    predic                         = predic - mean(predic);
    data.trial{1}(pupilchan, :)    = data.trial{1}(pupilchan, :) - mean(data.trial{1}(pupilchan, :));
    
    % take the predicted nuisances out
    residuals   = data.trial{1}(pupilchan, :) - predic';
    
    % plot
    clf; plot(data.trial{1}(pupilchan, :), 'b');
    hold on; plot(predic, 'r'); plot(residuals', 'g');
    
    % apply redefinetrial back to the residuals
    data.trial{1}(pupilchan, :)   = residuals;
    savefast(sprintf('~/Data/Commitment/Pupil/P%02d/P%02d_alleye_residuals.mat', sj, sj), 'data', 'IRF');
end