% for each SJ, plot an overview of raw timelock, deconv and resid

clear all; close all; clc;
addpath(genpath('~/code/Commitment/'));
cd('~/Data/Commitment');
dbstop if error;

subjects = [0:5 7:15];
warning off;

for sj = unique(subjects),
    
    fprintf('SUBJECT %d \n\n', sj);
    clearvars -except sj subjects grandavg_raw grandavg_irf grandavg_resid grandavg_nuisances;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% RAW
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    load(sprintf('~/Data/Commitment/CleanPupil/P%02d_alleye_filtered.mat', sj));
    pupilchan           = find(strcmp(data.label, 'EyePupil')==1);
    
    % epoch
    cfg         = [];
    cfg.trl     = data.trialinfo;
    lockeddata  = ft_redefinetrial(cfg, data);
    
    % plot the shiftlock for each event
    lockings(1).name    = 'blink';
    lockings(1).offset  = data.blinksmp(:,2);
    lockings(1).trials  = [];
    
    lockings(2).name    = 'int1';
    lockings(2).offset  = lockeddata.trialinfo(:,4) - lockeddata.trialinfo(:,1);
    lockings(2).prestim = 0.1;
    lockings(2).poststim = 1;
    lockings(2).trials  = [];
    
    lockings(3).name    = 'int2';
    lockings(3).offset  = lockeddata.trialinfo(:, 11) - lockeddata.trialinfo(:,1);
    lockings(3).trials  = find(~isnan(lockeddata.trialinfo(:,11)));
    lockings(3).prestim = 0.5;
    lockings(3).poststim = 1;
    
    lockings(4).name    = 'estimation';
    lockings(4).offset  = lockeddata.trialinfo(:, 15) - lockeddata.trialinfo(:,1);
    lockings(4).trials  = find(~isnan(lockeddata.trialinfo(:,15)));
    lockings(4).prestim = 0.5;
    lockings(4).poststim = 2;
    
    lockings(5).name    = 'feedback';
    lockings(5).offset  = lockeddata.trialinfo(:, 4) - lockeddata.trialinfo(:,1) + round(lockeddata.fsample*2);
    lockings(5).trials  = find(lockeddata.trialinfo(:,5) == -1);
    lockings(5).prestim = 0.5;
    lockings(5).poststim = 2;
    
    lockings(6).name    = 'choice1';
    lockings(6).offset  = lockeddata.trialinfo(:,9) - lockeddata.trialinfo(:,1);
    lockings(6).trials  = find(lockeddata.trialinfo(:,5) == 1);
    lockings(6).prestim = 0.5;
    lockings(6).poststim = 0.75;
    
    lockings(7).name    = 'choice0';
    lockings(7).offset  = lockeddata.trialinfo(:,9) - lockeddata.trialinfo(:,1);
    lockings(7).trials  = find(lockeddata.trialinfo(:,5) == 0);
    lockings(7).prestim = 0.5;
    lockings(7).poststim = 0.75;
    
    lockings(8).name    = 'choice-1';
    lockings(8).offset  = lockeddata.trialinfo(:,9) - lockeddata.trialinfo(:,1);
    lockings(8).trials  = find(lockeddata.trialinfo(:,5) == -1);
    lockings(8).prestim = 0.5;
    lockings(8).poststim = 0.75;
    
    cols = linspecer(length(lockings), 'qualitative');
    
    % skip blinks
    for l = 2:length(lockings),
        
        switch lockings(l).name
            case 'blink'
                % shiftlock manually
                cfg         = [];
                cfg.trl     = [data.blinksmp(:,2) data.blinksmp(:,2)+5*data.fsample zeros(length(data.blinksmp))];
                locked      = ft_redefinetrial(cfg, data);
                locked      = ft_timelockanalysis([], locked);
                cfg         = [];
                cfg.channel = pupilchan;
                locked      = ft_selectdata(cfg, locked);
            otherwise
                % no baseline correction
                locked = shiftoffset_timelock(lockeddata, lockings(l).trials, ...
                    lockings(l).offset, lockings(l).prestim, lockings(l).poststim, data.fsample, 0);
        end
     
        grandavg_raw.lockings(l).time(find(sj==subjects), :) = locked.time;
        grandavg_raw.lockings(l).avg(find(sj==subjects), :)  = locked.avg;
        grandavg_raw.lockings(l).var(find(sj==subjects), :)  = locked.var;
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DECONVOLVED IRFs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    clearvars -except lockings sj subjects locked cols grandavg_raw grandavg_irf grandavg_resid grandavg_nuisances
    load(sprintf('~/Data/Commitment/CleanPupil/P%02d_alleye_residuals.mat', sj));
    pupilchan           = find(strcmp(data.label, 'EyePupil')==1);
    grandavg_nuisances(find(sj==subjects), :, :) = IRF.predic;

    if 0,
        
        for l = 1:size(IRF.mean, 2),
            
            if l <= 6, subplot(length(lockings),length(lockings),l + length(lockings)); else hold on; end
            plot(0:1/data.fsample:5, IRF.mean(:, l), 'color', cols(l, :));
            % title(IRF.name{l});
            set(gca, 'XLim', [0 5], 'box', 'off', 'TickDir', 'out', 'YLim', [-2 2]);
            if l == 1,   ylabel('IRF');  end
            
            grandavg_irf.lockings(l).time(find(sj==subjects), :) = 0:1/data.fsample:5;
            grandavg_irf.lockings(l).avg(find(sj==subjects), :)  = IRF.mean(:, l);
            
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % RESIDUALS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % epoch
    cfg         = [];
    cfg.trl     = data.trialinfo;
    lockeddata  = ft_redefinetrial(cfg, data);
    
    % lockings structure is still there, same as for raw data
    for l = 2:length(lockings),
        
        switch lockings(l).name
            case 'blink'
                % shiftlock myself
                cfg         = [];
                cfg.trl     = [data.blinksmp(:,2) data.blinksmp(:,2)+5*data.fsample zeros(length(data.blinksmp))];
                locked      = ft_redefinetrial(cfg, data);
                locked      = ft_timelockanalysis([], locked);
                cfg         = [];
                cfg.channel = pupilchan;
                locked      = ft_selectdata(cfg, locked);
            otherwise
                % no baseline correction
                locked = shiftoffset_timelock(lockeddata, lockings(l).trials, ...
                    lockings(l).offset, lockings(l).prestim, lockings(l).poststim, data.fsample, 0);
        end

        grandavg_resid.lockings(l).time(find(sj==subjects), :) = locked.time;
        grandavg_resid.lockings(l).avg(find(sj==subjects), :)  = locked.avg;
        grandavg_resid.lockings(l).var(find(sj==subjects), :)  = locked.var;
    end
    
end

savefast('~/Data/Commitment/Pupil/GA_allpupil_example.mat', 'grandavg_raw', 'grandavg_resid', 'lockings', 'grandavg_nuisances');
