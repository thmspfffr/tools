% for each SJ, plot an overview of raw timelock, deconv and resid

clear all; close all; clc;
addpath(genpath('~/code/Commitment/'));
cd('~/Data/Commitment');
warning off % to avoid fieldtrip's redefinetrial filling the command line

subjects = [0:15];
figure; set(gcf, 'DefaultAxesFontSize', 4);

if 0,
    for sj = unique(subjects),
        
        clf;
        clearvars -except sj subjects grandavg_raw grandavg_irf grandavg_resid;
        clc; disp(sj);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% RAW
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        load(sprintf('~/Data/Commitment/Pupil/P%02d/P%02d_alleye_filtered.mat', sj, sj));
        pupilchan           = find(strcmp(data.label, 'EyePupil')==1);
        
        for t = 1:length(data.trial), assert(~any(any(isnan(data.trial{t}))));   end
        
        % epoch
        cfg         = [];
        cfg.trl     = data.trialinfo;
        lockeddata  = ft_redefinetrial(cfg, data);
        
        % make sure there are no trials with nans
        cfg = [];
        cfg.trials = true(size(lockeddata.trial));
        for t = 1:length(lockeddata.trial),
            if any(isnan(lockeddata.trial{t}(pupilchan, :))),
                cfg.trials(t) = false;
            end
        end
        lockeddata = ft_selectdata(cfg, lockeddata);
        
        % plot the shiftlock for each event
        lockings(1).name    = 'blink';
        lockings(1).offset  = data.blinksmp(:,2);
        lockings(1).trials  = [];
        lockings(2).name    = 'int1';
        lockings(2).offset  = lockeddata.trialinfo(:,4) - lockeddata.trialinfo(:,1);
        lockings(2).trials  = [];
        lockings(3).name    = 'int2';
        lockings(3).offset  = lockeddata.trialinfo(:, 11) - lockeddata.trialinfo(:,1);
        lockings(3).trials  = find(~isnan(lockeddata.trialinfo(:,11)));
        lockings(4).name    = 'estimation';
        lockings(4).offset  = lockeddata.trialinfo(:, 15) - lockeddata.trialinfo(:,1);
        lockings(4).trials  = find(~isnan(lockeddata.trialinfo(:,15)));
        lockings(5).name    = 'feedback';
        lockings(5).offset  = lockeddata.trialinfo(:, 4) - lockeddata.trialinfo(:,1) + round(lockeddata.fsample*2.75);
        lockings(5).trials  = find(lockeddata.trialinfo(:,5) == -1);
        lockings(6).name    = 'choice1';
        lockings(6).offset  = lockeddata.trialinfo(:,9) - lockeddata.trialinfo(:,1);
        lockings(6).trials  = find(lockeddata.trialinfo(:,5) == 1);
        lockings(7).name    = 'choice0';
        lockings(7).offset  = lockeddata.trialinfo(:,9) - lockeddata.trialinfo(:,1);
        lockings(7).trials  = find(lockeddata.trialinfo(:,5) == 0);
        lockings(8).name    = 'choice-1';
        lockings(8).offset  = lockeddata.trialinfo(:,9) - lockeddata.trialinfo(:,1);
        lockings(8).trials  = find(lockeddata.trialinfo(:,5) == -1);
        
        cols = linspecer(length(lockings), 'qualitative');
        
        for l = 1:length(lockings),
            disp(lockings(l).name);
            
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
                        lockings(l).offset, 0, 3.5, data.fsample, 0);
            end
            
            % plot
            if l <= 6, subplot(length(lockings),length(lockings),l); else hold on; end
            
            bh(l) = boundedline(locked.time, locked.avg, locked.var, 'cmap', cols(l, :));
            title(lockings(l).name);
            set(gca, 'XLim', [0 5], 'box', 'off', 'TickDir', 'out');
            if l == 1,   ylabel('raw');  end
            if l == length(lockings),
                title('choice');
                lh = legend([bh(end-2:end)], {lockings(end-2:end).name});
                sh = subplot(length(lockings),length(lockings),l); spos = get(sh, 'Position');
                axis off;
                set(lh, 'Position', spos, 'box', 'off');
            end
            
            grandavg_raw.lockings(l).time(find(sj==subjects), :) = locked.time;
            grandavg_raw.lockings(l).avg(find(sj==subjects), :)  = locked.avg;
            grandavg_raw.lockings(l).var(find(sj==subjects), :)  = locked.var;
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % DECONVOLVED IRFs
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        clearvars -except lockings sj subjects locked cols grandavg_raw grandavg_irf grandavg_resid
        load(sprintf('~/Data/Commitment/Pupil/P%02d/P%02d_alleye_residuals.mat', sj, sj));
        pupilchan           = find(strcmp(data.label, 'EyePupil')==1);
        
        for l = 1:size(IRF.mean, 2),
            
            if l <= 6, subplot(length(lockings),length(lockings),l + length(lockings)); else hold on; end
            plot(0:1/data.fsample:6, IRF.mean(:, l), 'color', cols(l, :));
            % title(IRF.name{l});
            set(gca, 'XLim', [0 5], 'box', 'off', 'TickDir', 'out');
            if l == 1,   ylabel('IRF');  end
            
            grandavg_irf.lockings(l).time(find(sj==subjects), :) = 0:1/data.fsample:6;
            grandavg_irf.lockings(l).avg(find(sj==subjects), :)  = IRF.mean(:, l);
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % RESIDUALS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        clearvars -except lockings sj subjects locked cols grandavg_raw grandavg_irf grandavg_resid
        load(sprintf('~/Data/Commitment/Pupil/P%02d/P%02d_alleye_residuals.mat', sj, sj));
        pupilchan           = find(strcmp(data.label, 'EyePupil')==1);
        
        % epoch
        cfg         = [];
        cfg.trl     = data.trialinfo;
        lockeddata  = ft_redefinetrial(cfg, data);
        
        cfg = [];
        cfg.trials = true(size(lockeddata.trial));
        % make sure there are no trials with nans
        for t = 1:length(lockeddata.trial),
            if any(isnan(lockeddata.trial{t}(pupilchan, :))),
                cfg.trials(t) = false;
            end
        end
        lockeddata = ft_selectdata(cfg, lockeddata);
        
        % lockings structure is still there, same as for raw data
        for l = 1:length(lockings),
            
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
                        lockings(l).offset, 0, 3.5, data.fsample, 0);
            end
            
            if l <= 6, subplot(length(lockings),length(lockings),l + 2*length(lockings)); else hold on; end
            boundedline(locked.time, locked.avg, locked.var, 'cmap', cols(l, :));
            set(gca, 'XLim', [0 5], 'box', 'off', 'TickDir', 'out');
            % title(lockings(l).name);
            if l == 1,   ylabel('residuals');  end
            
            grandavg_resid.lockings(l).time(find(sj==subjects), :) = locked.time;
            grandavg_resid.lockings(l).avg(find(sj==subjects), :)  = locked.avg;
            grandavg_resid.lockings(l).var(find(sj==subjects), :)  = locked.var;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % LAYOUT
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        suplabel(sprintf('P%02d', sj), 't');
        print(gcf, '-depsc', '-painters', sprintf('~/Data/Commitment/Figures/P%02d_allpupil.eps', sj));
        
    end
    
    savefast('~/Data/Commitment/Pupil/GA_allpupil.mat', 'grandavg_raw', 'grandavg_irf', 'grandavg_resid', 'lockings');

end


if length(subjects) > 5,
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% GRAND AVERAGE PLOT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cd /Users/anneurai/Dropbox/code/Commitment;
    
    clear all;
    load('~/Data/Commitment/Pupil/GA_allpupil.mat');
    
    % get significant timewindow
    
    cols = linspecer(length(lockings), 'qualitative');
    cols(6, :) = cols(2, :);
    cols(7, :) = cols(1, :);
    
    figure; set(gcf, 'DefaultAxesFontSize', 5);
    
    GA(1) = grandavg_raw;
    GA(2) = grandavg_irf;
    GA(3) = grandavg_resid;
    
    for g = 1:3, 
        for l = 1:7,
            
            if l <= 6, subplot(length(lockings),length(lockings),l + (g-1)*length(lockings)); else hold on; end
            
            bh(l) = boundedline(squeeze(mean(GA(g).lockings(l).time)), ...
                squeeze(nanmean(GA(g).lockings(l).avg)), ...
                squeeze(nanstd(GA(g).lockings(l).avg)) ./ sqrt(size(GA(g).lockings(l).avg, 1)), ...
                'cmap', cols(l, :));
            
            set(gca, 'XLim', [0 4], 'box', 'off', 'TickDir', 'out');
            
            hold on;
            
            if l > 6,
            % add significant timewin
            [significanttimewin] = stats_choiceVnochoice(GA(g));
            axis tight;
            
            ylims = get(gca, 'YLim');
            lh = line([significanttimewin(1) significanttimewin(1)], ylims); set(lh, 'color', 'k', 'linestyle', ':');
            lh = line([significanttimewin(end) significanttimewin(end)], ylims); set(lh, 'color', 'k', 'linestyle', ':');
            end
            
            % labels
            if l == 1,
                switch g
                    case 1
                        ylabel({'Pupil size (z)'});
                    case 2
                        ylabel('IRF');
                    case 3
                        ylabel('Residuals');
                end
            end
            % titles
            if g == 1,
                if l < 7,
                    title(lockings(l).name);
                elseif l == 7,
                    title('choice');
                end
            end
        end
        
    end
    
    xlabel('Time from response (s)');
     
    %lh = legend([bh(end-2:end)], {lockings(end-2:end).name});
    lh = legend([bh(6:7)], {lockings(6:7).name});
    
    sh = subplot(length(lockings),length(lockings),l+ (g-1)*length(lockings));
    spos = get(sh, 'Position'); axis off;
    set(lh, 'Position', spos, 'box', 'off');
    
    suplabel(sprintf('GA, n = %d', size(GA(g).lockings(l).avg, 1)), 't');
    print(gcf, '-depsc', '-painters', sprintf('~/Dropbox/Writing/Commitment/Figures/GA_allpupil.eps'));
end