function designM = makeDeconvDesignM(cfg, data)
% based on a datafile and config, return a design matrix to be used for
% deconvolution
% cfg.regressors = {'blink, 'int1', 'int2', 'feedback', 'estimation'}
% cfg.regressorscols = [Nan 4 6 9]
% cfg.impulse = in seconds;

% preallocate designM
designM = zeros(length(data.time{1}), cfg.impulse*data.fsample*length(cfg.regressors));
whos
disp('making design matrix');

for r = 1:length(cfg.regressors),
    disp(cfg.regressors{r});
    
    % get all samples of these events
    switch cfg.regressors{r}
        case 'blink'
            sample = data.blinksmp(:,2);
        case 'int1'
            sample = data.trialinfo(:,7);
        case 'choice'
            sample = data.trialinfo(:,12);
        case 'int2'
            sample = data.trialinfo(:,14);
        case 'estimation'
            sample = data.trialinfo(:,18);
        case 'feedback'
            sample = round(data.trialinfo(:, 7) + (2.75 * data.fsample));
            sample = sample(find(data.trialinfo(:,8) == -1));
        case 'choice1'
            sample = data.trialinfo((data.trialinfo(:,8)== 1) ,12);
        case 'choice0'
            sample = data.trialinfo((data.trialinfo(:,8)== 0) ,12);
        case 'choice-1'
            sample = data.trialinfo((data.trialinfo(:,8)== -1) ,12);
        otherwise
    end
    
    % remove bad samples
    sample(isnan(sample)) = [];
    
    % create a logical vector to speed up the analyses
    samplelogical           = zeros(length(data.time{1}), 1);
    samplelogical(sample)   = 1; % first sample of this regressor
    
    % put samples in design matrix at the right spot
    begincol = r + (r-1)*cfg.impulse*data.fsample;
    for c = begincol : begincol + cfg.impulse*data.fsample,
        % for each col, put ones at the next sample valuess
        designM(:, c)   = samplelogical;
        samplelogical   = [0; samplelogical(1:end-1)]; % shift
    end
end

end