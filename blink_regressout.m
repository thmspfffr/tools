function [newpupil] = blink_regressout(pupildata, fsample, blinksmp, saccsmp, plotme, addBackSlowDrift)
% method by Knapen, de Gee et al. (2016)
% http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0155574 
% estimate canonical responses to blinks and saccades, then take those out
% of the pupil timecourse
%
% Anne Urai, 2016

% initialize settings
if ~exist('plotme', 'var'); plotme = true; end % plot all this stuff
if ~exist('addBackSlowDrift', 'var'); addBackSlowDrift = true; 
disp('adding back slow drift at the end'); end % plot all this stuff

% get the stuff we need
dat.pupil       = pupildata;
dat.time        = [1: length(pupildata)] ./ fsample; % time axis

% make sure there are no NaNs left in the pupil timecourse
dat.pupil(isnan(dat.pupil)) = interp1(find(~isnan(dat.pupil)), ...
    dat.pupil(~isnan(dat.pupil)), find(isnan(dat.pupil)), ...
    'nearest', 'extrap');
      
% ====================================================== %
% STEP 1: BAND-PASS FILTER
% ====================================================== %

if plotme,
    clf;  subplot(611); plot(dat.time,dat.pupil);
    axis tight; box off; ylabel('Interp');
    set(gca, 'xtick', []);
end

% filter the pupil timecourse twice
% first, get rid of slow drift
[b,a]      = butter(2, 0.01 / fsample, 'high');
dat.hpfilt = filtfilt(b,a, dat.pupil); % filter with zero lag
% get residuals for later
dat.lowfreqresid = dat.pupil - dat.hpfilt;

% also get rid of fast instrument noise
[b,a] = butter(2, 10 / fsample, 'low');
dat.bpfilt = filtfilt(b,a, dat.hpfilt);

if plotme,
    subplot(612); plot(dat.time,dat.bpfilt);
    axis tight; box off; ylabel('Bandpass');
    set(gca, 'xtick', []);
end

% ====================================================== %
% STEP 2: DOWNSAMPLING
% ====================================================== %

% remove idx outside range
% blinksmp(sum((blinksmp > length(dat.pupil)), 2) > 0, :) = [];
% saccsmp(sum((saccsmp > length(dat.pupil)), 2) > 0, :) = [];

newFs = 10;
downsmp = resample(dat.bpfilt, newFs, fsample);

% also downsample the sample idx for blinks and saccades
newblinksmp = round(blinksmp * (newFs/fsample));
newsaccsmp  = round(saccsmp * (newFs/fsample));

% ====================================================== %
% STEP 3: DECONVOLUTION
% ====================================================== %

% dont continue if the subject didnt blink
if isempty(newblinksmp),
    newpupil = dat.bpfilt;
    return;
end

clear designM
colcnt = 1;
for r = 1:2, % two regressors
    
    % create a logical vector to speed up the analyses
    switch r
        case 1
            thissmp = newblinksmp(:, 2);
        case 2
            thissmp = newsaccsmp(:, 2);
    end
    
    % put samples in design matrix at the right spot
    imp = [-1 5];
    thissmp = thissmp + newFs * imp(1); % shift by the offset we're interested in
    thissmp(thissmp < 1) = []; % remove too early samples
    samplelogical = zeros(length(downsmp), 1);
    samplelogical(thissmp) = 1; % put 1s in the design matrix
    
    % shift the starting points so the deconvolution catches -500 ms
    for c = 1 : range(imp)*newFs,
        % for each col, put ones at the next sample values
        designM(:, colcnt)   = samplelogical;
        samplelogical   = [0; samplelogical(1:end-1)]; % shift
        colcnt = colcnt + 1;
    end
end

% deconvolve to get IRFs
deconvolvedPupil       = pinv(designM) * downsmp(:); % pinv more robust than inv
deconvolvedPupil       = reshape(deconvolvedPupil, range(imp)*newFs, 2);

% ====================================================== %
% STEP 4: FIT CANONICAL IRF GAMMA FUNCS
% ====================================================== %

% double Erlang gamma function from Hoeks and Levelt, Wierda
x = linspace(0, range(imp), numel(deconvolvedPupil(:, 1)))';

% curves should start at 0 before fitting
deconvolvedPupil(:, 1) = deconvolvedPupil(:, 1) - deconvolvedPupil(1,1);
deconvolvedPupil(:, 2) = deconvolvedPupil(:, 2) - deconvolvedPupil(1,2);

blinkIRF = doublegamma_fit(x, deconvolvedPupil(:, 1), 'blink');
saccIRF  = doublegamma_fit(x, deconvolvedPupil(:, 2), 'sacc');

% check if the fits look good
if plotme,
    subplot(6,3,7);
    % first, blink stuff
    plotx = linspace(imp(1), imp(2), numel(deconvolvedPupil(:, 1)));
    plot(plotx, deconvolvedPupil(:, 1), '.b', plotx, blinkIRF, 'r-');
    legend off; box off; axis tight; ylabel('Blink');
    
    subplot(6,3,8);
    % first, blink stuff
    plot(plotx, deconvolvedPupil(:, 2), '.b', plotx, saccIRF, 'r-');
    legend off;
    box off; axis tight; ylabel('Saccade');
end

% ====================================================== %
% STEP 5: UPSAMPLE AND MAKE REGRESSORS
% ====================================================== %

% upsample to the sample rate of the data
blinkIRFup = resample(blinkIRF, fsample, newFs);

% convolve with timepoints of events in original data
samplelogical = zeros(length(dat.pupil), 1);
offset = round(blinksmp(:, 2) + imp(1)*fsample);
offset(offset<1) = []; % remove those that we cant catch so early
samplelogical(offset)   = 1; % put 1s at these events

% convolve
reg1 = cconv(samplelogical, blinkIRFup);
reg1 = reg1(1:length(dat.pupil))';

% SAME FOR SACCADES
% upsample to the sample rate of the data
saccIRFup = resample(saccIRF, fsample, newFs);

% convolve with timepoints of events in original data
samplelogical = zeros(length(dat.pupil), 1);
offset = round(saccsmp(:, 2) + imp(1)*fsample);
offset(offset<1) = []; % remove those that we cant catch so early
samplelogical(offset)   = 1; % put 1s at these events

% convolve
reg2 = cconv(samplelogical, saccIRFup);
reg2 = reg2(1:length(dat.pupil))';

% ====================================================== %
% STEP 6: REGRESS OUT THOSE RESPONSES FROM DATA
% ====================================================== %

% make design matrix
designM = [ones(size(reg1))' reg1' reg2'];

% estimate glm weights
% [b, ~, resid] = regress(dat.bpfilt', designM);

b = designM \ dat.bpfilt;
prediction = designM * b;
resid = dat.bpfilt - prediction;

% use residuals
dat.residualpupil = resid;

if plotme,
    subplot(614); plot(dat.time,prediction);
    axis tight; box off; ylabel('Prediction');
    set(gca, 'xtick', []);
    
    subplot(615);
    plot(dat.time, dat.bpfilt', 'color', [0.5 0.5 0.5]); hold on;
    plot(dat.time,dat.residualpupil);
    axis tight; box off; ylabel('Residual');
    set(gca, 'xtick', []);
end

% ====================================================== %
% STEP 7: ADD BACK THE SLOW DRIFT
% ====================================================== %

if addBackSlowDrift,
    newpupil = dat.lowfreqresid + dat.residualpupil;
else
    newpupil = dat.residualpupil;
end

if plotme,
    subplot(616);
    plot(dat.time, dat.pupil', 'color', [0.5 0.5 0.5]); hold on;
    plot(dat.time,newpupil);
    axis tight; box off; ylabel('+ lowfreq');
end

end

function irf = doublegamma_fit(x, y, type)
% return the sum of squared error of the fit

switch type
    case 'blink'
        startpt = [-1, 1, 10, 10, 1, 2.5];
        lb = [-Inf, 1e-25, 9, 8, 0.5, 1.5];
        ub = [ -1e-25, Inf, 11, 12, 2, 5];
        disp('fitting blink double gamma');
    case 'sacc'
        startpt = [-1, 1, 10, 10, 1, 2.5];
        lb = [-Inf, 1e-25, 9, 8, 0.5, 1.5];
        ub = [ -1e-25, Inf, 11, 12, 3, 5];
        disp('fitting saccade double gamma');
end

doublegamma = @(s1, s2, n1, n2, tmax1, tmax2, x) ...
    s1 * (x.^n1) .* exp((-n1.*x) ./ tmax1) + s2 * (x.^n2) .* exp((-n2.*x) ./ tmax2);
doublegamma = fittype(doublegamma);
fitobj = fit(x, y, doublegamma, ...
    'startpoint', startpt, 'lower', lb, 'upper', ub);
irf = feval(fitobj, x);

end