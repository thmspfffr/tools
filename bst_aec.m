function [c, LowFreqs, HighFreqs] = bst_aec(F1,F2, sRate, bandNesting, bandNested, isUseParallel, isUseMex, numfreqs)
% BST_PAC: Calculate the DirectPAC metric for all the input time series.
%
% USAGE:  sPAC = bst_pac(F, sRate, bandNesting, bandNested, isUseParallel=0, isUseMex=0, numfreqs=[])
% 
% INPUTS:
%    - F             : Signal time series [nSignals x nTime]
%    - sRate         : Signal sampling rate (in Hz)
%    - bandNesting   : Candidate frequency band of phase driving oscillations e.g., [0.5 48] Hz
%                      Note that cycle minimal frequency bandNesting(1) needs to be
%                      at least 10 times smaller than signal length (duration)    
%    - bandNested    : Candidate frequency band of nested oscillatiosn e.g., [48,300] Hz
%    - isUseParallel : If 1, use parallel processing toolbox
%    - isUseMex      : If 1, use mex file instead of matlab loop
%    - numfreqs      : Number of frequency bins to use (if 0 or empty, use round(sRate/9))
% 
% OUTPUTS:
%    - DirectPAC   : Full array of direct PAC measures for all frequyency pairs
%    - LowFreqs    : List of nesting frequency for the DirectPAC maps (frequency for phase)
%    - NestedFreq  : List of nested frequency for the DirectPAC maps (frequency for amplitude)
%
% DOCUMENTATION:
%    - For more information, please refer to the method described in the following article:
%         ?zkurt TE, Schnitzler A, 
%         "A critical note on the definition of phase-amplitude cross-frequency coupling" 
%         J Neurosci Methods. 2011 Oct 15;201(2):438-43
%    - The current code is inspired from Ryan Canolty's code provided originally with the article:
%         Canolty RT, Edwards E, Dalal SS, Soltani M, Nagarajan SS, Kirsch HE, Berger MS, Barbaro NM, Knight RT,
%         "High gamma power is phase-locked to theta oscillations in human neocortex",
%         Science, 2006 Sep 15;313(5793):1626-8.

% @=============================================================================
% This function is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2018 University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Ryan Canolty, 2006
%          Esther Florin, Soheila Samiee, Sylvain Baillet, 2011-2013
%          Francois Tadel, 2013-2014

% Parse inputs
if (nargin < 7) || isempty(numfreqs) || isequal(numfreqs,0)
    numfreqs = [];
end
if (nargin < 6) || isempty(isUseMex)
    isUseMex = 0;
end
if (nargin < 5) || isempty(isUseParallel)
    isUseParallel = 0;
end
if (nargin < 4) || isempty(bandNesting) || isempty(bandNested)
    bandNesting = [];
    bandNested  = [];
end
% Input dimensions
[nSignals, nTime] = size(F1);
% If there is only one signal: disable parallel processing
if isUseParallel && (nSignals == 1)
    isUseParallel = 0;
end

% ===== DEFINE LOW/HIGH FREQ BINNING =====
% Initial code from Esther (deprecated)
if isempty(bandNesting)
    % === LOW FREQ ===
    centerLow = [2:.5:12, 14:2:48];
    %centerLow = [0.5:.2:2, 2.5:.5:12, 14:2:48];
    % === HIGH FREQ ===
    % Definitions
    fmin = 1;
    fmax = 250;
    numfreqs = 70;
    fstep = 0.75;
    % Calculate center frequencies
    temp1 = (0:numfreqs-1) * fstep;
    temp2 = logspace(log10(fmin), log10(fmax), numfreqs);
    temp2 = (temp2-temp2(1)) * ((temp2(end)-temp1(end)) / temp2(end)) + temp2(1);
    centerHigh = temp1 + temp2;
    % Taking only frequencies from 80 to 150 Hz
    centerHigh = centerHigh(51:62);
    % Group both
    chirpCenterFreqs = [centerLow, centerHigh];
    lfreq = 1:length(centerLow);
    hfreq = length(centerLow) + (1:length(centerHigh));
% New Esther/Soheila code
else
    % Definitions
    fmin = min(bandNesting);
    fmax = sRate/3;
    if isempty(numfreqs)
        numfreqs = round(sRate/9);
        fstep = 0.75;
    else
        fstep = 0.75 / numfreqs * fmax;
    end
    % Calculate center frequencies
    temp1 = (0:numfreqs-1) * fstep;
    temp2 = logspace(log10(fmin), log10(fmax), numfreqs);
    temp2 = (temp2-temp2(1)) * ((temp2(end)-temp1(end)) / temp2(end)) + temp2(1);
    chirpCenterFreqs = temp1 + temp2;
    % Remove unused frequencies
    chirpCenterFreqs(chirpCenterFreqs > max(bandNested)) = [];      %%% ESTHER
    chirpCenterFreqs((chirpCenterFreqs < min(bandNested)) & (chirpCenterFreqs >= max(bandNesting))) = [];      %%% ESTHER
    % Indices of center frequencies in the upper frequency range
    hfreq = find( chirpCenterFreqs >= min(bandNested) );
    % Number of cf bins to evaluate for PAC with lower-frequency oscillations
    % lfreq = find(chirpCenterFreqs < min(bandNested));
    lfreq = find(chirpCenterFreqs < max(bandNesting));   %%% ESTHER
end

% ===== CALCULATE CHIRPLETS =====
% Calculate chirplets
[chirpF, Freqs] = bst_chirplet(sRate, nTime, chirpCenterFreqs);

% ===== BAND PASS INPUT SIGNAL =====
% NOT IN ESTHER'S ORIGINAL CODE
if ~isempty(bandNesting)
    % Apply the band-pass filter to eliminate all the frequencies that we are not interested in
    bandLow  = min(bandNesting(1)*.5, sRate/2);
    bandHigh = min(bandNested(end)*1.5, sRate/2 - min(20, sRate/2 * 0.2) - 1);
    F1 = bst_bandpass_fft(F1, sRate, bandLow, bandHigh, 1, 1);
    F2 = bst_bandpass_fft(F2, sRate, bandLow, bandHigh, 1, 1);
end

% ===== FFT OF SIGNALS =====
% Transform sensor time series into analytic signals
F_fft1 = fft(F1, length(Freqs), 2);
% This step scales analytic signal such that: real(analytic_signal) = raw_signal
% but note that analytic signal energy is double that of raw signal energy
F_fft1(:,Freqs<0) = 0;
F_fft1(:,Freqs>0) = 2 * F_fft1(:,Freqs>0);
clear F1;

F_fft2 = fft(F2, length(Freqs), 2);
% This step scales analytic signal such that: real(analytic_signal) = raw_signal
% but note that analytic signal energy is double that of raw signal energy
F_fft2(:,Freqs<0) = 0;
F_fft2(:,Freqs>0) = 2 * F_fft2(:,Freqs>0);
clear F2;

% Define minimal frequency support with non-zeros chirplet coefficients
[row,scol] = find(F_fft1 ~= 0);
scol = max(scol)+1;
[chirprow,chirpcol] = find(squeeze(chirpF(1,:,:)) ~= 0);
chirprow = max(chirprow)+1;
% Minimal number of frequency coefficients
nfcomponents = min(chirprow,scol); 
clear row scol chirprow chirpcol


% ===== CALULATE AEC =====
% Filter signal in frequency domain
F_fft1 = bst_bsxfun(@times, F_fft1(:, 1:nfcomponents, ones(1,length(chirpCenterFreqs))), ...
                           chirpF(1,1:nfcomponents,:));
% Convert back to time domain
fs1 = ifft(F_fft1, length(Freqs), 2);
clear F_fft;

[row,scol] = find(F_fft2 ~= 0);
scol = max(scol)+1;
[chirprow,chirpcol] = find(squeeze(chirpF(1,:,:)) ~= 0);
chirprow = max(chirprow)+1;
% Minimal number of frequency coefficients
nfcomponents = min(chirprow,scol); 
clear row scol chirprow chirpcol

% Filter signal in frequency domain
F_fft2 = bst_bsxfun(@times, F_fft2(:, 1:nfcomponents, ones(1,length(chirpCenterFreqs))), ...
                           chirpF(1,1:nfcomponents,:));
% Convert back to time domain
fs2 = ifft(F_fft2, length(Freqs), 2);
clear F_fft;

% Magnitude and phase about each chirplet center frequency
AMP1 = abs( fs1(:, 1:nTime, lfreq ) );
AMP2 = abs( fs2(:, 1:nTime, hfreq ) );

for ilfreq = 1 : length(lfreq)
  for ilhfreq = 1 : length(hfreq)
%     c(ilfreq,ilhfreq) = ([AMP1(1,:,ilfreq)-mean(AMP1(1,:,ilfreq))]*[AMP2(1,:,ilhfreq)'-mean(AMP2(1,:,ilhfreq))])/ (std(AMP1(1,:,ilfreq))*std(AMP2(1,:,ilhfreq)));
    c(ilfreq,ilhfreq) = corr(AMP1(1,:,ilfreq)',AMP2(1,:,ilhfreq)');
  end
end

clear fs1 fs2;

LowFreqs  = chirpCenterFreqs(lfreq);
HighFreqs = chirpCenterFreqs(hfreq);
