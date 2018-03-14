% ChangeLog - see version control log at NBT website for details.
%
% Copyright (C) 2001  Klaus Linkenkaer-Hansen  (Neuronal Oscillations and Cognition group,
% Department of Integrative Neurophysiology, Center for Neurogenomics and Cognitive Research,
% Neuroscience Campus Amsterdam, VU University Amsterdam)
%
% Part of the Neurophysiological Biomarker Toolbox (NBT)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
% See Readme.txt for additional copyright information.
% -


function [DFACobject] = DetrendedFluctuationAmplitudeCorrelation(Signal, Fs, FitInterval, CalcInterval, DFA_Overlap, noBins, DFA_Plot, ChannelToPlot, res_logbin)

%noBins, number of bins to allow equal comparison across subjects for
%fluctuation scaling analysis

DFACobject = nbt_DFAC(size(Signal,2));

%% Get or set default parameters...
if (isempty(res_logbin))
    res_logbin = DFACobject.res_logbin;	% number of bins pr decade, i.e., the spacing of the logarithmic scale.
else
    DFACobject.res_logbin = res_logbin;
end

% force fixed values for FitInterval and CalcInterval
% FitInterval=[1 20];
% CalcInterval=[0.1 110];

% get parameters from Signalobject
% set parameters in DFAobject
DFACobject.FitInterval = FitInterval;
DFACobject.CalcInterval = CalcInterval;
DFACobject.Overlap = DFA_Overlap;
DFACobject.Fs = Fs;

DFACobject.noBins = noBins;


%******************************************************************************************************************
%% Begin analysis
% Find DFA_x

% Defining window sizes to be log-linearly increasing.
d1 = floor(log10(CalcInterval(1)*Fs));
d2 = ceil(log10(CalcInterval(2)*Fs));
DFA_x_t = round(logspace(d1,d2,(d2-d1)*res_logbin));	% creates vector from 10^d1 to 10^d2 with N log-equidistant points.
DFA_x = DFA_x_t((CalcInterval(1)*Fs <= DFA_x_t & DFA_x_t <= CalcInterval(2)*Fs));	% Only include log-bins in the time range of interest!
DFACobject.DFA_x = DFA_x;


%% Do check of FitInterval
if ((DFA_x(1)/Fs)>FitInterval(1) || (DFA_x(end)/Fs)<FitInterval(2))
    disp([ DFA_x(1) DFA_x(end)]/Fs)
    error('Scaling_DFA:WrongCalcInterval','The CalcInterval is smaller than the FitInterval')
end

%% The DFAC algorithm...

for ChannelID = 1:(size(Signal,2)) % loop over channels
    disp(ChannelID);
    
    if (isempty(DFACobject.DFA_y{ChannelID,1}))
        DFA_y = nan(size(DFA_x,2),1);
        DFAC_AmpR = nan(size(DFA_x,2),1); %Correlation of detrended normalized fluctuations and amp
        DFAC_AmpP = nan(size(DFA_x,2),1); %Correlation of detrended normalized fluctuations and amp
        DFAC_AmpR_Binned = nan(size(DFA_x,2),1); %Correlation of detrended normalized fluctuations and slope (Binned into noBins on sorted Amp)
        DFAC_AmpP_Binned = nan(size(DFA_x,2),1); %Correlation of detrended normalized fluctuations and slope (Binned into noBins on sorted Amp)
        try
        y = Signal(:,ChannelID) ;
        othery = Signal(:,ChannelID);       %original signal
        y = y-mean(y);                      %demeaned signal
        y = cumsum(y);                      % Integrate the above fluctuation time series ('y').
        
        for i = 1:size(DFA_x,2);			%'DFA_x(i)' is the window size, which increases uniformly on a log10 scale!
            D = zeros(floor(size(y,1)/(DFA_x(i)*(1-DFA_Overlap))),1);		% initialize vector for temporarily storing the root-mean-square of each detrended window.
            DAmp = zeros(floor(size(y,1)/(DFA_x(i)*(1-DFA_Overlap))),1);		% initialize vector for temporarily storing the average amplitude for each window
            
            tt = 0;
            for nn = 1:round(DFA_x(i)*(1-DFA_Overlap)):size(y,1)-DFA_x(i);	% we are going to detrend all windows in steps of 'n*(1-DFA_Overlap)'.
                tt=tt+1;
                DAmp(tt) = mean(othery(nn:nn+DFA_x(i)));    %mean of original signal
                temp = fastdetrend(y(nn:nn+DFA_x(i)) / DAmp(tt)); %To take account of fluctuations being bigger for larger amplitudes, divide each window by the original amplitude before detrending    
                D(tt) = (mean(temp.^2,1))^(1/2);		% the square of the fluctuation around the local trend (of the window).
            end
            DFA_y(i) = mean(D(1:tt),1);						% the root-mean-square fluctuation of the integrated and detrended time series
            [DFAC_AmpR(i),DFAC_AmpP(i)]  = corr(D(1:tt),DAmp(1:tt));    % the Detrended fluctuation amplitude correlation
            
            
            myFit = LinearModel.fit(zscore(DAmp(1:tt)),zscore(D(1:tt)), 'y~x1:x1 + x1');
            DFAC_SlopeR(i,1) = myFit.Coefficients.Estimate(3);
            DFAC_SlopeP(i,1) = myFit.Coefficients.pValue(3);
            DFAC_SlopeR(i,2) = myFit.Coefficients.Estimate(2);
            DFAC_SlopeP(i,2) = myFit.Coefficients.pValue(2);
            
            %To allow comparisons between different window sizes and
            %signals of different lengths, bin the windows based on
            %amplitude
            try
            binWindows = floor(tt/noBins);
            [~,sortedAmps] = sort(DAmp(1:tt));
            fluctBin = zeros(noBins,1); % average fluctuation for Bin
            ampBin= zeros(noBins,1);  %Average Amp for Bin
            
            for kk = 1:noBins;
                binAmp = sortedAmps(1 + (kk-1)*binWindows:kk*binWindows);
                fluctBin(kk) = mean(D(binAmp));
                ampBin(kk) = mean(DAmp(binAmp));
            end
            [DFAC_AmpR_Binned(i),DFAC_AmpP_Binned(i)]  = corr(fluctBin(:),ampBin(:));
            
            
            myFit = LinearModel.fit(zscore(ampBin(:)),zscore(fluctBin(:)), 'y~x1:x1 +x1');
            DFAC_SlopeR_Binned(i,1) = myFit.Coefficients.Estimate(3);
            DFAC_SlopeP_Binned(i,1) = myFit.Coefficients.pValue(3);
            
            DFAC_SlopeR_Binned(i,2) = myFit.Coefficients.Estimate(2);
            DFAC_SlopeP_Binned(i,2) = myFit.Coefficients.pValue(2);
            
            catch me
                DFAC_SlopeR_Binned(i,1:2) = nan;
                DFAC_SlopeP_Binned(i,1:2) = nan;
            end
        end  					  	       			% -- the F(n) in eq. (1) in Peng et al. 1995.
        
        DFACobject.DFA_y{ChannelID,1} = DFA_y;
        DFACobject.DFAC_AmpR{ChannelID,1} = DFAC_AmpR;
        DFACobject.DFAC_AmpP{ChannelID,1} = DFAC_AmpP;
        DFACobject.DFAC_AmpR_Binned{ChannelID,1} = DFAC_AmpR_Binned;
        DFACobject.DFAC_AmpP_Binned{ChannelID,1} = DFAC_AmpP_Binned;
        DFACobject.DFAC_SlopeR_Binned{ChannelID,1} = DFAC_SlopeR_Binned(:,1);
        DFACobject.DFAC_SlopeR{ChannelID,1} = DFAC_SlopeR(:,1);
        DFACobject.DFAC_SlopeP_Binned{ChannelID,1} = DFAC_SlopeP_Binned(:,1);
        DFACobject.DFAC_SlopeP{ChannelID,1} = DFAC_SlopeP(:,1);
        
        DFACobject.DFAC_SlopeR_Binned{ChannelID,2} = DFAC_SlopeR_Binned(:,2);
        DFACobject.DFAC_SlopeR{ChannelID,2} = DFAC_SlopeR(:,2);
        DFACobject.DFAC_SlopeP_Binned{ChannelID,2} = DFAC_SlopeP_Binned(:,2);
        DFACobject.DFAC_SlopeP{ChannelID,2} = DFAC_SlopeP(:,2);
        catch me
        end
        
        
    end
end


%% Fitting power-law
% for ChannelID = 1:(size(Signal,2)) % loop over channels
%     try
%     DFA_y = DFACobject.DFA_y{ChannelID,1};
%     
%     DFA_SmallTime_LogSample = min(find(DFA_x>=CalcInterval(1)*Fs));		%
%     DFA_LargeTime_LogSample = max(find(DFA_x<=CalcInterval(2)*Fs));
%     DFA_SmallTimeFit_LogSample = min(find(DFA_x>=FitInterval(1)*Fs));
%     DFA_LargeTimeFit_LogSample = max(find(DFA_x<=FitInterval(2)*Fs));
%     X = [ones(1,DFA_LargeTimeFit_LogSample-DFA_SmallTimeFit_LogSample+1)' log10(DFA_x(DFA_SmallTimeFit_LogSample:DFA_LargeTimeFit_LogSample))'];
%     Y = log10(DFA_y(DFA_SmallTimeFit_LogSample:DFA_LargeTimeFit_LogSample));
%     DFA_exp = X\Y; %least-square fit
%     
%     DFA_exp = DFA_exp(2);
%     
%     DFACobject.markerValues(ChannelID,1) = DFA_exp;
%     catch me
%         DFACobject.markerValues(ChannelID,1) = nan;
%     end
% end




%% Plotting
if (DFA_Plot ~=0)
    if ~ishandle(DFA_Plot)		%see if any figure handle is set
        figure(DFA_Plot)
        DFA_Plot = axes;
    end
    ChannelID = ChannelToPlot;
    DFA_y = DFACobject.DFA_y{ChannelID,1};
    disp('Plotting Channel')
    disp(ChannelID)
    try
        axes(DFA_Plot)
    catch
        figure(DFA_Plot)
        axes(gca)
    end
    hold on
    plot(log10(DFA_x(DFA_SmallTimeFit_LogSample:DFA_LargeTimeFit_LogSample)/Fs),log10(DFA_y(DFA_SmallTimeFit_LogSample:DFA_LargeTimeFit_LogSample)),'ro')
    delete(findobj(DFA_Plot,'Type','Line','-not','Marker','o')) % delete any redundant lines
    LineHandle=lsline;
    try % delete any fits to the black points if they exist
        BlackHandle=findobj(DFA_Plot,'Color','k');
        for i=1:length(BlackHandle)
            delete(LineHandle(LineHandle == BlackHandle(i)))
        end
    catch
    end
    plot(log10(DFA_x(DFA_SmallTime_LogSample:DFA_LargeTime_LogSample)/Fs),log10(DFA_y(DFA_SmallTime_LogSample:DFA_LargeTime_LogSample)),'k.')
    grid on
    %     zoom on
    axis([log10(min(DFA_x/Fs))-0.1 log10(max(DFA_x/Fs))+0.1 log10(min(DFA_y(3:end)))-0.1 log10(max(DFA_y))+0.1])
    xlabel('log_{10}(time), [Seconds]','Fontsize',12)
    ylabel('log_{10} F(time)','Fontsize',12)
    title(['DFA-exp=', num2str(DFACobject.markerValues(ChannelID,1))],'Fontsize',12)
end


end

function [signal,sn] = fastdetrend(signal)
% fast detrending of "signal"
n = size(signal,1);
if n == 1,
    signal = signal(:); % make signal a row vector
end

% set up linear fitting
N = size(signal,1);
a = [zeros(N,1) ones(N,1)];
a(1:N) = (1:N)'/N;
sn = a\signal;
signal = signal - a*(sn); % remove best fit
sn = sn(1);
if(n==1)
    signal = signal.'; % return correct dimensions
end
end
