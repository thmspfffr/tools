% nbt_DFA(NumChannels) - Creates a DFA biomarker object for 'NumSubjects' number of
% subjects
%
% Usage:
%   >>  DFAobject = DFA( NumSubjects, NumChannels);
%
% Inputs:
%   NumChannels - Number of Channels
%
% Outputs:
%   DFAobject     - DFA Biomarker object
%
% Example:
%
% References:
%
%
% See also:
% NBT_SCALING_DFA

%------------------------------------------------------------------------------------
% Originally created by Simon-Shlomo Poil (2008), see NBT website (http://www.nbtwiki.net) for current email address
%------------------------------------------------------------------------------------
%
% ChangeLog - see version control log at NBT website for details.
%
% Copyright (C) 2008 Simon-Shlomo Poil  (Neuronal Oscillations and Cognition group,
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

classdef nbt_DFAC %< nbt_SignalBiomarker
    properties
        markerValues
        DFA_y
        DFA_x
        DFA_yPos
        DFA_yNeg
        
        Fs
        DFAC_SlopeR
        DFAC_SlopeP
        DFAC_AmpR
        DFAC_AmpP
        DFAC_SlopeR_Binned
        DFAC_SlopeP_Binned
        DFAC_AmpR_Binned
        DFAC_AmpP_Binned
        
        DFA_expPos;
        DFA_expNeg;
        
        FitInterval
        CalcInterval
        Overlap
        res_logbin
        noBins
    end
    %     properties (Constant)
    %         biomarkerType = {'nbt_SignalBiomarker'};
    %         units = {' '};
    %     end
    methods
        function DFACobject= nbt_DFAC(NumChannels)
            if nargin == 0
                NumChannels = 1;
            end
            %% Define DFA_x and DFA_y fields
            DFACobject.DFA_y  = cell(NumChannels, 1);
            DFACobject.DFA_yPos  = cell(NumChannels, 1);
            DFACobject.DFA_yNeg  = cell(NumChannels, 1);
            
            DFACobject.noBins = nan;
            
            
            DFACobject.DFAC_SlopeR  = cell(NumChannels, 1);
            DFACobject.DFAC_SlopeP  = cell(NumChannels, 1);
            DFACobject.DFAC_AmpR  = cell(NumChannels, 1);
            DFACobject.DFAC_AmpP  = cell(NumChannels, 1);
            DFACobject.DFAC_SlopeR_Binned  = cell(NumChannels, 1);
            DFACobject.DFAC_SlopeP_Binned  = cell(NumChannels, 1);
            DFACobject.DFAC_AmpR_Binned  = cell(NumChannels, 1);
            DFACobject.DFAC_AmpP_Binned  = cell(NumChannels, 1);
            
            DFACobject.DFA_expPos  = cell(NumChannels, 1);
            DFACobject.DFA_expNeg  = cell(NumChannels, 1);
            
            
            DFACobject.DFA_x  = [];
            
            %% Define DFA exponent field
            DFACobject.markerValues = nan(NumChannels, 1);
            %% Define fields for additional information
            %  DFACobject.primaryBiomarker = 'markerValues';
            %  DFACobject.biomarkers = {'markerValues'};
            DFACobject.FitInterval = nan(2,1);
            DFACobject.CalcInterval = nan(2,1);
            DFACobject.Overlap = NaN;
            % DFACobject.Condition = NaN;
            % DFACobject.lastUpdate = NaN;
            DFACobject.Fs = nan;
            DFACobject.res_logbin = 10;
            %     DFACobject = setUniqueIdentifiers(DFACobject);
        end
        %
        %         function LastUpdate=GetLog(DFAobject)
        %
        %             LastUpdate = DFAobject.lastUpdate;
        %             disp('The DFAobject was last updated')
        %             disp(LastUpdate)
        %         end
        %
        %         function BiomarkerObject = setUniqueIdentifiers(BiomarkerObject)
        %             BiomarkerObject.uniqueIdentifiers = {'frequencyRange','FitInterval','CalcInterval','Overlap'};
        %         end
        %
        %         function plot(DFAobject, ChannelID, DFA_Plot)
        %             if ~ishandle(DFA_Plot)		%see if any figure handle is set
        %                 figure(DFA_Plot)
        %                 DFA_Plot = axes;
        %             end
        %
        %             DFA_x = DFAobject.DFA_x;
        %             DFA_y = DFAobject.DFA_y(ChannelID,1);
        %             Fs = DFAobject.Fs;
        %             disp('Plotting Channel')
        %             disp(ChannelID)
        %             try
        %                 axes(DFA_Plot)
        %             catch
        %                 figure(DFA_Plot)
        %                 axes(gca)
        %             end
        %             hold on
        %             plot(log10(DFA_x(DFA_SmallTimeFit_LogSample:DFA_LargeTimeFit_LogSample)/Fs),log10(DFA_y(DFA_SmallTimeFit_LogSample:DFA_LargeTimeFit_LogSample)),'ro')
        %             delete(findobj('Type','Line','-not','Marker','o')) % delete any redundant lines
        %             LineHandle=lsline;
        %             try % delete any fits to the black points if the exist
        %                 BlackHandle=findobj('Color','k');
        %                 for i=1:length(BlackHandle)
        %                     delete(LineHandle(LineHandle == BlackHandle(i)))
        %                 end
        %             catch
        %             end
        %             plot(log10(DFA_x(DFA_SmallTime_LogSample:DFA_LargeTime_LogSample)/Fs),log10(DFA_y(DFA_SmallTime_LogSample:DFA_LargeTime_LogSample)),'k.')
        %             grid on
        %             zoom on
        %             axis([log10(min(DFA_x/Fs))-0.1 log10(max(DFA_x/Fs))+0.1 log10(min(DFA_y(3:end)))-0.1 log10(max(DFA_y))+0.1])
        %             xlabel('log_{10}(time), [Seconds]','Fontsize',12)
        %             ylabel('log_{10} F(time)','Fontsize',12)
        %             title(['DFA-exp=', num2str(DFAobject.DFAexp(ChannelID,SubjectID))],'Fontsize',12)
        %         end
        %
        
    end
end