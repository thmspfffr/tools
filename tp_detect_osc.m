function osc = tp_detect_osc(x)
% Detects oscillations in simulated time series. See
% pmod_wc_wholebrain_detosc.m for details of the simulation.

sig = single(zscore(x));
siz = size(x,1);

seglen = 300;

nseg = floor(size(sig,1) / seglen);

for iseg = 1 : nseg
  
  tmp_max(iseg) = max(sig( (iseg-1)*seglen+1:iseg*seglen));
  tmp_min(iseg) = min(sig( (iseg-1)*seglen+1:iseg*seglen));
  
  dif(iseg) = tmp_max(iseg)-tmp_min(iseg);
  
end


if any(dif==0)
  osc = 0;
else
  if all(diff(dif)<0)
    osc = 0;
  else
    osc = 1;
  end
end
% end

% if any(dif==0)
%   osc = 0;
% else
%   if (sum(diff(dif)<0)/length(diff(dif)))>0.9
%     osc = 0;
%   else
%     osc = 1;
%   end
% end


% clust = bwlabel(idx);
% 
% for i = 1 : max(clust)
%   
%   ampl(i) = mean(sig(clust==i));
%   
% end

% [h,p]=ttest(ampl(60:79)', ampl(80:99)');
% [h,p]=ttest(ampl(61:90)', ampl(91:100)');

% osc = ~h;

% if var(sig(floor(size(sig,1)/2:end)))>0.1
%   
%   [~,lks]=findpeaks(sig(round(siz(1)/2):end),'MinPeakHeight',0,'MinPeakDistance',20);
%   kk = std(diff(lks));
% 
%   if kk < 5 && length(lks)>10
%     osc1 = 1;
%   else 
%     osc1 = 0;
%   end
%   
% else
%   osc1 = 0;
% end

% osc2 = ;
% osc3 = mean(x(size(x,1)/2:end))>0;
% 
% figure; plot(sig); hold on
% line([1000 1000],[min(sig) max(sig)],'color','k','linewidth',2,'linestyle','--')
% line([2000 2000],[min(sig) max(sig)],'color','k','linewidth',2,'linestyle','--')
% line([5000 5000],[min(sig) max(sig)],'color','r','linewidth',2,'linestyle','--')
% line([5200 5200],[min(sig) max(sig)],'color','r','linewidth',2,'linestyle','--')
% axis([1 length(sig) min(sig) max(sig)])
% 
% x = sig(900:2100);
% figure; plot(x); hold on
% line([100 100],[min(x) max(x)],'color','k','linewidth',2,'linestyle','--')
% line([1000 1000],[min(x) max(x)],'color','k','linewidth',2,'linestyle','--')
% axis([1 length(x) min(x)-0.1*(abs(min(x))) max(x)+0.1*(abs(max(x)))])
% 
% x = sig(4900:5300);
% figure; plot(x); hold on
% line([100 100],[min(x) max(x)],'color','r','linewidth',2,'linestyle','--')
% line([300 300],[min(x) max(x)],'color','r','linewidth',2,'linestyle','--')
% axis([1 length(x) min(x)-0.01*(abs(min(x))) max(x)+0.01*(abs(max(x)))])
