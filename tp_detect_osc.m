function [osc1] = tp_detect_osc(x)

siz = size(x);
sig = zscore(x);

if var(x(size(x,1)/2:end))>0.001
  
  [~,lks]=findpeaks(sig(round(siz(1)/2):end),'MinPeakHeight',0,'MinPeakDistance',20);
  kk = std(diff(lks));

  if kk < 20 && length(lks)>20
    osc1 = 1;
  else 
    osc1 = 0;
  end
  
else
  osc1 = 0;
end

% osc2 = ;
% osc3 = mean(x(size(x,1)/2:end))>0;


