function y = tp_interp_blinks(sig,dist,thresh)


% -----------------------------------------------------------------
% INTERPOLATE BLINKS
% -----------------------------------------------------------------
x = abs(diff(zscore(sig)));
[~,idx]=findpeaks(double(x>thresh),'MinPeakDistance',dist);

for iidx = 1 : size(idx,1)
  if idx(iidx)-50>0 && idx(iidx)+600<length(sig)
    sig(idx(iidx)-50:idx(iidx)+600)=NaN;
  elseif idx(iidx)-50<0
    sig(2:idx(iidx)+600)=NaN;
  elseif idx(iidx)+600>length(sig)
    sig(idx(iidx)-50:end)=NaN;
  end
end

clear idx x

sig = fixgaps(sig,'pchip');

if isnan(sig(1))
  tt = find(~isnan(sig),1,'first');
  sig(1) = sig(tt);
elseif isnan(sig(end))
  tt = find(~isnan(sig),1,'last');
  sig(end) = sig(tt);
end
y = fixgaps(sig,'pchip');