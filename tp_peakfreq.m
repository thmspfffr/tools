function peak = tp_peakfreq(pow,para)
% tp_peakfreq.m
% Computes peak frequency based on different methods
% ----- para.method -----
% 'com' (center of mass)
% 'peak' (peak of spectrum)
% 'peakcom' (mixture of com and peak)
% ----- Other settings -----
% para.f:           Frequencies 
% para.win:         How many values around peak for peak estimation?
% para.detrend:     Detrend power spectrum?
% para.detrendfreq: Frequency range that should be used for detrending
% -----
% thms.pfffr@gmail.com, 03/2018 
% -----

if size(pow,1) < size(pow,2)
  pow = pow';
end
if size(para.f,1) > size(pow,2)
  para.f = para.f';
end

if ~isfield(para,'detrend')
  para.detrend = 0;
elseif isfield(para,'detrend') && para.detrend == 1
  detrend_idx = [find(para.f>=para.detrendfreq(1),1,'first') find(para.f<=para.detrendfreq(2),1,'last')];
  detrend_pow = log10(pow(detrend_idx(1):detrend_idx(2)));
  
  X = [ones(length(para.f(detrend_idx(1):detrend_idx(2))),1) log10(para.f(detrend_idx(1):detrend_idx(2)))'];
  Y = detrend_pow;
  regre = X\Y;
  
  regre(2)*log10(para.f(1:length(pow)))+regre(1)
  
  
  pow = detrend(pow);
end

switch para.method
  
  case 'com'    
    peak = para.f(para.f_select)*pow(para.f_select) / sum(pow(para.f_select));
  
  case 'peak'
    [~,i] = max(pow(para.f_select));
    peak = para.f(para.f_select(i));
  
  case 'peakcom'    
    [~,i] = max(pow(para.f_select));
    i = para.f_select(1)+i-1;
    peak = para.f(i-para.win:i+para.win)*pow(i-para.win:i+para.win) / sum(pow(i-para.win:i+para.win)); 
  
  case 'gaussian'    
    a=fit(para.f(para.f_select)',pow(para.f_select),'gauss1');
    fun = a.a1*exp(-((para.f(para.f_select)-a.b1)/a.c1).^2);
    [~,i] = max(fun);
    peak = para.f(para.f_select(i));

end

