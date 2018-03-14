function [sig,freqs] = tp_kuramoto(dat,para)

% Computes Kuramoto order parameter of time series in dat
% para.epleng
% para.segleng
% para.segshift
% para.fsample

[ndat,nchan]=size(dat);

mywindow = repmat(hanning(para.segleng),1,nchan);

if isfield(para,'mywindow');
  mywindow=repmat(para.mywindow,1,nchan);
end

if isfield(para,'zeropad') && para.zeropad>0
  mywindow=[mywindow;zeros(para.zeropad,nchan)];
end

if ~isfield(para,'segshift')
  para.segshift = para.segleng/2;
end

nave  = 0;

nep   = floor(ndat/para.epleng);
nseg  = floor((para.epleng-para.segleng)/para.segshift)+1; 

freqs = 0 : para.fsample/para.segleng : para.fsample/2; 

idx = find(freqs>=8 & freqs<=12);

for j = 1 : nep;
  
  dataep=dat((j-1)*para.epleng+1:j*para.epleng,:);
  
  for i = 1 : nseg %average over all segments;
    
    dataloc=dataep((i-1)*para.segshift+1:(i-1)*para.segshift+para.segleng,:);
    
    if isfield(para,'mydetrend') && para.mydetrend==1;
      dataloc=detrend(dataloc);
    end
    
    if isfield(para,'zeropad') && para.zeropad>0
      dataloc=[dataloc;zeros(zeropad,nchan)];
    end
    
    datalocfft = fft(dataloc.*mywindow);
    sig(:,i) = mean(datalocfft(idx,:));
    
  end
end

