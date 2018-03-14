function y = tp_powspec(x,win,Fs,isplot)


if ~isreal(x)
  error('Input signal is complex valued.')
end

if isempty(win)
  segleng = length(x);
  xf = abs(fft(x));
  xf = xf(1:length(x)/2+1);

else
  
  segleng  = length(win);
  segshift = segleng / 2;
	nseg     = floor((length(x)-segleng)/segshift+1);

  xf = zeros(1,length(segleng));
  
  for iseg = 1 : nseg
  
    dloc = x((iseg-1)*segshift+1:(iseg-1)*segshift+segleng);
    dloc = dloc.*win;
    
    xf = xf + abs(fft(dloc));
   
  end
   
  xf = nanmean(xf,2);
  xf = xf(1:segleng/2+1);

end   

f= 0 : Fs/segleng:Fs/2;

if isplot
  plot(f(2:end),xf(2:end))
  xlabel('Frequency [Hz]');
  ylabel('Amplitude');
end


