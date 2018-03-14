
function powspec(dat,fs)

% COMPUTES POWER SPECTRUM 
% AVERAGE OVER ALL CHANNELS

if size(dat,1) == 1
  
else
  for ichan = 1 : size(dat,1)
    N = length(dat);
    xdft = fft(dat(ichan,:));
    xdft = xdft(1:N/2+1);
    psdx(ichan,:,:) = (1/(fs*N)).*abs(xdft).^2;
    psdx(ichan,2:end-1) = 2*psdx(ichan,2:end-1);
  end
  
  psdx = squeeze(squeeze(nanmean(psdx,1)));

  freq = linspace(0,fs/2,length(psdx)); 
  strt = find(freq > 4,1,'first');
  ends = find(freq > 60,1,'first');

  subplot(2,1,1);
  plot(log10(freq(strt:end)),log10(psdx(strt:end))); grid on;
  set(gca,'TickDir','out','XTick',[0:0.2:2.5],'XTickLabel',round([10.^(0:0.2:2.5)]))
  xlabel('Frequency (Hz)'); ylabel('(dB/Hz)');

  strt = find(freq > 15,1,'first');
  ends = find(freq > 18,1,'first');

  subplot(2,1,2);
  plot(log10(freq(strt:ends)),log10(psdx(strt:ends))); grid on;
  set(gca,'TickDir','out','XTick',[log10(15):0.01:log10(18)],'XTickLabel',[10.^(log10(15):0.01:log10(18))])
  xlabel('Frequency (Hz)'); ylabel('(dB/Hz)');
  linrng = [min(psdx(strt:ends)) max(psdx(strt:ends))];
  line([log10(16.7) log10(16.7)],[log10(linrng(1))-0.5 log10(linrng(2))+0.5],'color',[0 0 0])
  axis tight
end