% Simulate time series with long-range temporal correlations
% Add oscillation on top
% See how this affects DFA exponent and fluctuation function

ampl = 0.2;
foi_range = 0.05:0.05:0.5;
figure_w;
for ifoi = 1 : 9
  len = 100000;
  fs = 100;
  dt = 1/fs;
  t = dt:dt:len/fs;
  F = foi_range(ifoi);
  osc = sin(2*pi*F*t);


  
x = ffGn(len,0.9,0.5,0);

  x1 = x+ampl.*osc;

  dfa = tp_dfa(x',[1 100],100,0.5,25);
  dfa1 = tp_dfa(x1',[1 100],100,0.5,25);

  subplot(3,3,ifoi)
  
  dfa.y{1} = (log10(dfa.y{1}) - min(log10(dfa.y{1})))/(max(log10(dfa.y{1}))-min(log10(dfa.y{1})));
  dfa1.y{1} = (log10(dfa1.y{1}) - min(log10(dfa1.y{1})))/(max(log10(dfa1.y{1}))-min(log10(dfa1.y{1})));
  
  plot(log10(dfa.win), dfa.y{1},'-'); hold on
  plot(log10(dfa1.win), dfa1.y{1},'-')
  title(sprintf('Freq: %.3f Hz',foi_range(ifoi)))
  
  line([log10(1/foi_range(ifoi)) log10(1/foi_range(ifoi))],[0 1],'linestyle',':','color','k')
  set(gca,'xtick',log10(dfa1.win(1:4:end)),'xticklabel',round(dfa1.win(1:4:end)))
  
  axis([log10(dfa1.win(1)) log10(dfa1.win(end)) 0 1])
  tp_editplots
  
  if ifoi >6 
    xlabel('Window length [s]')
  end
  if ifoi == 1 || ifoi == 4 || ifoi == 7
    ylabel('Fluctuation (norm.)')
  end
  
end


