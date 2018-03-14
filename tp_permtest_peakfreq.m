function peak_freq = tp_permtest_peakfreq(dat,para)

for iperm = 1 : para.nperm
  
  fprintf('Permutation #%d\n', iperm);
      
  idx1 = randi(2,[28,1]);
  idx2 = 3-idx1;
  
  for i = 1 : 28
    
    x(:,i,1) = dat(:,i,idx1(i));
    x(:,i,2) = dat(:,i,idx2(i));
    
  end
  
  peak_freq(iperm,1) = tp_peakfreq(nanmean(x(:,:,1),2),para);
  peak_freq(iperm,2) = tp_peakfreq(nanmean(x(:,:,2),2),para);
  
end
