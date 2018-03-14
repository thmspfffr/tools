function dcca = tp_dcca(data,para)
% Takes the following inputs:
% para.fs:      sampling rate of input signal
% para.overlap: window overlap for covariance estimation
% para.w_calc:  window limits over which dcca is computed
% para.w_fit:   window limits over which dcca exp is estimated

% This function implements detrended cross-correlation analysis (DCCA)
% as described in: Podobnik, B. & Stanley, H. E. (2008). Detrended
% Cross-Correlation Analysis: A New Method For Analyzing Two Nonstationary
% Time-Series. Phys. Rev. Lett. Review

addpath ~/pconn/matlab/

wins = floor(logspace(0,2,10)*para.fs);

[t nchan] = size(data);

% integrate data / signal profile
data = cumsum(data);

for iwin = 1 : length(wins)
  
  fprintf('Estimating fluctuation in window %d of %d ...\n',iwin,length(wins));
  
  winlen = wins(iwin);
  nwin   = floor(size(data,1)/winlen);
  
  covdat = zeros(nchan,nchan);
  
  for twin = 1 : nwin

    % cut out signal
    dataseg = data(floor((twin-1)*winlen*para.overlap)+1:floor((twin-1)*winlen*para.overlap)+winlen,:);
    % detrend signal
    dataseg = tp_fastdetrend(dataseg);
    % covariance
    covdat  = covdat + cov(dataseg);
    
  end
  
%   covdat = covdat ./ nwin;
  
  allcov(:,:,iwin) = abs(covdat);
  
end


codfa = nan(nchan,nchan);

for i = 1 : nchan
  fprintf('Estimating scaling for node %d of %d ...\n',i,nchan)
  for j = 1 : nchan
    if i == j
      continue
    end
    
      s = pconn_regress(log10(wins),log10(squeeze(allcov(i,j,:)))); 
      dcca(i,j) = s(2);
      
  end
end





    
    
    


