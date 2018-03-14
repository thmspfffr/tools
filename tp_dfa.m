function dfa = tp_dfa(x,win,Fs,overlap,binnum)

% tp_dfa(x,win,Fs,overlap,binnum)
% Computes scaling exponent of time series using 
% detrended fluctuation analysis (Hardstone et al., 2012).
% Uses the following inputs:
% x:        signal to be analyzed
% win:      length of fitting window in seconds (default: [1 50])
% Fs:       sampling rate
% overlap:  overlap of windows (default: 0.5)
% binnum:   number of time bins for fitting (default: 10)
%

if isempty(win)
  win = [1 50];
end
if isempty(overlap)
  overlap = 0.5;
end
if isempty(binnum)
  binnum = 10;
end

d1      = log10(win(1)*Fs);
d2      = log10(win(2)*Fs);

wins    = unique(round(logspace(d1,d2,binnum)));

dfa.binnum  = length(wins);

for ichan = 1 : size(x,2) 
  
    fprintf('Processing channel %d ...\n',ichan);
    
    fluct = nan(size(wins,2),1);
    
    y = x(:,ichan)./mean(x(:,ichan));
    y = y-mean(y);
    y = cumsum(y);         		
        
    for i = 1 : size(wins,2)	
      
      lsq = zeros(floor(size(y,1)/(wins(i)*(1-overlap))),1);		
      cnt = 0;
      
      for nn = 1:round(wins(i)*(1-overlap)):size(y,1)-wins(i)
        cnt=cnt+1;
        % RMSE
        lsq(cnt) = (mean(fastdetrend(y(nn:nn+wins(i))).^2,1))^(1/2);		
      end
      
      fluct(i) = mean(lsq(1:cnt),1);						%
      
    end
    
   	dfa.y{ichan,1} = fluct;    
end

for ichan = 1: size(x,2) 

    X = [ones(1,length(wins))' log10(wins)'];
    Y = log10(dfa.y{ichan,1});
    tmp = X\Y;   
    dfa.exp(ichan)= tmp(2);
        
end

dfa.win = wins./Fs;
end

function signal = fastdetrend(signal)

  n = size(signal,1);
  if n == 1
    % If a row, turn into column vector
   signal = signal(:);			
  end


  N = size(signal,1);
  a = [zeros(N,1) ones(N,1)];
  a(1:N) = (1:N)'/N;

  signal = signal - a*(a\signal);

  if(n==1)
      signal = signal.';
  end
end