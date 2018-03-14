function dfa = tp_dfa(x,win,Fs,overlap)

% res_logbin = 15;

d1      = log10(win(1)*Fs);
d2      = log10(win(2)*Fs);

wins    = round(logspace(d1,d2,floor((d2-d1)*10)));	

for ichan = 1 : size(x,2) 
  
    fprintf('Processing channel %d ...\n',ichan);
    
    fluct = nan(size(wins,2),1);
    
    y = x(:,ichan)./mean(x(:,ichan));
    y = y-mean(y);
    y = cumsum(y);         		
        
    for i = 1:size(wins,2);				
      
      lsq = zeros(floor(size(y,1)/(wins(i)*(1-overlap))),1);		
      tt = 0;
      
      for nn = 1:round(wins(i)*(1-overlap)):size(y,1)-wins(i);	% we are going to detrend all windows in steps of 'n*(1-DFA_Overlap)'.
        tt=tt+1;
        lsq(tt) = (mean(detrend(y(nn:nn+wins(i))).^2,1))^(1/2);		% the square of the fluctuation around the local trend (of the window).
      end
      
      fluct(i) = mean(lsq(1:tt),1);						% the root-mean-square fluctuation of the integrated and detrended time series
    	fluct_tmp{i} = lsq(1:tt);
      
    end
    
   	dfa.y{ichan,1} = fluct;
   	dfa.y_all{ichan,1} = fluct_tmp{i};
    
end



for ichan = 1: size(x,2) 
%     dfa.plf(ichan) = plfit(dfa.y{ichan,1})/10;

    X = [ones(1,length(wins))' log10(wins)'];
    Y = log10(dfa.y{ichan,1});
    tmp = X\Y;   
    dfa.exp(ichan)= tmp(2);
        
end

dfa.win = wins;