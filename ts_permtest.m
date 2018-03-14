function [h p] = permtest(x,y,para);
% 
% if ~exist(para)
%   para.alpha  = 0.05;
%   para.tail   = 0;
%   para.paried = 0;
%   para.nperm  = 10000;
% end

d_emp = mean(x)-mean(y);

if ~para.paired
  
  dat = [x(:,1);x(:,2)];

  for iperm = 1 : para.nperm

    idx = randperm(size(dat,1));

    x = dat(idx(1:size(dat,1)/2));
    y = dat(idx(size(dat,1)/2+1:end));
    
    d(iperm) = mean(x)-mean(y);
    
  end
    
elseif para.paired
  
  dat = [x y]; clear x y
 
  for iperm = 1 : para.nperm

    disp(sprintf('Perm #%d',iperm));
    
    idx1 = randi(2,[size(dat,1),1]);
    idx2 = 3-idx1;
                 
    for i = 1 : length(idx1)
      
      x(i) = dat(i,idx1(i));
      y(i) = dat(i,idx2(i));
      
    end
    
    d_perm(iperm) = mean(x)-mean(y); clear x y
    
  end
  
end

switch para.tail
  case 0
    p = 1-sum(abs(d_emp)>abs(d_perm))/para.nperm;
  case 1
    p = 1-sum(d_emp>d_perm)/para.nperm;
  case -1
    p = 1-sum(d_emp<d_perm)/para.nperm;
end
 

if p<para.alpha
  h = 1;
else
  h = 0;
end
      
    
    

  
  
  