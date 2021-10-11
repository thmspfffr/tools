function p = tp_permtest(x,y,nperm,tail)

dat = [x(:) y(:)];

if sum(any(isnan(dat),2))>0
  warning(sprintf('Ignoring %d NaN entry/entries',sum(any(isnan(dat),2))))
  dat(any(isnan(dat),2),:)=[];
end

 
% para.nperm = 1000;
% para.alpha = 0.05;

tmp = clock;
seed = ((tmp(1)+tmp(2)*tmp(3))/tmp(4)+tmp(5))*tmp(6);
rng(seed,'twister')

all_idx1 = randi(2,[size(dat,1),nperm]);


for iperm = 1 : nperm
     
    % within subjects permutation test
%     fprintf('Perm #%d\n',iperm);
    
    idx1 = all_idx1(:,iperm);
    idx2 = 3-idx1;
    
    for i = 1 : length(idx1)
      permdat(i,1) = dat(i,idx1(i));
      permdat(i,2) = dat(i,idx2(i));
    end
    
    perm_diff(iperm) = mean(permdat(:,1)-permdat(:,2));
    
end

emp_diff = mean(dat(:,1)-dat(:,2));

if tail == 0
  % two-tailed
  p = 1-sum(abs(emp_diff)>abs(perm_diff))/nperm; 
elseif tail == 1 
  % "right": empirical mean is greater than permutation mean
  p = 1-sum(emp_diff>perm_diff)/nperm; 
elseif tail == -1
  % "left": empirical mean is less than permutation mean
  p = 1-sum(emp_diff<perm_diff)/nperm; 
else
  error('para.tail has to be -1, 0 or 1')
end
