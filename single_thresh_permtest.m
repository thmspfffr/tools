function [p,max_t] = w

para.nperm = 10000;
para.alpha = 0.05;
para.tail = 'both'; % alternatives left/right
% X is the input variable, e.g. N_Times x N_Freqs x N_Subj

for iperm = 1 : para.nperm
  
  fprintf('Permutation %d out of %d...\n',iperm,para.nperm)
  signs = sign(rand(size(X,3),1)-0.5);
  
  permdat = X.*permute(repmat(signs,[1 size(X,1) size(X,2)]),[2 3 1]);
  
  [~,~,~,s]=ttest(zeros(size(permdat)),permdat,'dim',3,'tail',para.tail,'alpha',para.alpha);
  
  max_t(iperm) = max(abs(s.tstat(:)));
  
end

[~,~,~,s]=ttest(zeros(size(permdat)),X,'dim',3,'tail',para.tail,'alpha',para.alpha);

emp_t = s.tstat;

if strcmp(para.tail,'both')
  emp_t = abs(emp_t);
  for i1 = 1:size(X,1)
    for i2 = 1:size(X,2)
      p(i1,i2)=1-sum(emp_t(i1,i2) > max_t)/para.nperm;
    end
  end
  
else
  error('Not supported')
end







