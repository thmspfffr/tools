function pout=fdr1(p,q,method)
%
% calculates pvalue for a false discovery rate (fdr)
%
% usage: pout=fdr(p,q,method)
%
% input:  p is a column vector of p-values; q is the desired false 
% detection rate (i.e. the proportion of false positives out of all 
% positives; typically q=0.05); and method=0 for conservative (pessimistic) 
% assumptions about correlations between voxels and method=1 for optimistic 
% assumptions. 
%
% outpout: pout is the p value
%
% This is based on the paper "Thresholding of statistical maps in functional
% neuroimaging using the false discovery rate" by C.R. Genovese at al, 
% NeuroImage 15, pp. 870-878 (2002). 
%
% Guido Nolte 4/30/03 
%

if nargin<3;method=1;end

[n,m]=size(p);

if m>1
    error('p must be a column vector')
end


if q<0 | q>1
    error('second argument must be between 0 and 1');
end


if method==1
    %disp('method is optimistic');
    c=1;
elseif method==0
     %   disp('method is pessimistic');

    c=sum (1./(1:n)');
else
    error('third argument must be either 0 (pessimistic) or 1 (optimistic)')
end

ps=sort(p);
x=(1:n)'*q/(n*c);
ii=(1:n)';

iip=ii(ps<= x);

[n1,m1]=size(iip);

if n1==0
    %disp('warning: there was no solution, pout set to bonferoni level')
    pout=q/n;
else
   pout=ps(max(iip));
   if pout<q/n
       pout=q/n;
       %disp(' pout set to bonferoni level')
   end
   
end



return;
