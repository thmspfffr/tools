function [signsout,yout]=adjustsign(y)
% adjusts signs of signals in each channel, such 
% the  correlation between pairs of signals is 
% (mostly) positive
% usage: [signsout,yout]=adjustsign(y)
% 
% input: 
% y:  nxm matrix for n time points and m channels
%
% output: 
% signsout: mx1 vector of signs
% yout: sign adjusted data

[n,nchan]=size(y);
signsout=zeros(nchan,1);
y2=detrend(y,'constant');
yout=y;
c=y2'*y2;
[u,s,v]=svd(c);

y2m=y*u(:,1);

for i=1:nchan;
    signsout(i)=sign(y2m'*y2(:,i));
    yout(:,i)=y(:,i)*signsout(i);
end

return;
