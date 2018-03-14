
function [lpc coh]=tp_lpc(data,para,A1,A2)

epleng    = para.epleng;
segleng   = para.segleng;
fsample   = para.fsample;
segshift  = para.segshift;
f         = para.foi;




scale=sqrt(mean(mean(data.^2)));
data=data/scale;
nn=(1:segleng)'-segleng/2;

mywin=hanning(segleng);
s1=cos(nn*f*2*pi/fsample).*mywin;
s2=sin(nn*f*2*pi/fsample).*mywin;
ss=s1-sqrt(-1)*s2;
if nargin<8;
    A2=A1;
end
[nchan ns1]=size(A1)
[nchan ns2]=size(A2)

[n nchan]=size(data);
ss=repmat(ss,1,nchan);

nep=floor(n/epleng);
nseg=(epleng-segleng)/segshift+1;

res1=zeros(ns1,ns1);
if nargin>7
res2=zeros(ns2,ns2);
res3=zeros(ns1,ns2);
end

kk=0;
for i=1:nep;
    dloc=data((i-1)*epleng+1:i*epleng,:);
    for j=1:nseg
        kk=kk+1;
        dloc2=dloc((j-1)*segshift+1:(j-1)*segshift+segleng,:);
        dataf=sum(dloc2.*ss);
        datasf1=transpose(dataf*A1);
         datasf2=transpose(dataf*A2);
         if kk==1;
             res1=datasf1*datasf1';
             if nargin>7
             res2=datasf2*datasf2';
             res3=datasf1*datasf2';
             end
         else
             res1=res1+datasf1*datasf1';
             if nargin>7
             res2=res2+datasf2*datasf2';
             res3=res3+datasf1*datasf2';
             end
         end
   
        
    end
end

res1=res1/kk;
if nargin>7
res2=res2/kk;
res3=res3/kk;
end

if nargin>7
coh=res3./sqrt(diag(res1)*diag(res2)');
else
    coh=res1./sqrt(diag(res1)*diag(res1)');
end
lpc=imag(coh)./sqrt(1-real(coh).^2+eps);



return;