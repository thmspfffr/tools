
function [lpc stdlpc]=tp_data2lpc_jackknife(data,para,A1,A2)
% calculates lagged phase coherence (lpc) and its standard deviation
% between all combinations of given voxels
%
% usage:  
% [lpc stdlpc]=data2lpc_jackknife(data,segleng,segshift,epleng,f,fsample,A1)
% input: 
% data:   NxM matrix for N time points and M channels
% segleng:  length of segment (in samples)
% segshift:  shift of segments (in samples)
% epleng: length of epoch (or trial) (in samples)
% f:   frequency of interest (in Hz)
% fsample: sampling frequeny (in Hz)
% A1: MxP matrix of one.dimensional filters for P voxels
% 
% outpot: 
% lpc:  PxP matrix
% stdlpc: PxP matrix, estimated standard deviation of lpc
% 
% remark: lpc is defined as abs(imag(coh)/sqrt(1-real(coh)^2)))
% for complex coherence coh. 

epleng    = para.epleng;
segleng   = para.segleng;
fsample   = para.fsample;
segshift  = para.segshift;
f         = para.foi;


nn=(1:segleng)'-segleng/2;
mywin=hanning(segleng);
s1=cos(nn*f*2*pi/fsample).*mywin;
s2=sin(nn*f*2*pi/fsample).*mywin;
ss=s1-sqrt(-1)*s2;

[nchan ns]=size(A1);
[n nchan]=size(data);
ss=repmat(ss,1,nchan);

nep=floor(n/epleng);
nseg=(epleng-segleng)/segshift+1;

cs=zeros(nchan,nchan);

  csloc=zeros(nchan,nchan,nep);
kk=0;
for i=1:nep;
   dloc=data((i-1)*epleng+1:i*epleng,:);
 
    for j=1:nseg
        kk=kk+1;
        dloc2=dloc((j-1)*segshift+1:(j-1)*segshift+segleng,:);
        dataf=transpose(sum(dloc2.*ss));
        xx=dataf*dataf';
        csloc(:,:,i)=csloc(:,:,i)+xx;
        cs=cs+xx;
    end
    csloc(:,:,i)=csloc(:,:,i)/nseg;
    
end

cs=cs/kk;
css=A1'*cs*A1/kk;
coh=css./sqrt(diag(css)*diag(css)');
lpc=abs(imag(coh)./sqrt(1-real(coh).^2));

lpclocm=zeros(ns,ns);
lpclocvar=zeros(ns,ns);

for i=1:nep;
    csx=A1'*((nep*cs-csloc(:,:,i))/(nep-1))*A1;
    cohloc=csx./sqrt(diag(csx)*diag(csx)');
    lpcloc=abs(imag(cohloc)./sqrt(1-real(cohloc).^2));
    lpclocm=lpclocm+lpcloc;
    lpclocvar=lpclocvar+lpcloc.^2;
 end
lpclocm=lpclocm/nep;
lpclocvar=lpclocvar/nep;

stdlpc=sqrt(lpclocvar-lpclocm.^2)*sqrt(nep);




return;