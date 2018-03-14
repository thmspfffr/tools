function [cs coh ssout]=tp_data2cs_wavelet(data,segleng,segshift,epleng,f,fsample,wavelet)
% calculates cross-spectrum and coherence based on a wavelet using
% a Hanning window for given frequeny
%
% usage:  
% [cs coh]=data2cs_wavelet(data,segleng,segshift,epleng,f,fsample);
% input: 
% data:   NxM matrix for N time points and M channels
% segleng:  length of segment (in samples)
% segshift:  shift of segments (in samples)
% epleng: length of epoch (or trial) (in samples)
% f:   frequency of interest (in Hz)
% fsample: sampling frequeny (in Hz)
%
% outpot: 
% cs: cross-spectrum 
% coh: coherency (complex)
% ss: the complex wavelet
    [n nchan]=size(data);

switch wavelet
  case 'hanning'

    nn=(1:segleng)'-segleng/2;
    mywin=hanning(segleng);
    s1=cos(nn*f*2*pi/fsample).*mywin;
    s2=sin(nn*f*2*pi/fsample).*mywin;
    ss=s1-sqrt(-1)*s2;
    ssout=ss;

    ss=repmat(ss,1,nchan);
    
  case 'ft'
    
    para.fsample = fsample;
    para.freqoi = f;
    para.gwidth = 3;
    para.width  = 7;
    
    w = tp_mkwavelet(para);
    ssout = w;
    ss=repmat(w,1,nchan);

    segleng = length(w);
    segshift = floor(segleng/2);
    
end

nep=floor(n/epleng);
nseg=(epleng-segleng)/segshift+1;

cs=zeros(nchan,nchan);
coh=cs;


kk=0;
for i=1:nep;
   dloc=data((i-1)*epleng+1:i*epleng,:);
    for j=1:nseg
        kk=kk+1;
        dloc2=dloc((j-1)*segshift+1:(j-1)*segshift+segleng,:);
        dataf=transpose(sum(dloc2.*ss));
        
         if kk==1;
             cs=dataf*dataf';
         else
              cs=cs+dataf*dataf';
         end
   
        
    end
end

cs=cs/kk;
coh=cs./sqrt(diag(cs)*diag(cs)');




return;