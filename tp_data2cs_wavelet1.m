
function [dataf, cs]=tp_data2cs_wavelet1(data,f,fsample)
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
[~,ff,opt]=tp_mkwavelet(f, .5, fsample);

segleng = round(400/((ff(2)-ff(1))/2));
segshift = floor(segleng / 2);
% segleng = opt.n_win;
% segshift = opt.n_shift;
epleng=size(data,1);

nn=(1:segleng)'-segleng/2;
mywin=hanning(segleng);
s1=cos(nn*f*2*pi/fsample).*mywin;
s2=sin(nn*f*2*pi/fsample).*mywin;
ss=s1-sqrt(-1)*s2;
ssout=ss;



[n nchan]=size(data);
ss=repmat(ss,1,nchan);

nep=floor(n/epleng);
nseg=(epleng-segleng)/segshift+1;

cs=zeros(nchan,nchan);
coh=cs;


kk=0;
for i=1:nep;
   dloc=data((i-1)*epleng+1:i*epleng,:);
    for j=1:nseg
        
        dloc2=dloc((j-1)*segshift+1:(j-1)*segshift+segleng,:);
        
        if ~any(isnan(dloc2(:,1)))
        
        dataf(:,j)=transpose(sum(dloc2.*ss));
         kk=kk+1;
         if kk==1;
             cs=dataf*dataf(:,j)';
         else
              cs=cs+dataf(:,j)*dataf(:,j)';
         end
         
        else
          
%           dataf(:,j) = nan(size(dloc2,1),size(dloc2,2));
          
        end
   
        
    end
end

cs=cs/kk;
coh=cs./sqrt(diag(cs)*diag(cs)');




return;