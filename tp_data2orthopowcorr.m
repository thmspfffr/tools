
function [resout]=tp_data2orthopowcorr(data,segleng,segshift,epleng,f,fsample,filt)
% calculates power correlation after orthogonalization between two sets of brain
% voxels. Usually, these two sets are either identical (then the rsult is
% for all pairs) or one set consists of a single voxel which is taken as
% seed. Data are epoched, with all epeochs (i.e. trials) stacked on top of each other.
% Each epoch is further divided into segments. Coefficients in the Fourier-domain
% are calculated for each segment. Thus the segment length determines the
% frequency resolution.
%
% usage:
% [resout]=data2orthopowercorr(data,segleng,segshift,epleng,f,fsample,A1,A2)
%
% input:
% data: TxN data matrix for T time points and N channels
% segleng: length of a segment (in number of samples)
% segshift: amount of shift (in number of samples) (explanation: segments
%           are shifted within but not across epochs. E.g.
%           segshift=floor(segleng/2) leads to segments of  50% overlap
%            (or slightly less if segleng is odd))
% epleng: length of epoch (or trial) in number of trials)
% f:      frequency of interest (in Hz)
% fsample: sampling rate

% output
% resout: M1xM2 matrix of power correlations after orthogonalization.

scale=sqrt(nanmean(nanmean(data.^2)));
data=data/scale;
nn=(1:segleng)'-segleng/2;

mywin=hanning(segleng);
s1=cos(nn*f*2*pi/fsample).*mywin;
s2=sin(nn*f*2*pi/fsample).*mywin;
ss=s1-sqrt(-1)*s2;
% if nargin<8;
%     A2=A1;
% end
[nchan ns1]=size(filt);
% [nchan ns2]=size(A2);

[n nchan]=size(data);
ss=repmat(ss,1,nchan);

nep=floor(n/epleng);
nseg=(epleng-segleng)/segshift+1;

res1=zeros(ns1,ns1);
res2=zeros(ns1,ns1);
res3=zeros(ns1,ns1);
res4=zeros(ns1,ns1);
res5=zeros(ns1,ns1);

kk=0;

%datarawf=[];
for i=1:nep
%   i=i
  dloc=data((i-1)*epleng+1:i*epleng,:);
  
  % idenfity segments containing NaNs and exclude from analysis
  % ---------------------------
  for i=1:nseg
    tmp=dloc((i-1)*segshift+1:(i-1)*segshift+segleng,1);
    nanidx(i) = any(isnan(tmp));
  end
  idx = find(~nanidx);
  % ---------------------------
  
  for j=1:length(idx)
    j
    kk=kk+1;
    dloc2=dloc((idx(j)-1)*segshift+1:(idx(j)-1)*segshift+segleng,:);
    
    if any(isnan(dloc2(:,1)))
      warning('NaN detected')     
      continue
    end
    
    dataf=sum(dloc2.*ss);
    datasf1=dataf*filt;
    %datarawf=[datarawf;datasf1]
    
    for i1=1:ns1;
      
      x1=datasf1(i1);
      x2=imag(datasf1*conj(x1)./abs(x1));
      y1=abs(x1)^2;
      y2=abs(x2).^2;
      if kk==1;
        res1(i1,:)=y1*y2;
        res2(i1,:)=y1;
        res3(i1,:)=y2;
        res4(i1,:)=y1^2;
        res5(i1,:)=y2.^2;
      else
        res1(i1,:)=res1(i1,:)+y1*y2;
        res2(i1,:)=res2(i1,:)+y1;
        res3(i1,:)=res3(i1,:)+y2;
        res4(i1,:)=res4(i1,:)+y1^2;
        res5(i1,:)=res5(i1,:)+y2.^2;
      end
    end
    
  end
end

res1=res1/kk;
res2=res2/kk;
res3=res3/kk;
res4=res4/kk;
res5=res5/kk;

resout=(res1-res2.*res3)./(sqrt((res4-res2.^2).*(res5-res3.^2))+eps);

return;