
function [resout]=tp_orthopowercorr(data,para,A1,A2)

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
% A1: NxM1 matrix, spatial filter for N sensors and M1 voxels for first set of voxels
% A2: NxM2 matrix, spatial filter for N sensors and M2 voxels for second set of voxels
%
% output
% resout: M1xM2 matrix of power correlations after orthogonalization.
%
% TP, 11/07/2014 
% keepsegs: keep correlation coefficients for each segment

if ~exist('para','var')
  para.segave = 1;
end

epleng    = para.epleng;
segleng   = para.segleng;
fsample   = para.fsample;
segshift  = para.segshift;
f         = para.foi;

scale = sqrt(mean(mean(data.^2)));
data  = data/scale;


  [~, ns1] = size(A1);
  [~, ns2] = size(A2);


[n nchan] = size(data);

  
%   case 'hanning'
    nn    = (1:segleng)'-segleng/2;
    mywin = hanning(segleng);
    s1 = cos(nn*f*2*pi/fsample).*mywin;
    s2 = sin(nn*f*2*pi/fsample).*mywin;
    ss = s1-sqrt(-1)*s2;
    ss = repmat(ss,1,nchan);

%   case 'wavelet'
%     win = pconn_wavelet(fsample,f,3);
%     segleng = size(win,1);
%     ss = repmat(win,1,nchan);
%     nn    = (1:segleng)'-segleng/2;
% end

[idx_triu(:,1),idx_triu(:,2)]=find(triu(ones(ns1,ns1),1));
[idx_tril(:,1),idx_tril(:,2)]=find(tril(ones(ns1,ns1),-1));

if nargin<8;
  A2 = A1;
end

nep  = floor(n/epleng);
nseg = floor((epleng-segleng)/segshift+1);


%datarawf=[];
perc = 0;
resout = zeros(ns1,ns1,nep);

for i = 1 : nep

  res1 = zeros(ns1,ns2);
  res2 = zeros(ns1,ns2);
  res3 = zeros(ns1,ns2);
  res4 = zeros(ns1,ns2);
  res5 = zeros(ns1,ns2);

  kk = 0;

  if nep ~= 1
    dloc = data((i-1)*epleng+1:i*epleng,:);
  else 
    dloc = data((i-1)*epleng+1:i*epleng,:); clear data
  end
  
  for j = 1 : nseg
      
    fprintf('Ep %d/%d / Seg %d/%d ...\n',i,nep,j,nseg)

    kk = kk + 1;
    dloc2 = dloc((j-1)*segshift+1:(j-1)*segshift+segleng,:);
    dataf = sum(dloc2.*ss);
    
    datasf1 = dataf*A1;
    datasf2 = dataf*A2;
    
    % CONVERT TO AAL
%     if pars.aal
%       datasf1 = dataf;
%       datasf2 = dataf;
%     
    
  end


  res1 = res1/kk;
  res2 = res2/kk;
  res3 = res3/kk;
  res4 = res4/kk;
  res5 = res5/kk;

  tmp = (res1-res2.*res3)./(sqrt((res4-res2.^2).*(res5-res3.^2))+eps);
  
  % average over the two directions, fisher z-transform first
  resout(:,:,i) = tanh(triu(atanh(tmp)./2,1)+tril(atanh(tmp)./2,-1)')+[tanh(triu(atanh(tmp)./2,1)+tril(atanh(tmp)./2,-1)')]';
  
  
  clear tmp 
  
end


return;