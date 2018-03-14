function powcorr = tp_orthopowcorr_fast(dat,para)

[n ns1] = size(dat);

% input: complex data, TxN (time x voxel)
kk = 0;

epleng    = para.epleng;
segleng   = para.segleng;
segshift  = para.segshift;

nep  = floor(n/epleng);
nseg = floor((epleng-segleng)/segshift+1);


for j = 1 : nseg
  
  fprintf('Seg %d/%d ...\n',i,nep,j,nseg)
  
  kk = kk + 1;
  dataf = dat((j-1)*segshift+1:(j-1)*segshift+segleng,:);
  
  
  for i1 = 1 : ns1
    
    x1 = dataf(i1);
    x2 = imag(dataf*conj(x1)./abs(x1));
    
    y1 = log(abs(x1)^2);
    y2 = log(abs(x2).^2);
    
    if kk == 1
      res1(i1,:) = y1*y2;
      res2(i1,:) = y1;
      res3(i1,:) = y2;
      res4(i1,:) = y1^2;
      res5(i1,:) = y2.^2;
    else
      res1(i1,:) = res1(i1,:) + y1*y2;
      res2(i1,:) = res2(i1,:) + y1;
      res3(i1,:) = res3(i1,:) + y2;
      res4(i1,:) = res4(i1,:) + y1^2;
      res5(i1,:) = res5(i1,:) + y2.^2;
    end
    
  end
  
  res1 = res1/kk;
  res2 = res2/kk;
  res3 = res3/kk;
  res4 = res4/kk;
  res5 = res5/kk;
  
  tmp = (res1-res2.*res3)./(sqrt((res4-res2.^2).*(res5-res3.^2))+eps);
  
  % average over the two directions, fisher z-transform first
  powcorr = tanh(triu(atanh(tmp)./2,1)+tril(atanh(tmp)./2,-1)')+[tanh(triu(atanh(tmp)./2,1)+tril(atanh(tmp)./2,-1)')]';
  
end
end