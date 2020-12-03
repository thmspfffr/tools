function [Itrue,Iboot] = ckmi_gg_dfi_ak(x, y,yshuffle)

% HACKED MI FUNCTION WITHOUT BIAS AND RANDOMIZATION INPUT FOR Y


ln2 = log(2);

Ntrl = size(x,1);
Nvarx = size(x,2);
Nvary = size(y,2);

x = bsxfun(@minus,x,sum(x,1)/Ntrl);
y = bsxfun(@minus,y,sum(y,1)/Ntrl);

%----------------------------------------------------
% joint variable
xy = [x y];
Cxy = (xy'*xy) / (Ntrl - 1);
% submatrices of joint covariance
Cx = Cxy(1:Nvarx,1:Nvarx);
ystart = Nvarx + 1;
Nvarxy = Nvarx + Nvary;
Cy = Cxy(ystart:Nvarxy,ystart:Nvarxy);

chCx = chol(Cx);
chCy = chol(Cy);
chCxy = chol(Cxy);

% entropies in nats
HX = sum(log(diag(chCx))); % + 0.5*Nvarx*log(2*pi*exp(1));
HY = sum(log(diag(chCy))); % + 0.5*Nvary*log(2*pi*exp(1));
HXY = sum(log(diag(chCxy))); % + 0.5*(Nvarx+Nvary)*log(2*pi*exp(1));
% convert to bits
Itrue = (HX + HY - HXY) / ln2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bootstrap
% Nboot = size(yshuffle,2);
matsize = size(yshuffle);
Nboot = matsize(end);
Iboot = zeros(1,Nboot);

if Nboot
  for r=1:Nboot
    % loop each element of yshuffle
    if numel(matsize)==2
        y = yshuffle(:,r);
    elseif numel(matsize)==3
        y = yshuffle(:,:,r);
    end
    xy = [x y];
    Cxy = (xy'*xy) / (Ntrl - 1);
    % submatrices of joint covariance
    Cx = Cxy(1:Nvarx,1:Nvarx);
    ystart = Nvarx + 1;
    Nvarxy = Nvarx + Nvary;
    Cy = Cxy(ystart:Nvarxy,ystart:Nvarxy);
    
    chCx = chol(Cx);
    chCy = chol(Cy);
    chCxy = chol(Cxy);
    
    % entropies in nats
    HX = sum(log(diag(chCx))); % + 0.5*Nvarx*log(2*pi*exp(1));
    HY = sum(log(diag(chCy))); % + 0.5*Nvary*log(2*pi*exp(1));
    HXY = sum(log(diag(chCxy))); % + 0.5*(Nvarx+Nvary)*log(2*pi*exp(1));
    % convert to bits
    Iboot(r) = (HX + HY - HXY) / ln2;
    
  end
end

return;

