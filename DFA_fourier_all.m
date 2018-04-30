function [fluct,slope]=DFA_fourier_all(DATA,win_length,mymethod)
% calculates fluctions for DFA   and slopes for specified window lengths
% usng the stationary method in the Fourier domain.  
% Input:
% Data: nxm matrox for n time points and m channels 
% win_length: kx1 vector of window lengths. Unit is number of samples. 
%             Note, that the values don't have to be integers. 
% mymethod (optional). mythod is either 'boxcar' (the default), which
%                     corresponds to the classical approach of detrending in a window,
%                     or 'gauss', which means that detrending is done with
%                     a gaussian weight function of a width which related
%                     to win_length
% output:
% fluct: kxm matrix of fluctuations for k window lengths and m channels
% slope: kxm matrix of slopes for k window lengths and m channels
% 
% trivial example for white noise in 20 channels: 
% data=randn(1000000,20);
% w=10*2.^(0:10);
% [fluct,slope]=DFA_fourier_all(data,w); %using boxcar
% figure; loglog(w,fluct)
% figure;semilogx(w,slope)
% 
% % or
% [fluct,slope]=DFA_fourier_all(data,w,'gauss'); %using gauss window
% figure; loglog(w,fluct)
% figure;semilogx(w,slope)


if nargin<3
    mymethod='boxcar';
end
 [n,nchan]=size(DATA);
 isodd=bitget(n,1); %is 1 if n is odd, and zero if n is even
nw=length(win_length); 
 for i=1:nchan;
     DATA(:,i)=DATA(:,i)-mean(DATA(:,i));
 end
 
output=zeros(nw,nchan);
slope=zeros(nw,nchan);
       
    %   DATA=cumsum(DATA);

    
      DATAF=fft(DATA);

      if isodd==1
          nx=(n+1)/2;
      else
          nx=n/2+1;
      end
      
      
      DATAp=2*abs(DATAF(2:nx,:)).^2;
      if isodd==0;
      DATAp(nx-1,:)=DATAp(nx-1,:)/2;
      end
   
      
      if strcmp(mymethod,'boxcar')
      ff=(1:nx-1)';
      g1=sin(pi*ff/n);
    
      
      for k=1:nw;
          wl=win_length(k);
         
           hsin=sin(pi*ff*wl/n);
       hcos=cos(pi*ff*wl/n);
      hx=(1-hsin./(wl*g1));
      h=hx./(2*g1);
      h2=h.^2; 
        h2=repmat(h2,1,nchan);
      F2=sum(h2.*DATAp);
      output(k,:)=sqrt(F2)/n;
       
      
          
      hy=-hx.*(hcos*pi.*ff/n-hsin/wl)./(wl*g1);
      h3=hy./(4*g1.^2);
      h3=repmat(h3,1,nchan);
    
    
       slope(k,:)=sum(h3.*DATAp)./F2*wl;
       
      end
      
      elseif strcmp(mymethod,'gauss');
        % DATAp=1;          
            ff=(1:nx-1)'/n;
            ff2=ff.^2;
            df=1/n;
      for k=1:nw;
          wl=win_length(k);
         sigma=wl/sqrt(12);
         h=exp(-2*pi^2*sigma^2*ff2);
         h2=(1-h).^2./(ff2);

         h2x=repmat(h2,1,nchan);
         F2=sum(h2x.*DATAp)*df/(4*pi^2);

         hx=(1-h).*h;
         hx=repmat(hx,1,nchan);
         D=df*sum(hx.*DATAp)/(4*pi^2);

         slope(k,:)=sigma^2*4*pi^2.*D./F2;
         output(k,:)=sqrt(F2)/sqrt(n);

       
      end
       elseif strcmp(mymethod,'hanning');
               ff=(1:nx-1)';
      g1=sin(pi*ff/n);
    
      
      for k=1:nw;
          wl=win_length(k);
         
           hsin=sin(pi*ff*wl/n);
       hcos=cos(pi*ff*wl/n);
      hx=(1-hsin./(wl*g1));
      h=hx./(2*g1);
      h2=h.^2; 
        h2=repmat(h2,1,nchan);
      F2=sum(h2.*DATAp);
      output(k,:)=sqrt(F2)/n;
%        
%       
%           
%       hy=-hx.*(hcos*pi.*ff/n-hsin/wl)./(wl*g1);
%       h3=hy./(4*g1.^2);
%       h3=repmat(h3,1,nchan);
%     
%     
%        slope(k,:)=sum(h3.*DATAp)./F2*wl;
       
      end
          
      end
      
  fluct=output;      
       return;