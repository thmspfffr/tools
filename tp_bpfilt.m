function y = tp_bpfilt(x,Fbp,Fs,N,filtertype)

Fn = Fs/2;

switch filtertype
  case 'fir'
      if isempty(N)
        N = 3*fix(Fs / Fbp(1));
      end
      if N > floor( (size(x,1) - 1) / 3)
        N=floor(size(x,1)/3) - 1;
      end
      [B, A] = fir1(N, [min(Fbp)/Fn max(Fbp)/Fn]);
  case 'butter'
    if isempty(N)
      N = 4;
    end
    [B, A] = butter(N, [min(Fbp)/Fn max(Fbp)/Fn]);
end
    


meandat = mean(x);
x = bsxfun(@minus, x, meandat);

y = filter_with_correction(B,A,x,'twopass-average');
