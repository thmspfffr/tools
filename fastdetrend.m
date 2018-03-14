function signal = fastdetrend(signal)

n = size(signal,1);
if n == 1
  % If a row, turn into column vector
 signal = signal(:);			
end


N = size(signal,1);
a = [zeros(N,1) ones(N,1)];
a(1:N) = (1:N)'/N;

signal = signal - a*(a\signal);

if(n==1)
    signal = signal.';
end
end