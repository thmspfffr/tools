% wilson cowan

% simulates single wilson cowan unit
% parameters from deco et al. (2009) pnas
close all
clear all

fsample    = 1000;
time_steps = 100; % in seconds
dt         = 1/fsample;
% --------------------
% response function
% --------------------
b     = 40;
c     = 50;
a_exc = 0.1;
a_inh = 0.3;
% --------------------
% model parameters
% --------------------
te  = 1;    % tau(exc.)
ti  = 0.2;  % tau(inh.)
I_b = 30;   % background input
we  = 1.5;  % self-exc.
wi  = 1.5;  % inh. weight
sig = 0.2;    % noise
g   = 1;    % global coupling
v   = 2.5;
% --------------------
% get coupling
% --------------------
load    ~/Downloads/SC_AAL_90.mat    
areas = 1;
C     = SC/max(max(SC))*0.2;
% --------------------
% initial values
% --------------------
x   = 0.1*randn(1,areas);
y   = 0.1*randn(1,areas);
inp = 1;
% --------------------
%

% compute distances
d = compute_distance;
v = d./v;
T = round(v);

%%
fprintf('Initializing model ...\n')

for it = 1 : time_steps/dt
  
  ice = we*x - y + I_b;
  ici = wi*x;
  % excitatory pool
  x   = x + dt/te * (-x + inp * resp_fun(ice,a_exc,b,c) ) + sqrt(dt)*sig*randn(1,1);
  % inhibitory pool
  y   = y + dt/ti * (-y + inp * resp_fun(ici,a_inh,b,c) ) + sqrt(dt)*sig*randn(1,1);
  
  xs(it) = x;
  ys(it) = y;
  
end

figure_white; hold on
subplot(2,1,1); plot(xs)

segleng   = 4000;
segshift  = 200;
nseg      = floor((length(xs)-segleng)/segshift+1);


for iseg = 1 : nseg
  [p(:,iseg),f]=pwelch(xs((iseg-1)*segshift+1:(iseg-1)*segshift+segleng),hanning(4000),50,1:.1:20,400);
end

subplot(2,1,2); imagesc(flipud(log10(p))); colormap(jet);
xlabel('Time [ms]'); ylabel('Freq [Hz]');
set(gca,'YTick',[1:40:191],'YTickLabel',[20:-4:2])

%%
xs = xs(20000:end);
siginfo = nbt_Info;

siginfo.converted_sample_frequency = 400;

x=abs(hilbert(nbt_filter_fir(xs',8,12,siginfo.converted_sample_frequency,2/8)));

%%
% RUN SIMULATION
% T(i,j) = ;

for it = 1 : time_steps/dt
  it
  for i = 1 : areas
    
    C(i,i) = we/g;
    T(i) = 0;
        
    ice = I_b + sum(g*C(i,:).*x(it-T(i),i)') - y(it,i);
    ici = wi*x(it,i);

    % excitatory pool
    x(it+1,i)   = x(it,i) + dt/te * (-x(it,i) + inp * resp_fun(ice,a_exc,b,c) ) + sqrt(dt)*sig*randn(1,1);
    % inhibitory pool
    y(it+1,i)   = y(it,i) + dt/ti * (-y(it,i) + inp * resp_fun(ici,a_inh,b,c) ) + sqrt(dt)*sig*randn(1,1);

%     xs(it,i) = x(i);
%     ys(it,i) = y(i);
    
  end
  
end

fc = corrcoef(x);
%% 

win  =hanning(400);

segleng   = 400;
segshift  = 200;
nseg      = floor((length(xs)-segleng)/segshift+1)

for iseg = 1 : nseg
  iseg
  p(:,iseg)=pwelch(xs((iseg-1)*segshift+1:(iseg-1)*segshift+segleng),win,50,1:1:100,fsample);

end

figure;
imagesc(log10(p))
  
  


