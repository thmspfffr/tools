clear all;


% fixed params:

N = 4;

% wII=4;
% wIE=16;
% wEI=10;
% wEE=12;

% Connectivity:
wII=4;
wIE=16;
wEI=12;
wEE=12;

g = 1;

tauE = 1;
tauI = 2; %1
tau = [tauE;tauI;tauE;tauI];

Ii = -4;
Io=zeros(N,1);
Io(2) = Ii;
Io(4) = Ii;

dt=0.001;
tmax = 500; %5000;
tspan=0:dt:tmax;
L = length(tspan);

ds = 2;
Tds = length(0:ds*dt:tmax)-1;

sigma = 0.002;
sdt = sqrt(dt)*sigma;


% Connectivity:

W11 = [wEE -wEI; wIE -wII];
W22 = W11;
W12 = [g 0; 0 0];
W21 = [g 0; 0 0];


W = [W11 W12; W21 W22];

Gain = .20;
p = 1.5;

Ies = -4:.2:0;

numIes = length(Ies);

Rc     = zeros(numIes,4);
RcAmpl = zeros(numIes,4);

nTrials = 10;

uu = -4:.1:4;
Transf = zeros(length(uu),2,4);

for Q = 1:4  % conditions (no gain change, equivalent gain change, daE<daI, daE>daI)
    
% transfer functions:
% gains are given by 1/aE and 1/aI
if Q == 1
gE = 1;
gI = 1;
elseif Q == 2
gE = 1+Gain;
gI = 1+Gain;    
elseif Q == 3
gE = 1+p*Gain;
gI = 1+Gain;
elseif Q == 4
gE = 1+Gain;
gI = 1+p*Gain;    
end

aE = 1/gE;
aI = 1/gI;

Fe = @(x) 1./(1 + exp(-x/aE) );
Fi = @(x) 1./(1 + exp(-x/aI) );
F = @(x) [feval(Fe,x(1));feval(Fi,x(2));feval(Fe,x(3));feval(Fi,x(4))];

% store functions:
Transf(:,1,Q) = feval(Fe,uu);
Transf(:,2,Q) = feval(Fi,uu);
    

 for k = 1:numIes
    
    Ie = Ies(k);
    Io(1) = Ie;
    Io(3) = Ie;                       
    
    for trial = 1:nTrials
        
        r = rand(N,1);
        R = zeros(Tds,N);
        tt = 0;
        
        % Warm-up:
        for t = 1:5000
            u = W*r + Io;            
            K = feval(F,u);            
            r = r + dt*(-r + K)./tau + sdt*randn(N,1);
        end            
        
        for t = 1:L            
            u = W*r + Io;            
            K = feval(F,u);            
            r = r + dt*(-r+ K)./tau + sdt*randn(N,1);            
            if mod(t,ds)==0
                tt=tt+1;
                R(tt,:)=r;
            end            
        end
        
         RstatE = R(:,[1 3]);
         Ampl(:,1)=abs(R(:,1) + 1i*R(:,2));
         Ampl(:,2)=abs(R(:,3) + 1i*R(:,4));
        
        rc = corr(RstatE);
        Rc(k,Q) = Rc(k,Q) + rc(2)/nTrials;
        
        rc = corr(RstatE);
        RcAmpl(k,Q) = RcAmpl(k,Q) + rc(2)/nTrials;
        
        
    end
        
    
 end

end

figure
xSize = 15; ySize = 4.8;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperUnits','centimeters')
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[50 50 xSize*50 ySize*50],'Color','w')

axes('position',[.07 .2 .22 .69])
plot([Ies(1) Ies(end)],[0 0],'k:')
hold on
plot(Ies,Rc,'.-')
xlabel('Control param. \itI_{E}','fontname','times')
ylabel('corr. EE (c_{E-E })','fontname','times')
l = {'no gain','\deltaG_{E} = \deltaG_{I}',...
    '\deltaG_{E} > \deltaG_{I}','\deltaG_{E} < \deltaG_{I}'};
col=lines(4);
set(gca,'ylim',[-.4 1.001],'fontsize',9)
xlim=get(gca,'xlim');
ylim=get(gca,'ylim');
for i=1:4
   plot(xlim(1)+.83*(xlim(2)-xlim(1)),ylim(1)+(.38-(i-1)*.1)*(ylim(2)-ylim(1)),...
       'o','color',col(i,:),'markerfacecolor',col(i,:))  
   text(.9,.4-(i-1)*.11,l{i},'units','normalized','backgroundcolor','w','fontsize',8)
end
plot(-1.9*[1 1],ylim,'k:')
box off
text(-.22,1.05,'A','units','normalized','fontsize',14)


Q = 2;
Iopt = -3.2;
% Iopt = find( Rc(:,1)-Rc(:,Q)>.05 & Rc(:,1)<.7 & Rc(:,Q)>0 ); 

axes('position',[.405 .6 .19 .3])
plot(uu,Transf(:,1,1),'b')
hold on
plot(uu,Transf(:,1,Q),'k')
ylabel('F_{E}  ','fontname','times','rotation',0)
set(gca,'fontsize',9)

text(-.25,1.1,'B','units','normalized','fontsize',14)

axes('position',[.405 .2 .19 .3])
plot(uu,Transf(:,2,1),'b')
hold on
plot(uu,Transf(:,2,Q),'k')
set(gca,'fontsize',9)
xlabel('Input','fontname','times')
ylabel('F_{I}  ','fontname','times','rotation',0)


aE = 1;
aI = 1;

Fe = @(x) 1./(1 + exp(-x/aE) );
Fi = @(x) 1./(1 + exp(-x/aI) );
F = @(x) [feval(Fe,x(1));feval(Fi,x(2));feval(Fe,x(3));feval(Fi,x(4))];

    Ii = Iopt; %Ies(Iopt(1));
    Io(1) = Ie;
    Io(3) = Ie;                       

        r = rand(N,1);
        R = zeros(Tds,N);
        tt = 0;      
        for t = 1:L           
            u = W*r + Io;            
            K = feval(F,u);            
            r = r + dt*(-r./tau + K) + sdt*randn(N,1);            
            if mod(t,ds)==0
                tt=tt+1;
                R(tt,:)=r;
            end            
        end
        
         RstatE = R(tt-1999:tt,[1 3]);
         
axes('position',[.68 .6 .3 .3])

time = ds*dt*(0:size(RstatE,1)-1);
plot(time,RstatE)
text(-.22,1.1,'C','units','normalized','fontsize',14)

aE = 1/(1+p*Gain);
aI = 1/(1+Gain);    

Fe = @(x) 1./(1 + exp(-x/aE) );
Fi = @(x) 1./(1 + exp(-x/aI) );
F = @(x) [feval(Fe,x(1));feval(Fi,x(2));feval(Fe,x(3));feval(Fi,x(4))];

    Ie = Iopt; %Ies(Iopt(1));
    Io(1) = Ie;
    Io(3) = Ie;                       

        r = rand(N,1);
        R = zeros(Tds,N);
        tt = 0;      
        for t = 1:L           
            u = W*r + Io;            
            K = feval(F,u);            
            r = r + dt*(-r./tau + K) + sdt*randn(N,1);            
            if mod(t,ds)==0
                tt=tt+1;
                R(tt,:)=r;
            end            
        end
        
         RstatE = R(tt-1999:tt,[1 3]);
         
axes('position',[.68 .2 .3 .3])

time = ds*dt*(0:size(RstatE,1)-1);
plot(time,RstatE)

text(.4,.8,sprintf('I_{E} = %2.2f',Ii),'units','normalized')

% export_fig Fig_TwoWC_Correlations_Ie.tiff -r300


