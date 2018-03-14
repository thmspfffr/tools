function cmap = cmap_pm3d(x)

% x=nvals

x=linspace(0,1,x);
r = sqrt(x);
g = x.^3;
b=sin(x*2*pi);
b(b<0)=0;
cmap =[r;g;b]';
