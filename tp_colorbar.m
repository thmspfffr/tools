function c = tp_colorbar(varargin);  

fig = get(gca);

c = colorbar; 
c.Box = 'off'; 

c.Position=[fig.Position(1)+fig.Position(3)-0.02 c.Position(2) 0.7*c.Position(3) c.Position(4)]; 
c.FontSize = 7;
if c.Limits(1)<0
  c.Ticks = [c.Limits(1) 0 c.Limits(2)];
  c.TickLabels = {num2str(c.Limits(1)); num2str(0); num2str(c.Limits(2))};
  c.Label.String = varargin;
c.Label.Rotation = 270;
else
   c.Ticks = [c.Limits(1) c.Limits(2)];
   c.TickLabels = {num2str(c.Limits(1));  num2str(c.Limits(2))};
end
