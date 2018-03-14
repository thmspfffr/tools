function tp_signbar(x,y,ticklength,sign)
% plot significance bar
% tp_signbar(x,y,ticklength,sign)
% use: x = x(1)/x(2) coordinates
% y = y(1)/y(2) coordinates
% ticklength: ticks in percent of scale

g = get(gca);

ticklength = ticklength*mean(g.YLim)./2;

line([x(1) x(2)],[y y],'linewidth',1,'color','k')

line([x(1) x(1)],[y-ticklength y],'linewidth',1,'color','k')
line([x(2) x(2)],[y-ticklength y],'linewidth',1,'color','k')

if strcmp(sign,'*')
  text((x(1)+x(2))./2, y+0.015*mean(g.YLim)./2,'*','horizontalalignment','center')
elseif strcmp(sign,'**')
  text((x(1)+x(2))./2, y+0.015*mean(g.YLim)./2,'**','horizontalalignment','center')

elseif strcmp(sign,'***')
    text((x(1)+x(2))./2, y+0.015*mean(g.YLim)./2,'***','horizontalalignment','center')

else
  warning('Really???')
end
