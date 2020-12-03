function y = tp_drawsquare(x1,x2,y1,y2)

line(gca,[x1 x2],[y2 y2],'color','k')
line(gca,[x1 x2],[y1 y1],'color','k')
line(gca,[x1 x1],[y1 y2],'color','k')
line(gca,[x2 x2],[y1 y2],'color','k')
