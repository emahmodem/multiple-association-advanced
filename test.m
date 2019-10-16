figure
x = 0:0.01:20;
y1 = 200*exp(-0.05*x).*sin(x);
y2 = 0.8*exp(-0.5*x).*sin(10*x);
[AX,H1,H2] = plotyy(x,y1,x,y1/100,'plot');
set(get(AX(1),'Ylabel'),'String','Signal') 
set(get(AX(2),'Ylabel'),'String','Signal/100') 
set(H1,'LineStyle','--','Color','r')
set(H2,'LineStyle','none','Color','r')
set(AX,{'ycolor'},{'b';'b'})      % ...and to adjust the axis color