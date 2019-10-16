xv = [-300   0   0 -300 -300] ;
yv = [-300 -300  0  0   -300] ;

xq = [-20.7153];
yq = [-147.1237];

plot(xv,yv,'LineWidth',2) ; hold on;
plot(xq,yq,'ro')
[in,on] = inpolygon(xq,yq,xv,yv);