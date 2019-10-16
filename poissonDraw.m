ve  = N_points(2,:)' / (300/4) ;

 h = histogram(ve,'Normalization','probability');
 
 hold on;
 
 pd = fitdist(ve,'Poisson');
 
 x_values = h.BinLimits(1):h.BinLimits(2);
y = pdf(pd,x_values);
plot(x_values,y,'LineWidth',2)