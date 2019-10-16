a = 0.5 ;

t = 1000 ;

la_s = [1e-3:1e-3:1e-2 2e-2:1e-2:1e-1 2e-1:1e-1:1] ;

la_u = 600e-6;

k = la_s ./ la_u;

pa = 1 - (1 + (3.5.*k).^-1).^-3.5 ;
%Pc = exp(- a^-1 * pi .* pa .* la_s .* log(1 + t));
Pc = exp(- a^-1 * pi .* la_s .* log(1 + t));

plot(la_s , Pc)