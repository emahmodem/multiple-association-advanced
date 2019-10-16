a = 0.05 ;
a = 1.037;
t_dB= 5;

t = 10^(t_dB/10);

la_s = [1e-3:1e-3:1e-2 2e-2:1e-2:1e-1 2e-1:1e-1:1] ;

la_u = 600e-6;

k = la_s ./ la_u;

pa = 1 - (1 + (3.5.*k).^-1).^-3.5 ;
%Pc = exp(- a^-1 * pi .* pa .* la_s .* log(1 + t));
R_2 = a ./ (pi * pa .* la_s);

R = pa .* la_s * 1e6 .* log(1 + t).* exp(- pi/a .*  pa .* la_s .* log(1 + t)) ;
semilogx(la_s , R)
hold on;