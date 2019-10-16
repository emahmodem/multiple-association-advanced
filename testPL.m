aa = 4;
a = 0.94;
b = 0.5; 

r = 1:300;

PL1 = exp(-a .* r .^ b );
PL2 = r.^-aa;

a = 5e-5;
b   = 2;
PL3 = exp(-a .* r .^ b );

semilogy(r,PL1,'r-',r,PL3,'r-.',r,PL2,'k--')