M = 10
%la = [1e-3:1e-3:1e-2 2e-2:1e-2:1e-1 2e-1:1e-1:1] ;
la = [100:10:1000]*1e-6
f = zeros(numel(M),numel(la));
g = zeros(numel(M),numel(la));
for n = 1:M
    %f(n,:) = sqrt(n ./ (pi .* la))
    g(n,:) = (gamma(n+1/2) ./ gamma(n)) ./ sqrt(pi .* la)
end

%semilogx (la , f, 'r-' , la,g,'k--' )
figure;semilogx (la , g )