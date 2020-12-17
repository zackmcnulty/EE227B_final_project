function [ gx ] = gradiente( fname,x )
% Calcular el gradiente de fname: R^n --> R en el punto x.
% gx es el vector gradiente


n = length(x);
paso = 1.e-05;
fx = feval(fname,x);
gx = zeros(n,1);
for k = 1:n
    xa = x; xa(k) = x(k) + paso;
    fxa = feval(fname,xa);
    gx(k) = (fxa -fx)/paso;
end


end

