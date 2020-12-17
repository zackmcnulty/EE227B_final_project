function [ J ] = Jacobiana( hname,x )
% Cálculo de la matriz< jacobiana de hname: R^n --> R^m en el punto x.
%
%   J(k,j) = parcial de h_k / pracial en x_k

hx = feval(hname,x);
n = length(x);
m = length(hx);
paso = 1.e-05;
J = zeros(m,n);

for j = 1:n
   xa = x; xa(j) = xa(j) + paso;
   hxa = feval(hname,xa);
   J(:,j)  = ( hxa- hx)/paso; 
end


end

