function [ B ] = hessiandelag(fname,hname,x, lambda )
% Aproximación numérica a la matriz:
%  hessiana(fname(x)) + sum_{j=1}^m \lambda(j)* hessian(hname_j(x))
%   
% fname: R^n --> R
% hname: R^n --> R^m
% lambda vector columna de dimensión m.

fx = feval(fname,x);
hx = feval(hname,x);
paso = 1.e-05;
n = length(x);
m = length(hx);
B = zeros(n);

for k = 1:n
    x1 = x; x1(k) = x1(k) + paso;
    fx1 = feval(fname,x1);   hx1 = feval(hname,x1);
    for j = 1:k
        x2 = x1; x2(j) = x2(j) + paso;
        fx2 = feval(fname,x2);  hx2 = feval(hname,x2);
        x3 = x; x3(j) = x3(j) + paso;
        fx3 = feval(fname,x3);  hx3 = feval(hname,x3);
        B(k,j) = ( fx2 -fx1 - fx3 + fx)/(paso^2);
        w = ( hx2 -hx1 - hx3 + hx)/(paso^2);
        w = w.*lambda;
        B(k,j) = B(k,j) + sum(w);
        if (k~=j)
            B(j,k) = B(k,j);
        end
    end
    
end


















end

