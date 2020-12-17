function [x, iter ] = PCSBFGS( fname, hname,x )
% Programación cuadrática sucesiva con búsqueda de línea (función L1 de mérito)
% y actualización recortada de BFGS para el problema
%
% Min fname(x)
% S. A. hname(x) = 0 
%  
% donde
% fname: R^n --> R
% hname: R^n --> R^m
% son funciones no lineales

tol = 1.e-08;
maxiter = 50;
iter = 0;
C = 1;              % Primer parámetro para la función de mérito
gx = gradiente(fname,x);
hx = feval(hname,x);   J = Jacobiana(hname,x);
% Jacobiana es una matriz nxm
n = length(x);  m = length(hx);
r = rank(J);
if( r==m)
    lambda = - (J*J')\ (J*gx);
else
   lambda = ones(m,1);
end
B = eye(n);
vparo = gx + J'*lambda;
v = [vparo; hx];

while(norm(v) > tol && iter < maxiter)
    K = [B J'; J zeros(m,m)];
    ld = -[gx;hx];
    pp = K\ld;
    p = pp(1:n);
    plam = pp(n+1:n+m);
    %----------------------------------------------
    % Búsqueda de línea con la función de mérito L1
    if (norm(hx,1)>1.e-10)
        aux = (p'*B*p)/norm(hx,1);
        aux = abs(aux) + norm(plam,inf)+1;
        C = max([ aux   C ]);  % parámetro para la función de Mérito
    end
    phik = feval(fname,x)+C*norm(hx,1);
    dphi = -p'*B*p-hx'*plam-C*norm(hx,1); % derivada direccional
    alfa = 1;
    jter = 0;
    xaux = x + alfa*p;
    phiaux = feval(fname,xaux)+ C*norm(feval(hname,xaux),1);
    while( phiaux > phik + alfa*(1.e-04)*dphi && jter < 8)
        alfa = alfa/2;
        xaux = x + alfa*p;
        phiaux = feval(fname,xaux)+ C*norm(feval(hname,xaux),1);
        jter = jter + 1;
    end
%     if (jter == 8)
%         alfa = 1;
%     end
    %-------------------------------------
    % BFGS recortado
    xp = x + alfa*p;
    s = xp - x;
    gxp = gradiente(fname,xp);
    J2  = Jacobiana(hname,xp);
    w1 = gxp + J2'*plam;
    w2 = gx +J'*plam;
    y = w1 - w2;
    z = s'*B*s;
    

    if( s'*y >= (0.2)*z )
        theta = 1;
    else
        theta = ( (0.8)*z )/(z - (0.2)*s'*y) ;
    end
    r = theta*y + (1-theta)*(B*s);
    %
    B = B - ((B*s*s'*B)/z ) + (( r*r')/(s'*r));
    %--------------------------
    x = xp;
    J = J2;
    gx = gxp;
    hx = feval(hname,x);
    r = rank(J);
    if( r==m)
       lambda = - (J*J')\ (J*gx);
     else
    lambda =  plam;
    end
    iter = iter + 1;
      
    vparo = gx + J'*lambda;
    v = [vparo; hx];
    x = 2*(x>0)-1;
    disp(sprintf('%2.0f  %2.12f %2.12f' ,iter, norm(v),alfa) )
end
  
end

