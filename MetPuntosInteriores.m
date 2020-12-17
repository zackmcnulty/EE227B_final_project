function [ x, mu ] = MetPuntosInteriores( fname, gname, x )
%Método de puntos interiores para
% min f(x)
% s.t. g(x)>=0
%
%

n=length(x);       %Dimensión de la variable primal
gx= feval(gname,x);
p= length(gx);      %Número de restricciones
z=ones(p,1);        %Variable de holgura inicial
mu=ones(p,1);       %Multiplicador de Lagrange inicial
tau=1.0;

gf=gradiente(fname,x);
J=Jacobiana(gname,x);

v=[gf-J'*mu; -gx+z; diag(z)*diag(mu)*ones(p,1)];

maxiter=100;
tol=1.e-05;
damping=0.95;
iter=0;

while(norm(v)>tol && iter<maxiter)
    %RESOLVER EL SISTEMA LINEAL
    B=hessiandelag(fname,gname,x, -mu );
    K=[ B zeros(n,p) -J'; -J eye(p) zeros(p); zeros(p,n) diag(mu) diag(z)];
    ld=- v + [ zeros(n,1); zeros(p,1); tau*ones(p,1)];
   
    d=K\ld;
    dx=d(1:n);
    dz=d(n+1:n+p);
    dmu=d(n+p+1:n+p+p);
    
    
    %---------------------------------------------------------------------
    %RECORTE DEL PASO EN Z Y MU
    %   alfa in (0,1)
    
    alfa=1;
    ind_aff=find(dz<0);
    if (isempty(ind_aff)==0)
        alfa=min(alfa, min(-z(ind_aff)./dz(ind_aff)));
    end
    
    ind_aff=find(dmu<0);
    if (isempty(ind_aff)==0)
        alfa=min(alfa, min(-mu(ind_aff)./dmu(ind_aff)));
    end
    
    alfa=damping*alfa;
    
    
    %---------------------------------------------------------------------
    %ACTUALIZACIÖN
    x=x+alfa*dx;
    z=z+alfa*dz;
    mu=mu+alfa*dmu;
    
    tau=z'*mu/p;
    
    gf=gradiente(fname,x);
    J=Jacobiana(gname,x);
    v=[gf-J'*mu; -gx+z; diag(z)*diag(mu)*ones(p,1)];
    
    x = 2*(x>0)-1;
    iter=iter+1;
    
end

end
