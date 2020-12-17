function [ceq] = bnd(x)

n = length(x);

ceq = zeros(2*n,1);
for i= 1:n
ceq(i) = x(i)+1; 
ceq(n+i) = 1-x(i);
end

end

