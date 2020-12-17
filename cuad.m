function [ceq] = cuad(x)

ceq = zeros(length(x),1);
for i= 1:length(x)
ceq(i) = x(i)^2-1; % Compute nonlinear equalities at x.
end

end

