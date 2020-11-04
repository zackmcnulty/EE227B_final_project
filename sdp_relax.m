function [x, min_val] = sdp_relax(rho, H, y)
% Given a real channel matrix H, a SNR rho and y is the received channel
% output a feasible binary x and its associated optimal value

% Solving the SDP relaxation

n = size(H,1);

y = reshape(y, [n,1]);

A = rho/n * (H.')*H;
B = -sqrt(rho/n) * (H.') * y;
C = (y.')*y;

Q = [A, B; 
     B.', C];
 
cvx_begin SDP quiet
   variable X(n+1,n+1) symmetric
   
   minimize(trace(Q*X))
   subject to
        X >= 0
        diag(X) == 1  
cvx_end

% Generating a feasible (rank 1) x from SDP relaxation solution

% 1) compute largest eigenvector/eigenvalue
[V, D] = eig(X);
v = V(:, end);
lambda = D(end, end);


% 2) sample iid binary vector samples

L = 30; % number samples
prob = (1+v) / 2;
xls = 2*(rand(L, n+1) >= prob.' ) - 1;


% 3) find argmin
argmin = xls(1,:);
min_val = xls(1,:) * Q * (xls(1,:).'); 
for i=2:L
    
    val = xls(i,:) * Q * (xls(i,:).'); 
    if val < min_val
       min_val = val;
       argmin = xls(i,:);
    end
end

% feasible x
x = argmin(end) * argmin(1:end-1).';
end

