function [x] = sdp_relax(H, y, pskDemodulator)
% Given a real channel matrix H and y is the received channel
% output a feasible binary x (-1,1) using the SDP relaxation of the ML problem 

% Solving the SDP relaxation

n = size(H,1);
y = reshape(y, [n,1]);

A = H.' * H;
B = -(H.') * y;
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

%%%%% Generating a feasible (rank 1) x from SDP relaxation solution %%%%%

% 1) compute largest eigenvector/eigenvalue
[v, ~] = eigs(X, 1);

% 2) sample iid binary vector samples

L = 1000; % number samples
prob = (1+v) / 2;

% create L samples where each sample x has p(x(i) = 1) = (1+v_i)/2
xls = 2*(rand(L, n+1) <= reshape(prob, [1,n+1])) - 1; 


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
x = pskDemodulator(argmin(end) * argmin(1:end-1).');
end

