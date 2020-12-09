function [x_opt] = L1_norm_ball(H, y, epsilon, possible_signals)
% L1_norm_ball: Robust maximum likelihood with row-wise L1 norm-ball uncertainty
%   
% param:
%       * H : channel matrix
%       * y : received signal
%       * epsilon : radius of the norm ball (e.g. level of uncertainty)
%       * possible_signals: all signals in constellation map (e.g. {-1,1}^n)

m = size(H,1);
n = size(H, 2); 
y = reshape(y, [n,1]); % force y to be a column vector


A = repmat(y, [1, 2^n]) - H * possible_signals;
objective_values = vecnorm(A, 2) + 2*epsilon * vecnorm(A,1) + epsilon * m;  

% optimize over all possible signals
[opt_val, idx] = min(objective_values);

% find optimal x
x_opt = possible_signals(:, idx);
end

