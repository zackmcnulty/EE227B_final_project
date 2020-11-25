function estZF = zf(H, r, N, modOrd,pskModulator, pskDemodulator)
% Given a real channel matrix H, a SNR rho and y is the received channel
% output a feasible binary x and its associated optimal value

% ZF-SIC receiver
       
    % Initialization
    estZF = zeros(N*modOrd, 1);
    orderVec = 1:N;
    k = N+1;
    % Start ZF nulling loop
    for n = 1:N
        % Shrink H to remove the effect of the last decoded symbol
        H = H(:, [1:k-1,k+1:end]);
        % Shrink order vector correspondingly
        orderVec = orderVec(1, [1:k-1,k+1:end]);
        % Select the next symbol to be decoded
        G = (H'*H) \ eye(N-n+1); % Same as inv(H'*H), but faster
        [~, k] = min(diag(G));
        symNum = orderVec(k);

        % Hard decode the selected symbol
        decBits = pskDemodulator(G(k,:) * H' * r);
        estZF(modOrd * (symNum-1) + (1:modOrd)) = decBits;

%         % Subtract the effect of the last decoded symbol from r
%         if n < N
%             r = r - H(:, k) * pskModulator(decBits);
%         end
    end
end