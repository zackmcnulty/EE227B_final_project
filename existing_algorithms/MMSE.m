function estMMSE = MMSE(H, r, N, modOrd,pskModulator, pskDemodulator, snrLinear, one_shot)
% Given a real channel matrix H, a SNR rho and y is the received channel
% output a feasible binary x and its associated optimal value
if one_shot
% One shot estimator
    G = (H'*H + (N/snrLinear)*eye(N)) \ eye(N);
    estMMSE = pskDemodulator(G * H' * r);
else
% MMSE-SIC receiver
    % Initialization
    estMMSE = zeros(N*modOrd, 1);
    orderVec = 1:N;
    k = N+1;
    % Start MMSE nulling loop
    for n = 1:N
        H = H(:, [1:k-1,k+1:end]);
        orderVec = orderVec(1, [1:k-1,k+1:end]);
        % Order algorithm (matrix G calculation) is the only difference
        % with the ZF-SIC receiver
        G = (H'*H + ((N-n+1)/snrLinear)*eye(N-n+1)) \ eye(N-n+1);
        [~, k] = min(diag(G));
        symNum = orderVec(k);

        decBits = pskDemodulator(G(k,:) * H' * r);
        estMMSE(modOrd * (symNum-1) + (1:modOrd)) = decBits;

        if n < N
            r = r - H(:, k) * pskModulator(decBits);
        end
    end
end
end