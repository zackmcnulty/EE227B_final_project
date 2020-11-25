function estMMSE = MMSE(H, r, N, modOrd,pskModulator, pskDemodulator, snrLinear)
% Given a real channel matrix H, a SNR rho and y is the received channel
% output a feasible binary x and its associated optimal value

% One shot estimator
	estMMSE = sqrt(snrLinear / N) * H' * ((eye(N) + snrLinear / N * (H * H')) \ r);
	% The above is not boolean and needs to be booleanized
% 	% If booleanizing to {-1,1}
% 	estMMSE = -sign(estMMSE);
	% If booleanizing to {0,1}
% 	estMMSE = sign(estMMSE - 0.5);
% 	estMMSE(estMMSE == -1) = 0;
% 	estMMSE = 1 - estMMSE;
    estMMSE = pskDemodulator(estMMSE);

% MMSE-SIC receiver
%     % Initialization
%     estMMSE = zeros(N*modOrd, 1);
%     orderVec = 1:N;
%     k = N+1;
%     % Start MMSE nulling loop
%     for n = 1:N
%         H = H(:, [1:k-1,k+1:end]);
%         orderVec = orderVec(1, [1:k-1,k+1:end]);
%         % Order algorithm (matrix G calculation) is the only difference
%         % with the ZF-SIC receiver
%         G = (H'*H + ((N-n+1)/snrLinear)*eye(N-n+1)) \ eye(N-n+1);
%         [~, k] = min(diag(G));
%         symNum = orderVec(k);
% 
%         decBits = pskDemodulator(G(k,:) * H' * r);
%         estMMSE(modOrd * (symNum-1) + (1:modOrd)) = decBits;
% 
%         if n < N
%             r = r - H(:, k) * pskModulator(decBits);
%         end
%     end
end