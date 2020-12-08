function estML = ML(H, r, N, modOrd,allTxSig, allBits)
% Given a real channel matrix H, a SNR rho and y is the received channel
% output a feasible binary x and its associated optimal value

  % ML receiver
        [~, k] = min(sum(abs(repmat(r,[1,2^(modOrd*N)]) - H*allTxSig).^2));
        estML = allBits(:,k);
end