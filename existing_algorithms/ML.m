function estML = ML(H, r, N, modOrd,allTxSig, allBits)

  % ML receiver
        [~, k] = min(sum(abs(repmat(r,[1,2^(modOrd*N)]) - H*allTxSig).^2));
        estML = allBits(:,k);
end