function dec_round = sdp(Hn, rxSig, Nt, pskDemodulator)
M = 100;
chan_op = rxSig;
% SDP relaxation for ML detection
L = [Hn.'*Hn , -Hn.'*chan_op;-chan_op.'*Hn, chan_op.'*chan_op];

cvx_begin sdp  
variable S(Nt+1,Nt+1) symmetric
minimize (trace(L*S))
subject to
    diag(S)== 1;
    S>=0;
cvx_end

if (string(cvx_status) == "Solved")
    % Best rank 1 approximation
    [V,D] = eig(S);
    [maxim,index] = max(diag(D));
    s = V(:,index);
    
    % Rademacher distribution 
    prob = (1+s)/2;
    xls = 2*(rand(M,Nt+1) >= prob') - 1;
    xls = xls*xls(end);

    % find argmin
    argmin = xls(1,:);
    min_val = xls(1,:)*L*xls(1,:)'; 
    for k=2:M
        val = xls(k,:) * L * (xls(k,:).'); 
        if val < min_val
           min_val = val;
           argmin = xls(k,:);
        end
    end

    % feasible x
    x =  argmin(1:end-1)*argmin(end);
    x_round = pskDemodulator(x');

    %Nesterov rounding
    nest = mvnrnd(zeros(Nt,1),S(1:end-1,1:end-1)*S(end,end),1);
    nest_round = pskDemodulator(nest');

    % demapping s 
    dec_round = s(1:end-1)*s(end);
    dec_round = pskDemodulator(dec_round);

end

end