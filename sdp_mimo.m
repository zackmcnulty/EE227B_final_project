% SDP for MIMO ML detection
clear all
close all
clc
Nr = 32; % number of receive antennas
Nt = 32; % number of transmit antennas
fade_var = 1; % fade variance of the channel
rep = 1; % number of replications
M = 30; % number samples


%------------------------------------------------------
BER_round = zeros(1, 20); 
BER_rad = zeros(1, 20); 
BER_nor = zeros(1, 20); 

flag = 1;
for i = 1:rep
    j=1;

    % source
    a = randi([0 1],1,Nt);
    % bpsk mapper
    seq = 1-2*a;
    %------ channel----------------------------------------
    % fade channel matrix
    H = normrnd(0,sqrt(fade_var),Nr,Nt);
    % awgn
    noise = normrnd(0,1,Nr,1);

    for SNR_dB = 2:1:20
        % SNR parameters
        noise_var = (10^(0.1*SNR_dB)*Nt)/(1*Nt*Nr);
        Hn = sqrt(noise_var).*H;

        % channel output
        chan_op = Hn*seq.' + noise;
        %------------------------------------------------------
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

            % demapping s 
            dec_round = 2*(s>0)-1;
            dec_round(1:end) = dec_round*dec_round(end); dec_round(end) = [];

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
            x =  argmin(1:end-1);

            %Nesterov rounding
            nest = mvnrnd(zeros(Nt,1),S(1:end-1,1:end-1),1);
            nest_round = 2*(nest>0)-1;

            % bit error rate 
            BER_round(j) = BER_round(j) + nnz(seq-dec_round')/Nt;
            BER_rad(j) = BER_rad(j) + nnz(seq-x)/Nt;
            BER_nor(j) = BER_nor(j) + nnz(seq-nest_round)/Nt;

        end
        j = j+1
    end
end



