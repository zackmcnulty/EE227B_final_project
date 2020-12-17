% This code is modified base on https://www.mathworks.com/help/comm/ug/spatial-multiplexing.html
clear;
%% parameters
N =10;                  % Number of transmit antennas
M = 10;                  % Number of receive antennas
EbNoVec = 3:3:15;       % Eb/No in dB
modOrd = 1;             % constellation size = 2^modOrd
memory_length = 2;
num_symbol = 2;
num_error = 20;
%% setup simulation
% Create a local random stream to be used by random number generators for
% repeatability.
stream = RandStream('mt19937ar');

% Create PSK modulator and demodulator System objects
pskModulator   = comm.PSKModulator(...
            'ModulationOrder',  2^modOrd, ...
            'PhaseOffset',      0, ...
            'BitInput',         true);
pskDemodulator = comm.PSKDemodulator( ...
            'ModulationOrder',  2^modOrd, ...
            'PhaseOffset',      0, ...
            'BitOutput',        true);

% Create error rate calculation System objects for 3 different receivers
zfBERCalc = comm.ErrorRate;
mmseBERCalc = comm.ErrorRate;
zfsicBERCalc = comm.ErrorRate;
mmsesicBERCalc = comm.ErrorRate;
mlBERCalc = comm.ErrorRate;
sdpBERCalc = comm.ErrorRate;

% Get all bit  and symbol combinations for ML receiver
allBits = de2bi(0:2^(modOrd*N*num_symbol)-1, 'left-msb')';
allTxSig = reshape(pskModulator(allBits(:)), N*num_symbol, 2^(modOrd*N*num_symbol));

% Pre-allocate variables to store BER results for speed
[BER_ZF, BER_MMSE, BER_ZF_SIC, BER_MMSE_SIC, BER_ML, BER_SDP] = deal(zeros(length(EbNoVec), 3));

%% Evaluation
% Calculate SNR from EbNo for each independent transmission link
snrIndB = EbNoVec + 10*log10(modOrd);
snrLinear = 10.^(0.1*snrIndB);
% Set up a figure for visualizing BER results
fig = figure;
grid on;
hold on;
ax = fig.CurrentAxes;
ax.YScale = 'log';
xlim([snrIndB(1)-0.01 snrIndB(end)]);
ylim([1e-6 1e0]);
xlabel('SNR (dB)');
ylabel('BER');
fig.NumberTitle = 'off';
fig.Renderer = 'zbuffer';
fig.Name = 'Spatial Multiplexing';
title('BPSK System');
set(fig,'DefaultLegendAutoUpdate','off');

% Loop over selected EbNo points
for idx = 1:length(EbNoVec)
    % Reset error rate calculation System objects
    reset(zfBERCalc);
    reset(mmseBERCalc);
    reset(zfsicBERCalc);
    reset(mmsesicBERCalc);
    reset(mlBERCalc);
    reset(sdpBERCalc);

    while (BER_ZF(idx, 3) < 5e5) && ((BER_MMSE(idx, 2) < num_error) || ...
          (BER_ZF_SIC(idx, 2) < num_error) && (BER_MMSE_SIC(idx, 2) < num_error) || ...
          (BER_ZF(idx, 2) < num_error) ||  (BER_ML(idx, 2)   < num_error)||  (BER_SDP(idx, 2) < num_error))
        % Flat Rayleigh fading channel with independent links
        zeroChan = zeros(M,N);
        rayleighChan = zeros(M, N, memory_length);
        for i = 1:memory_length
            rayleighChan(:,:,i) = randn(stream, M, N);
%             rayleighChan(:,:,i) = (randn(stream, M, N) +  1i*randn(stream, M, N))/sqrt(2);
        end
        ChannelMatrix = zeros((num_symbol+memory_length-1)*M, num_symbol*N);
        for j = 1:num_symbol+memory_length-1
            for i = 1:num_symbol
                if (j-i>=0) && (j-i<memory_length)
                    ChannelMatrix((j-1)*M+1:j*M, (i-1)*N+1:i*N)=rayleighChan(:,:,j-i+1);
                else
                    ChannelMatrix((j-1)*M+1:j*M, (i-1)*N+1:i*N)=zeroChan;
                end
            end
        end
        txSig = zeros(N,num_symbol);
        msg = zeros(N,num_symbol);
        for i = 1:num_symbol
            % Create random bit vector to modulate
            msg(:,i) = randi(stream, [0 1], [N*modOrd, 1]);

            % Modulate data
            txSig(:,i) = real(pskModulator(msg(:,i)));
        end
        txSig = reshape(txSig, [numel(txSig), 1]);
        msg = reshape(msg, [numel(msg), 1]);

        % Add noise to faded data
        rxSig = awgn(ChannelMatrix*txSig, snrIndB(idx),'measured', stream);

        estZF = zf(ChannelMatrix, rxSig, N*num_symbol, modOrd, pskModulator, pskDemodulator, true);
        estMMSE = MMSE(ChannelMatrix, rxSig, N*num_symbol, modOrd, pskModulator, pskDemodulator, snrLinear(idx), true);
        estZF_SIC = zf(ChannelMatrix, rxSig, N*num_symbol, modOrd, pskModulator, pskDemodulator, false);
        estMMSE_SIC = MMSE(ChannelMatrix, rxSig, N*num_symbol, modOrd, pskModulator, pskDemodulator, snrLinear(idx), false);
%         [Chan_real, rxSig_real] = complex2real(rayleighChan, rxSig);
%         estMMSE_c = MMSE(Chan_real, rxSig_real, 2*N, modOrd, pskModulator, pskDemodulator, snrLinear(idx));
%         estMMSE = estMMSE_c(1:end/2);
        estML = ML(ChannelMatrix, rxSig, N*num_symbol, modOrd,allTxSig, allBits);
        estSDP = sdp(ChannelMatrix, rxSig, N*num_symbol, pskDemodulator);
%         estSDP = estSDP(1:end/2);

        % Update BER
        BER_ZF(  idx, :) = zfBERCalc(msg, estZF);
        BER_MMSE(idx, :) = mmseBERCalc(msg, estMMSE);
        BER_ZF_SIC(  idx, :) = zfsicBERCalc(msg, estZF_SIC);
        BER_MMSE_SIC(idx, :) = mmsesicBERCalc(msg, estMMSE_SIC);
        BER_ML(  idx, :) = mlBERCalc(msg, estML);
        BER_SDP(  idx, :) = sdpBERCalc(msg, estSDP);
    end

    % Plot results
    semilogy(snrIndB(1:idx), BER_ZF(  1:idx, 1), 'r*', ...
             snrIndB(1:idx), BER_MMSE(1:idx, 1), 'bo', ...
             snrIndB(1:idx), BER_ZF_SIC(  1:idx, 1), 'r*', ...
             snrIndB(1:idx), BER_MMSE_SIC(1:idx, 1), 'bo', ...
             snrIndB(1:idx), BER_ML(  1:idx, 1), 'gs', ...
             snrIndB(1:idx), BER_SDP(  1:idx, 1), 'kx');
    legend('ZF', 'MMSE', 'ZF-SIC', 'MMSE-SIC', 'ML', 'SDP');
    drawnow;
end
% Draw the lines
semilogy(snrIndB, BER_ZF(  :, 1), 'r-', ...
         snrIndB, BER_MMSE(:, 1), 'b-', ...
         snrIndB, BER_ZF_SIC(  :, 1), 'r--', ...
         snrIndB, BER_MMSE_SIC(:, 1), 'b--', ...
         snrIndB, BER_ML(  :, 1), 'g-', ...
         snrIndB, BER_SDP(  :, 1), 'k-');
hold off;
saveas(fig,'MIMO_mem.png');
