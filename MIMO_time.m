%% parameters
N = 2;                  % Number of transmit antennas
M = 2;                  % Number of receive antennas
EbNoVec = 2:4:14;        % Eb/No in dBtx
modOrd = 1;             % constellation size = 2^modOrd
memory_length = 2;
num_symbol = 5;
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
mlBERCalc = comm.ErrorRate;
% dspBERCalc = comm.ErrorRate;

% Get all bit and symbol combinations for ML receiver
allBits = de2bi(0:2^(modOrd*N)-1, 'left-msb')';
allTxSig = reshape(pskModulator(allBits(:)), N, 2^(modOrd*N));

% Pre-allocate variables to store BER results for speed
[BER_ZF, BER_MMSE, BER_ML, BER_DSP] = deal(zeros(length(EbNoVec), 3));

%% Evaluation
% Set up a figure for visualizing BER results
fig = figure;
grid on;
hold on;
ax = fig.CurrentAxes;
ax.YScale = 'log';
xlim([EbNoVec(1)-0.01 EbNoVec(end)]);
ylim([1e-12 1e0]);
xlabel('Eb/No (dB)');
ylabel('BER');
fig.NumberTitle = 'off';
fig.Renderer = 'zbuffer';
fig.Name = 'Spatial Multiplexing';
title('Uncoded BPSK System');
set(fig,'DefaultLegendAutoUpdate','off');

% Loop over selected EbNo points
for idx = 1:length(EbNoVec)
    % Reset error rate calculation System objects
    reset(zfBERCalc);
    reset(mmseBERCalc);
    reset(mlBERCalc);
%     reset(dspBERCalc);

    % Calculate SNR from EbNo for each independent transmission link
    snrIndB = EbNoVec(idx) + 10*log10(modOrd);
    snrLinear = 10^(0.1*snrIndB);

    while (BER_ZF(idx, 3) < length(EbNoVec)*1e6) && ((BER_MMSE(idx, 2) < 10) || ...
          (BER_ZF(idx, 2) < 10))
%     for loop=1:10000
        % Flat Rayleigh fading channel with independent links
        zeroChan = zeros(M,N);
        rayleighChan = zeros(M, N, memory_length);
        for i = 1:memory_length
            rayleighChan(:,:,i) = (randn(stream, M, N) +  1i*randn(stream, M, N))/sqrt(2);
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
            txSig(:,i) = pskModulator(msg(:,i));
        end
        txSig = reshape(txSig, [numel(txSig), 1]);
        msg = reshape(msg, [numel(msg), 1]);

        % Add noise to faded data
        rxSig = awgn(ChannelMatrix*txSig, snrIndB, 0, stream);

        % Estimation with ZF/ MMSE/ ML
        estZF = zf(ChannelMatrix, rxSig, N*num_symbol, modOrd, pskModulator, pskDemodulator);
        estMMSE = MMSE(ChannelMatrix, rxSig, N*num_symbol, modOrd, pskModulator, pskDemodulator, snrLinear);
%         estML = ML(ChannelMatrix, rxSig, N*num_symbol, modOrd,allTxSig, allBits);
%         estDSP = sdp_relax(rayleighChan, rxSig);

        % Update BER
        BER_ZF(  idx, :) = zfBERCalc(msg, estZF);
        BER_MMSE(idx, :) = mmseBERCalc(msg, estMMSE);
%         BER_ML(  idx, :) = mlBERCalc(msg, estML);
%         BER_DSP(  idx, :) = mlBERCalc(msg, estDSP);
    end

    % Plot results
    semilogy(EbNoVec(1:idx), BER_ZF(  1:idx, 1), 'r*', ...
             EbNoVec(1:idx), BER_MMSE(1:idx, 1), 'bo');
    legend('ZF-SIC', 'MMSE-SIC');
    drawnow;
end

% Draw the lines
semilogy(EbNoVec, BER_ZF(  :, 1), 'r-', ...
         EbNoVec, BER_MMSE(:, 1), 'b-');
hold off;
saveas(fig,'MIMO_time.png');
