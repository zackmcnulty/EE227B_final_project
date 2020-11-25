%% parameters
N = 2;                  % Number of transmit antennas
M = 2;                  % Number of receive antennas
EbNoVec = 2:4:14;        % Eb/No in dB
modOrd = 1;             % constellation size = 2^modOrd
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
ylim([1e-5 1e-1]);
xlabel('SNR (dB)');
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

    while (BER_ZF(idx, 3) < 1e5) && ((BER_MMSE(idx, 2) < 100) || ...
          (BER_ZF(idx, 2) < 100) ||  (BER_ML(idx, 2)   < 100))
        % Create random bit vector to modulate
        msg = randi(stream, [0 1], [N*modOrd, 1]);

        % Modulate data
        txSig = pskModulator(msg);

        % Flat Rayleigh fading channel with independent links
        rayleighChan = (randn(stream, M, N) +  1i*randn(stream, M, N))/sqrt(2);

        % Add noise to faded data
        rxSig = awgn(rayleighChan*txSig, snrIndB(idx), 0, stream);

        estZF = zf(rayleighChan, rxSig, N, modOrd, pskModulator, pskDemodulator);
        estMMSE = MMSE(rayleighChan, rxSig, N, modOrd, pskModulator, pskDemodulator, snrLinear(idx));
        estML = ML(rayleighChan, rxSig, N, modOrd,allTxSig, allBits);

        % Update BER
        BER_ZF(  idx, :) = zfBERCalc(msg, estZF);
        BER_MMSE(idx, :) = mmseBERCalc(msg, estMMSE);
        BER_ML(  idx, :) = mlBERCalc(msg, estML);
    end

    % Plot results
    semilogy(snrIndB(1:idx), BER_ZF(  1:idx, 1), 'r*', ...
             snrIndB(1:idx), BER_MMSE(1:idx, 1), 'bo', ...
             snrIndB(1:idx), BER_ML(  1:idx, 1), 'gs');
    legend('ZF-SIC', 'MMSE-SIC', 'ML');
    drawnow;
end
% Draw the lines
semilogy(snrIndB, BER_ZF(  :, 1), 'r-', ...
         snrIndB, BER_MMSE(:, 1), 'b-', ...
         snrIndB, BER_ML(  :, 1), 'g-');
hold off;
saveas(fig,'MIMO.png');
