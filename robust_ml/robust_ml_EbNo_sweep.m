%% parameters
clear all; close all; clc;

N = 2;                 % Number of transmit antennas
M = 2;                 % Number of receive antennas
EbNoVec = 2:2:16;      % Eb/No in dB
modOrd = 1;            % BPSK modulation (do not change); constellation size = 2^modOrd
ntrials = 1e5;         % number of samples to use to approximate BER


num_matrices = 1e2;
ntrials_per_matrix = ntrials/num_matrices;

epsilon = 5;         %  amount of uncertainty in channel matrix (e.g. rowwise L1 error <= epsilon)

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
L1_BERCalc = comm.ErrorRate;
L2_BERCalc = comm.ErrorRate;
mlBERCalc = comm.ErrorRate;

% Get all bit and symbol combinations for ML receiver
allBits = de2bi(0:2^(modOrd*N)-1, 'left-msb')';
allTxSig = reshape(pskModulator(allBits(:)), N, 2^(modOrd*N));

% Pre-allocate variables to store BER results for speed
[BER_L1, BER_L2, BER_ML] = deal(zeros(length(EbNoVec), 3));

 % Flat Rayleigh fading channel with independent links
% rayleighChan = (randn(stream, M, N) +  1i*randn(stream, M, N))/sqrt(2);
rayleighChan = randn(stream, M, N); % nominal channel matrix

%% Evaluation
% Set up a figure for visualizing BER results
fig = figure;
grid on;
hold on;
ax = fig.CurrentAxes;
ax.YScale = 'log';
xlim([EbNoVec(1)-0.01 EbNoVec(end)]);
ylim([1e-12 1e-1]);
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
    reset(L1_BERCalc);
    reset(L2_BERCalc);
    reset(mlBERCalc);

    % Calculate SNR from EbNo for each independent transmission link
    snrIndB = EbNoVec(idx) + 10*log10(modOrd);
    snrLinear = 10^(0.1*snrIndB);
   
    for jj=1:num_matrices
        
        % add uncertainty to channel
        p = 2;                                   % Lp rowwise error
        U = rand(stream, M,N) - 0.5;
        U = (epsilon * U ./ vecnorm(U,p,2) ) .* rand(stream, M,1);
        uncertainChannel = rayleighChan + U;     % true channel = nominal + uncertainty


        for j=1:ntrials_per_matrix


            % Create random bit vector to modulate
            msg = randi(stream, [0 1], [N*modOrd, 1]);

            % Modulate data
            txSig = pskModulator(msg);


            % Add noise to faded data
            rxSig = real(awgn(uncertainChannel*txSig, snrIndB, 'measured', stream));

            % Estimation with ML, L1-robust ML, L2-robust ML
            estL1 = L1_norm_ball(rayleighChan, rxSig, epsilon, allTxSig, allBits);
            estL2 = L2_norm_ball(rayleighChan, rxSig, epsilon, allTxSig, allBits);
            estML = ML(rayleighChan, rxSig, N, modOrd,allTxSig, allBits);

            % Update BER
            BER_L1(  idx, :) = L1_BERCalc(msg, estL1);
            BER_L2(idx, :) = L2_BERCalc(msg, estL2);
            BER_ML(  idx, :) = mlBERCalc(msg, estML);

        end
    end
    % Plot results
    semilogy(EbNoVec(1:idx), BER_L1(  1:idx, 1), 'r*', ...
             EbNoVec(1:idx), BER_L2(1:idx, 1), 'bo', ...
             EbNoVec(1:idx), BER_ML(  1:idx, 1), 'gs');
    legend('L1', 'L2', 'ML');
    drawnow;
end

% Draw the lines
semilogy(EbNoVec, BER_L1(  :, 1), 'r-', ...
         EbNoVec, BER_L2(:, 1), 'b-', ...
         EbNoVec, BER_ML(  :, 1), 'g-');
hold off;


