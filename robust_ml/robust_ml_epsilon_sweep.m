% For a fixed SNR, this code tests the Robust L1/L2 norm ball reformulations 
% against the standard Maximum Likelihood algorithm across a range of
% possible uncertainty levels. Namely, as described in the report it sweeps
% across possible choices of epsilon for the norm balls in the uncertainty
% set. 

% Each method is provided the nominal channel matrix 'rayleighChan', but
% signals are produced from the uncertain matrix 'uncertainChannel' that is
% selected randomly within the norm ball uncertain set. We use several
% uncertain channels to try to get an idea of the average performance of
% the algorithms

%% parameters
clear all; close all; clc;

N = 10;                 % Number of transmit antennas
M = 10;                 % Number of receive antennas
EbNo = 10;             % Eb/No in dB
modOrd = 1;            % BPSK modulation (do not change); constellation size = 2^modOrd
ntrials = 1e5;         % number of samples to use to approximate BER


num_matrices = 100;                          % number of uncertain matrices to draw from uncertainty set and test
ntrials_per_matrix = ntrials/num_matrices;   % number of signals to use each uncertain matrix for

p = 2;                      % Lp rowwise error
epsilons = 0.5:0.1:1.25;         %  amount of uncertainty in channel matrix (e.g. rowwise L1 error <= epsilon)

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
[BER_L1, BER_L2, BER_ML] = deal(zeros(length(epsilons), 3));

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
xlim([epsilons(1)-0.01, epsilons(end)+0.5]);
%ylim([1e-12 1e-1]);
xlabel('\epsilon');
ylabel('BER');
fig.NumberTitle = 'off';
fig.Renderer = 'zbuffer';
fig.Name = 'Spatial Multiplexing';
title('Uncoded BPSK System');
set(fig,'DefaultLegendAutoUpdate','off');


% Calculate SNR from EbNo for each independent transmission link
snrIndB = EbNo + 10*log10(modOrd);
snrLinear = 10^(0.1*snrIndB);


% Loop over selected epsilons points
for idx = 1:length(epsilons)
    % Reset error rate calculation System objects
    reset(L1_BERCalc);
    reset(L2_BERCalc);
    reset(mlBERCalc);
    
    epsilon = epsilons(idx);
   
    for jj=1:num_matrices
        
        % add uncertainty to channel
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
    semilogy(epsilons(1:idx), BER_L1(  1:idx, 1), 'r*', ...
             epsilons(1:idx), BER_L2(1:idx, 1), 'bo', ...
             epsilons(1:idx), BER_ML(  1:idx, 1), 'gs');
    legend('L1', 'L2', 'ML');
    drawnow;
end

% Draw the lines
semilogy(epsilons, BER_L1(  :, 1), 'r-', ...
         epsilons, BER_L2(:, 1), 'b-', ...
         epsilons, BER_ML(  :, 1), 'g-');
hold off;


