% SDP for MIMO ML detection
clear all
close all
clc
Nr = 32; % number of receive antennas
Nt = 32; % number of transmit antennas
fade_var = 1; % fade variance of the channel
% SNR parameters
SNR_dB = 40; % SNR per bit (dB)
noise_var = 1*fade_var*Nt*Nr/(2*10^(0.1*SNR_dB)*Nt);
%------------------------------------------------------
% source
a = randi([0 1],1,Nt);
% bpsk mapper
seq = 1-2*a;
%------ channel----------------------------------------
% fade channel matrix
H = normrnd(0,sqrt(fade_var),Nr,Nt);
% awgn
noise = normrnd(0,sqrt(noise_var),Nr,1);
% channel output
chan_op = H*seq.' + noise;
%------------------------------------------------------
% SDP relaxation for ML detection
L = [H.'*H -H.'*chan_op;-chan_op.'*H chan_op.'*chan_op];
cvx_begin sdp
variable S(Nt+1,Nt+1) symmetric
minimize (trace(L*S))
subject to
diag(S)== 1;
S>=0;
cvx_end
% Best rank 1 approximation
[V,D] = eig(S);
[maxim,index] = max(diag(D));
s = sqrt(maxim)*V(:,index);
s(end)=[];
% demapping s to bits
dec_a = transpose(s<0);
% bit error rate
BER = nnz(a-dec_a)/Nt
