% Fading channel
% Viterbi Decoding Algorithms
clear all; close all; clc;
% n = 2; K = 3; k=1
% L-bit message sequence
%% Transmitter
L = 10^6;    % short packet

%% Fading channel/ Rayleigh
EbN0_dB = [0:2:12];
EN0_dB = EbN0_dB - 10*log10(2); % (1-bit -> 2-bit); E/N0 = 1/2 * Eb/N0
sigma = 1; % Gaussian noise variance
errViterbi = [];   
for i = 1:length(EN0_dB)
    % Transmitter
    [m,c] = ConvolutionalEncoder(L); % message m, codeword c
    
    % BPSK modulation
    s = 1 - 2*c; % 0 -> 1; 1 -> -1 % *sqrt(Eb)
    
    % Noise addition
    n = sigma/sqrt(2)*[randn(size(c)) + j*randn(size(c))]; % Gaussian noise, 0 mean and sigma^2 variance 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Rayleigh channel fading
    h = sigma/sqrt(2)*[randn(1,length(c)) + j*randn(1,length(c))];  % Assume - constant during the transmission
    
    % Send over Gaussian Link to the receiver
    r = s + 10^(-EN0_dB(i)/20)*n; % additive white gaussian noise 
    r_fading = h.*s + 10^(-EN0_dB(i)/20)*n; % additive white gaussian noise

    % Equalization to remove fading effects. Ideal Equalization Considered
    r_fading = r_fading./h;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % BPSK demodulator at the Receiver
    r_fading = real(r_fading) < 0;
    r_hard = real(r)<0;
    
    % Decoding
    m_est1 = ViterbiDecoder(r_hard); 
    m_est2 = ViterbiDecoder(r_fading);
    
    % counting the errors
    errViterbi(i) = size(find([m - m_est1(1:L)]),2);
    errViterbi_fading(i) = size(find([m - m_est2(1:L)]),2);
end

%% BER graphs
BER_Viterbi = errViterbi / L;
BER_Viterbi_fading = errViterbi_fading / L;

BER_theoretical = 0.5*erfc(sqrt(10.^(EbN0_dB/10))); % theoretical ber uncoded AWGN

figure
semilogy(EbN0_dB,BER_theoretical,'LineWidth',1.5);
hold on
semilogy(EbN0_dB,BER_Viterbi,'LineWidth',1.5);
hold on
semilogy(EbN0_dB,BER_Viterbi_fading,'LineWidth',1.5);
axis([0 12 10^-7 0.5])
grid on
legend('BER-theoretical,uncoded', 'BER-Viterbi (n=2, K=3, k=1)','BER-Viterbi, fading');
xlabel('Eb/No, dB');
ylabel('Bit Error Rate');
title('BER with Viterbi decoding for BPSK in AWGN');


%% Convolutional Coding Encoder
% n = 2; K = 3; k=1
% L-bit message sequence
function [m,c] = ConvolutionalEncoder(L)
    m = rand(1,L) > 0.5;
    g1 = [1 1 1];
    g2 = [1 0 1];
    upper = mod(conv(g1,m),2);
    lower = mod(conv(g2,m),2);
    c = [upper; lower]; 
    c = c(:).';     % codeword
end

%% Viterbi Decoder
function m_est = ViterbiDecoder(r)%window_size
    state = [0 0;0 1;1 0;1 1];
    survivorPath = zeros(4,length(r)/2);
    SM_k = zeros(4,1);  % State Metric: sum of Hamming distance, BM_k
    
    survivorPath(:,1)=[1;1;1;1];%[1;0;1;0];
    r_12 = r(1:4);
    r_12 = [r_12;r_12;r_12;r_12]; % for comparing to 'state'
    BM_k = sum(r_12 ~= [[0 0;1 1;1 1;0 0], state],2); % branch metric
    % state 00
    SM_k(1) = BM_k(1);
    survivorPath(1,2) = 1;
    % state 01
    SM_k(2) = BM_k(3);
    survivorPath(2,2) = 3;
    % state 10
    SM_k(3) = BM_k(4);
    survivorPath(3,2) = 1;
    % state 11
    SM_k(4) = BM_k(2);
    survivorPath(4,2) = 3;
    
    for k = 3:length(r)/2
        r_k = r(2*k-1:2*k);
        r_k = [r_k;r_k;r_k;r_k]; % for comparing to 'state'
        BM_k = sum(r_k ~= state,2); % Branch Metric: Hamming Distance 

        SM = SM_k; % SM: k-1 th step
        % state 00
        [SM_k(1), idx] = min([SM(1)+BM_k(1),SM(2)+BM_k(4)]);
        survivorPath(1,k) = idx;
        % state 01
        [SM_k(2), idx] = min([SM(3)+BM_k(3),SM(4)+BM_k(2)]);
        survivorPath(2,k) = idx+2;
        % state 10
        [SM_k(3), idx] = min([SM(1)+BM_k(4),SM(2)+BM_k(1)]);
        survivorPath(3,k) = idx;
        % state 11
        [SM_k(4), idx] = min([SM(3)+BM_k(2),SM(4)+BM_k(3)]);
        survivorPath(4,k) = idx+2;    
    end

    % Traceback
    stateTable = [ 0   0   0   0; 0   0   0   0; 1   1   0   0; 0   0   1   1 ]; 
    currState = 1;% find(SM_k == min(SM_k));
    m_est = zeros(1,length(r)/2);
    for l = length(r)/2:-1:1
        prevState = survivorPath(currState,l); 
        m_est(l) = stateTable(currState,prevState);
        currState = prevState;
    end
end
