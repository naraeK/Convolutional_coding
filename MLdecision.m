% Maximum likelihood (ML) decision
% Viterbi Decoding Algorithms
clear all; close all; clc;
% n = 2; K = 3; k=1
% L-bit message sequence
%% Transmitter
L = 10^5;    % short packet
tic;
%% AWGN channel
EbN0_dB = [0:1:12];
EN0_dB = EbN0_dB - 10*log10(2); % (1-bit -> 2-bit); E/N0 = 1/2 * Eb/N0
sigma = 1; % Gaussian noise variance
errViterbi = [];   
for i = 1:length(EN0_dB)
    % Transmitter
    [m,c] = ConvolutionalEncoder(L); % message m, codeword c
    % BPSK modulation
    s = 1 - 2*c; % 0 -> 1; 1 -> -1 % *sqrt(Eb)
    
    n = sigma/sqrt(2)*[randn(size(c)) + j*randn(size(c))]; % Gaussian noise, 0 mean and sigma^2 variance 
    % Noise addition
    r = s + 10^(-EN0_dB(i)/20)*n; % additive white gaussian noise

    % receiver - hard decision decoding
    r_hard = real(r)<0;
    
    m_est1 = ViterbiDecoder(r_hard);
    m_est2 = Viterbi_ML(r);
    
    % counting the errors
    errViterbi_h(i) = size(find([m - m_est1(1:L)]),2);
    errViterbi_ML(i) = size(find([m - m_est2(1:L)]),2);
end

%% BER graphs
BER_Viterbi_h = errViterbi_h / L;
BER_Viterbi_ML = errViterbi_ML / L;

BER_theoretical = 0.5*erfc(sqrt(10.^(EbN0_dB/10))); % theoretical ber uncoded AWGN

figure
semilogy(EbN0_dB,BER_theoretical,'LineWidth',1.5);
hold on
semilogy(EbN0_dB,BER_Viterbi_h,'LineWidth',1.5);
hold on
semilogy(EbN0_dB,BER_Viterbi_ML,'LineWidth',1.5);
axis([0 12 10^-7 0.5])
grid on
legend('BER-theoretical,uncoded', 'BER-Viterbi (hard)', 'BER-Viterbi (ML)');
xlabel('Eb/No, dB');
ylabel('Bit Error Rate');
title('BER with Viterbi decoding for BPSK in AWGN');

toc;
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

%% Viterbi Decoder - hard decision
function m_est = ViterbiDecoder(r)
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

%% Viterbi Decoder - ML decision
function m_est = Viterbi_ML(r)
    state = [1 1;1 -1;-1 1;-1 -1];
    survivorPath = zeros(4,length(r)/2);
    SM_k = zeros(4,1);  % State Metric: sum of Hamming distance, BM_k
    
    survivorPath(:,1)=[1;0;1;0];
    r_12 = real(r(1:4));
    r_12 = [r_12;r_12;r_12;r_12]; % for comparing to 'state'

    BM_k = (r_12-[[1 1;-1 -1;-1 -1;1 1],state]).^2; % branch metric
    BM_k = sum(BM_k,2);
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

        BM_k = sum((r_k - state).^2,2); % Branch Metric: Hamming Distance 

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
