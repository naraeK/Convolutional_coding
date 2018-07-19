% Fading channel
% Soft decision with Maximum likelihood (ML) decision
% Viterbi Decoding Algorithms
clear all;  clc;
% n = 2; K = 3; k=1
% L-bit message sequence
%% Transmitter
L = 10^6;    % short packet

tic;
%% AWGN channel & Rayleigh fading channel
EbN0_dB = [0:1:15];
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
    r = s + 10^(-EN0_dB(i)/20)*n; % additive white gaussian noise
 
    % Rayleigh channel fading
    h = sigma/sqrt(2)*[randn(1,length(c)) + j*randn(1,length(c))];  % Assume - constant during the transmission
%     h = 1/sqrt(2)*[randn(1) + j*randn(1)];
    % Send over Gaussian Link to the receiver
    r_fading = h.*s + 10^(-EN0_dB(i)/20)*n; % additive white gaussian noise

    % Equalization to remove fading effects. Ideal Equalization Considered
    r_fading = r_fading./h;%real(r_fading)./real(h);
    r_fading = real(r_fading);
    
    % BPSK demodulator at the Receiver
    r_fading_hard = real(r_fading) < 0;
    r_hard = real(r)<0;
    
    % Decoding
    m_est1 = ViterbiDecoder(r_hard,'hard'); 
    m_est2 = ViterbiDecoder(r_fading_hard,'hard');
    m_est3 = ViterbiDecoder(r,'soft');
    m_est4 = ViterbiDecoder(r_fading,'soft');
    
    % counting the errors
    errViterbi_h(i) = size(find([m - m_est1(1:L)]),2);
    errViterbi_fading(i) = size(find([m - m_est2(1:L)]),2);
    errViterbi_soft(i) = size(find([m - m_est3(1:L)]),2);
    errViterbi_fading_soft(i) = size(find([m - m_est4(1:L)]),2);
end

%% BER graphs
BER_Viterbi_h = errViterbi_h / L;
BER_Viterbi_fading = errViterbi_fading / L;
BER_Viterbi_s = errViterbi_soft / L;
BER_Viterbi_f_s = errViterbi_fading_soft / L;

% Theoretical BER in AWGN
BER_theoretical = 0.5*erfc(sqrt(10.^(EbN0_dB/10))); % theoretical ber uncoded AWGN
% Rayleigh Tehoretical BER
SNR = 10.^(EbN0_dB/10);
BER_theor_fading = 0.5.*(1-sqrt(SNR./(SNR+1)));

figure
semilogy(EbN0_dB,BER_theoretical,'LineWidth',1.5);
hold on
semilogy(EbN0_dB,BER_theor_fading,'LineWidth',1.5);
hold on
semilogy(EbN0_dB,BER_Viterbi_fading,'LineWidth',1.5);
hold on
semilogy(EbN0_dB,BER_Viterbi_f_s,'LineWidth',1.5);
hold on
semilogy(EbN0_dB,BER_Viterbi_h,'LineWidth',1.5);
hold on
semilogy(EbN0_dB,BER_Viterbi_s,'LineWidth',1.5);
axis([0 15 10^-6 0.5])
grid on
legend('Theoretical,uncoded,AWGN', 'Theoretical,Rayleigh fading', 'Viterbi(hard),Rayleigh fading', 'Viterbi(soft),Rayleigh fading', 'Viterbi(hard),AWGN', 'Viterbi(soft),AWGN');
xlabel('Eb/No, dB');
ylabel('Bit Error Rate');
title('BER with Viterbi decoding for BPSK');

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

%% Viterbi Decoder - hard decision & soft decision (Maximum Likelihood decisions)
function m_est = ViterbiDecoder(r,mode)
    survivorPath = zeros(4,length(r)/2);
    SM_k = zeros(4,1);  % State Metric: sum of Hamming distance, BM_k

    survivorPath(:,1)=[1;1;1;1];%[1;0;1;0];
    r_12 = real(r(1:4));
    r_12 = [r_12;r_12;r_12;r_12]; % for comparing to 'state'

    if mode == 'hard'
        state = [0 0;0 1;1 0;1 1];
        BM_k = sum(r_12 ~= [[0 0;1 1;1 1;0 0], state],2); % branch metric
    elseif mode == 'soft'
        state = [1 1;1 -1;-1 1;-1 -1];
        BM_k = abs(r_12-[[1 1;-1 -1;-1 -1;1 1],state])%.^2; % branch metric
        BM_k = sum(BM_k,2);
    end
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
        r_k = real(r(2*k-1:2*k));
        r_k = [r_k;r_k;r_k;r_k]; % for comparing to 'state'

        if mode == 'hard'
        BM_k = sum(r_k ~= state,2); % Branch Metric: Hamming Distance 
        elseif mode == 'soft'           
        BM_k = sum(abs(r_k - state),2); % Branch Metric: Hamming Distance 
        end
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

