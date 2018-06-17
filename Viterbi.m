% Viterbi Decoding Algorithms
clear all; clc;
% n = 2; K = 3; k=1
% L-bit message sequence
%% Transmitter
L = 4; %256;    % short packet

%% AWGN channel
EbN0_dB = [0:10];
EN0_dB = EbN0_dB - 10*log10(2); % (1-bit -> 2-bit); E/N0 = 1/2 * Eb/N0
sigma = 1; % Gaussian noise variance
   
for i = 1:length(EN0_dB)
    % Transmitter
    [m,c] = ConvolutionalEncoder(L); % message m, codeword c
    % BPSK modulation
    s = 1 - 2*c; % 0 -> 1; 1 -> -1 % *sqrt(Eb)
    
    n = sigma/sqrt(2)*[randn(size(c)) + j*randn(size(c))]; % Gaussian noise, 0 mean and sigma^2 variance 
    % Noise addition
    r = s + 10^(-EN0_dB(i)/20)*n; % additive white gaussian noise

    % receiver - hard decision decoding
    r = real(r)>0;
    
    m_est = ViterbiDecoder(c); % test
end
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
    SM_k = zeros(4,1);  % State Metric
    
    survivorPath(:,1)=[1;0;1;0];
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
        BM_k = sum(r_k ~= state,2); % Hamming Distance: Branch Metric

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
    stateTable = [ 0   0   0   0; 0   0   0   0; 1   1   0   0; 0   0   1   1 ]; %%%%
    currState = 1;%find(SM_k == min(SM_k));
    m_est = zeros(1,length(r)/2);
    for l = length(r)/2:-1:1
        prevState = survivorPath(currState,l); 
        m_est(l) = stateTable(currState,prevState);
        currState = prevState;
    end
    m_est = m_est(1:length(r)/2 -2);
end
