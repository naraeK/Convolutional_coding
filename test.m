% Viterbi Decoding Algorithms
clear all; close all; clc;
% n = 2; K = 3; k=1
% L-bit message sequence
%% Transmitter
L = 10^6;    % short packet
depth = L + 2;
%% AWGN channel
EbN0_dB = [0:10];
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
    r = real(r)<0;
    
    m_est = Viterbi_window15(r,L); % test
    
    % counting the errors
    errViterbi(i) = size(find([m - m_est(1:L)]),2);
end

BER_Viterbi = errViterbi / L;

BER_theoretical = 0.5*erfc(sqrt(10.^(EbN0_dB/10))); % theoretical ber uncoded AWGN

figure
semilogy(EbN0_dB,BER_theoretical,'LineWidth',1.5);
hold on
semilogy(EbN0_dB,BER_Viterbi,'LineWidth',1.5);
axis([0 10 10^-5 0.5])
grid on
legend('BER-theoretical,uncoded', 'BER-Viterbi (n=2, K=3, k=1)');
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
function m_est = Viterbi_window15(r,L)
    window = 15;
    depth = L + 2; % n * (L + M -1) / 2
    outer = ceil(depth/window);
    state = [0 0;0 1;1 0;1 1];
    entirePath = zeros(4,depth);
    SM_k = zeros(4,1);  % State Metric: sum of Hamming distance, BM_k
    m_est = zeros(1,depth);
    
    for out = 1:outer
        
        SM_k = zeros(4,1);  % State Metric: sum of Hamming distance, BM_k
        
        if (out==1)||(depth <= 15)
            window = min(depth, window);
            survivorPath = [];

            [survivorPath, SM_k, BM_k] = initial_two_steps(r,survivorPath, SM_k, state);
            entirePath(:,1:2) = survivorPath;
            [survivorPath, SM_k, BM_k] = further_steps(3, window, r, survivorPath, SM_k, state);
            entirePath(:,3:window) = survivorPath;
 %           ls = find(SM_k == min(SM_k));
 %           part_est = traceback(ls, SM_k, survivorPath);
 %           m_est(1:window) = part_est;
 
        elseif out == outer
            survivorPath = zeros(4,depth-window*(out-1));
            [survivorPath, SM_k, BM_k] = further_steps(window*(out-1)+1, depth, r, survivorPath, SM_k, state);
            entirePath(:,window*(out-1)+1:depth) = survivorPath;
 %           ls = 1;
 %           part_est = traceback(ls, SM_k, survivorPath);
 %           m_est(window*(out-1)+1:depth) = part_est;
            
        else 
            survivorPath = zeros(4,window);
            [survivorPath, SM_k, BM_k] = further_steps(window*(out-1)+1, window*out, r, survivorPath, SM_k, state);
            entirePath(:,window*(out-1)+1:window*out) = survivorPath;
 %           ls = find(SM_k == min(SM_k));
 %           part_est = traceback(ls, SM_k, survivorPath);
 %           m_est(window*(out-1)+1:window*out) = part_est;
        end
    end
    ls = 1;
    m_est = traceback(ls, SM_k, entirePath);
end

%
function [survivorPath, SM_k, BM_k] = initial_two_steps(r,survivorPath, SM_k, state)
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
end

%
function [survivorPath, SM_k, BM_k] = further_steps(lef, rig, r, survivorPath, SM_k, state)
    for k = lef:rig
        r_k = r(2*k-1:2*k);
        r_k = [r_k;r_k;r_k;r_k]; % for comparing to 'state'
        BM_k = sum(r_k ~= state,2); % Branch Metric: Hamming Distance 

        SM = SM_k; % SM: k-1 th step
        % state 00
        [SM_k(1), idx] = min([SM(1)+BM_k(1),SM(2)+BM_k(4)]);
        survivorPath(1,k-lef+1) = idx;
        % state 01
        [SM_k(2), idx] = min([SM(3)+BM_k(3),SM(4)+BM_k(2)]);
        survivorPath(2,k-lef+1) = idx+2;
        % state 10
        [SM_k(3), idx] = min([SM(1)+BM_k(4),SM(2)+BM_k(1)]);
        survivorPath(3,k-lef+1) = idx;
        % state 11
        [SM_k(4), idx] = min([SM(3)+BM_k(2),SM(4)+BM_k(3)]);
        survivorPath(4,k-lef+1) = idx+2;    
    end
end

%
function part_est = traceback(ls, SM_k, survivorPath)
    stateTable = [ 0   0   0   0; 0   0   0   0; 1   1   0   0; 0   0   1   1 ]; %%%%
    currState = ls;%find(SM_k == min(SM_k));
    [x y] = size(survivorPath);
    part_est = zeros(1,y);
    for l = y:-1:1
        prevState = survivorPath(currState,l); 
        part_est(l) = stateTable(currState,prevState);
        currState = prevState;
    end
end