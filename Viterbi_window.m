% Viterbi Decoding Algorithms
clear all; close all; clc;
% n = 2; K = 3; k=1
% L-bit message sequence
%% The size of the sample or the packet length
L = 10^5;    % valuable for short packet

%% Binary Symmetric Channel (BSC) channel
% EbN0_dB = [0:10];
% EN0_dB = EbN0_dB - 10*log10(2); % (1-bit -> 2-bit); E/N0 = 1/2 * Eb/N0
% errViterbi = [];   
% for i = 1:length(EN0_dB)
%     % Transmitter
%     [m,c] = ConvolutionalEncoder(L); % message m, codeword c
%     
%     % Binary Symmetric Channel - bsc
%     p = 10^-2; % Error Probability in BSC: 0<=p<=1
%     r = bsc(c,p); 
%     
%     m_est1 = ViterbiDecoder(r);
%     m_est2 = Viterbi_window15(r,L); % test
%     
%     % counting the errors
%     errViterbi(i) = size(find([m - m_est1(1:L)]),2);
%     errViterbi_w(i) = size(find([m - m_est2(1:L)]),2);
% end

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
    
    n = sigma/sqrt(2)*[randn(size(c)) + 1i*randn(size(c))]; % Gaussian noise, 0 mean and sigma^2 variance 
    % Noise addition
    r = s + 10^(-EN0_dB(i)/20)*n; % additive white gaussian noise

    % receiver - hard decision decoding
    r = real(r)<0;
    
    m_est1 = ViterbiDecoder(r);
    m_est2 = Viterbi_window15(r,L); % test
    
    % counting the errors
    errViterbi(i) = size(find([m - m_est1(1:L)]),2);
    errViterbi_w(i) = size(find([m - m_est2(1:L)]),2);
end

%% BER graphs
BER_Viterbi = errViterbi / L;
BER_Viterbi_window = errViterbi_w / L;

BER_theoretical = 0.5*erfc(sqrt(10.^(EbN0_dB/10))); % theoretical ber uncoded AWGN

figure
semilogy(EbN0_dB,BER_theoretical,'LineWidth',1.5);
hold on
semilogy(EbN0_dB,BER_Viterbi,'-d','LineWidth',1.5);
hold on
semilogy(EbN0_dB,BER_Viterbi_window,'-*','LineWidth',1.5);
axis([0 10 10^-5 0.5])
grid on
legend('BER-theoretical,uncoded', 'BER-Viterbi (n=2, K=3, k=1)','BER-Viterbi with window (l=15)');
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
function m_est = ViterbiDecoder(r)
    state = [0 0;0 1;1 0;1 1];
    survivorPath = zeros(4,length(r)/2);
    SM_k = zeros(4,1);      % State Metric: sum of Hamming distance, BM_k
    
    survivorPath(:,1)=[1;0;1;0];
    r_12 = r(1:4);
    r_12 = [r_12;r_12;r_12;r_12];   % for comparing to 'state'
    BM_k = sum(r_12 ~= [[0 0;1 1;1 1;0 0], state],2);   % branch metric
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
        r_k = [r_k;r_k;r_k;r_k];     % for comparing to 'state'
        BM_k = sum(r_k ~= state,2);  % Branch Metric: Hamming Distance 

        SM = SM_k;  % SM: k-1 th step
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
    currState = 1; % find(SM_k == min(SM_k));
    m_est = zeros(1,length(r)/2);
    for l = length(r)/2:-1:1
        prevState = survivorPath(currState,l); 
        m_est(l) = stateTable(currState,prevState);
        currState = prevState;
    end
end

%% Viterbi Decoder by windowing size of 15
function m_est = Viterbi_window15(r,L)
    window = 15;    % size of the window
    depth = L + 2;  % n * (L + M -1) / 2
    outer = ceil(depth/window);     % # of windows
    state = [0 0;0 1;1 0;1 1];
    
    entirePath = zeros(4,depth);
    m_est = zeros(1,depth);
    SM_k = zeros(4,1);  % State Metric: sum of Hamming distance, BM_k

    for out = 1:outer   
        ss = find(SM_k == min(SM_k));   % Finding the minimum hamming weight at the end of window
        ss = ss(1);
       
        if out ~= 1
            currState = ss;
            for l = length(survivorPath):-1:window-1
                prevState = survivorPath(currState,l); 
                currState = prevState;
            end
            ss = currState;     % Starting state of the next window
        end
        
        SM_k = zeros(4,1);  % Reset in every window
        
        if (out==1)||(depth <= 15)      % first window or short block
            window = min(depth/2, window);
            survivorPath = [];
            ss = 1;
            [survivorPath, SM_k, BM_k] = initial_two_steps(ss,1,r,survivorPath, SM_k, state);
            entirePath(:,1:2) = survivorPath;
            [survivorPath, SM_k, BM_k] = further_steps(3, 2*window, r, survivorPath, SM_k, state);
            entirePath(:,3:window) = survivorPath(:,1:window-2);
 
        elseif out == outer     % final window
            survivorPath = [];
            [survivorPath, SM_k, BM_k] = initial_two_steps(ss,2*window*(out-1)+1,r,survivorPath, SM_k, state);
            
            entirePath(:,window*(out-1)+1:window*(out-1)+2) = survivorPath;
            [survivorPath, SM_k, BM_k] = further_steps(window*(out-1)+3, depth, r, survivorPath, SM_k, state);
            entirePath(:,window*(out-1)+3:depth) = survivorPath;
            
        else    % mid-process
            survivorPath = [];
            [survivorPath, SM_k, BM_k] = initial_two_steps(ss,2*window*(out-1)+1,r,survivorPath, SM_k, state);
            entirePath(:,window*(out-1)+1:window*(out-1)+2) = survivorPath;

            rig = min(depth, 2*window*out);
            [survivorPath, SM_k, BM_k] = further_steps(window*(out-1)+3, rig, r, survivorPath, SM_k, state);
            entirePath(:,window*(out-1)+3:window*out) = survivorPath(:,1:window-2);
        end
    end
    ls = 1;
    m_est = traceback(ls, SM_k, entirePath);
end

% Initialization for each window
function [survivorPath, SM_k, BM_k] = initial_two_steps(ss,start,r,survivorPath, SM_k, state)
    list = [0 0 1 1 1 0 0 1;1 1 0 0 0 1 1 0;1 1 0 0 0 1 1 0;0 0 1 1 1 0 0 1];
    if (ss == 1)||(ss == 2)     % found states at the end of window
        survivorPath(:,1)= ss * ones(4,1);
        r_12 = r(start:start+3);
        r_12 = [r_12;r_12;r_12;r_12];   % for comparing to 'state'
        BM_k = sum(r_12 ~= [list(:,2*ss-1:2*ss), state],2);     % branch metric
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

    elseif (ss == 3)||(ss == 4)     % found states at the end of window
        survivorPath(:,1)= ss * ones(4,1);
        r_12 = r(start:start+3);
        r_12 = [r_12;r_12;r_12;r_12];   % for comparing to 'state'
        BM_k = sum(r_12 ~= [list(:,2*ss-1:2*ss), state],2);     % branch metric
        % state 00
        SM_k(1) = BM_k(4);
        survivorPath(1,2) = 2;
        % state 01
        SM_k(2) = BM_k(2);
        survivorPath(2,2) = 4;
        % state 10
        SM_k(3) = BM_k(1);
        survivorPath(3,2) = 2;
        % state 11
        SM_k(4) = BM_k(3);
        survivorPath(4,2) = 4;
    end
end

% Computation from 3rd steps for each window
function [survivorPath, SM_k, BM_k] = further_steps(lef, rig, r, survivorPath, SM_k, state)
    for k = lef:rig     % window range: from lef-th element to rig-th element
        r_k = r(2*k-1:2*k);
        r_k = [r_k;r_k;r_k;r_k];    % for comparing to 'state'
        BM_k = sum(r_k ~= state,2);     % Branch Metric: Hamming Distance 

        SM = SM_k;  % SM: k-1 th step
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

% Tracking the most-likely message from the survivor path
function part_est = traceback(ls, SM_k, survivorPath)
    stateTable = [ 0   0   0   0; 0   0   0   0; 1   1   0   0; 0   0   1   1 ]; 
    currState = ls; %find(SM_k == min(SM_k));
    [x y] = size(survivorPath);
    part_est = zeros(1,y);
    for l = y:-1:1
        prevState = survivorPath(currState,l); 
        part_est(l) = stateTable(currState,prevState);
        currState = prevState;
    end
end