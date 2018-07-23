% Fading channel
% one relay; two hops
% Soft decision with Maximum likelihood (ML) decision
% Viterbi Decoding Algorithms
% n = 2; K = 3; k=1
% L-bit message sequence
clear all; clc;
%% Transmitter
L = 10^4-2;    % short packet

tic;
%% AWGN channel & Rayleigh fading channel
EbN0_dB = [0:1:15];
EN0_dB = EbN0_dB - 10*log10(2); % (1-bit -> 2-bit); E/N0 = 1/2 * Eb/N0
errViterbi = [];   
for i = 1:length(EN0_dB)
    % Transmitter
    m = [];
    [m,c] = ConvolutionalEncoder(L,m); % message m, codeword c
    % BPSK modulation
    s = 1 - 2*c; % 0 -> 1; 1 -> -1 % *sqrt(Eb)
    
    % AWGN noise addition
    n0 = 1/sqrt(2)*[randn(size(c)) + j*randn(size(c))]; % Gaussian noise, 0 mean and sigma^2 variance 
    n = 1/sqrt(2)*[randn(size(c)) + j*randn(size(c))];
 
    % Rayleigh channel fading
    h0 = 1/sqrt(2)*[randn(size(c)) + j*randn(size(c))];  % the channel gain of the first hop in the relay path
    h = 1/sqrt(2)*[randn(size(c)) + j*randn(size(c))];
    
% First Hop    
    % Send over Gaussian Link to the receiver
    r_fading0 = h0.*s + 10^(-EN0_dB(i)/20)*n0; % additive white gaussian noise

    % Equalization to remove fading effects. Ideal Equalization Considered
    r_fading0 = r_fading0./h0;
        
    % BPSK demodulator at the Receiver
    r_fading_hard0 = real(r_fading0) < 0;
    
    % Decoding
%     m_est02 = ViterbiDecoder(r_fading_hard0,'hard');
    m_est02 = Viterbi_window200(r_fading_hard0,L,'hard');
    m_est04 = Viterbi_window200(r_fading0,L,'soft');%ViterbiDecoder(r_fading0,'soft');
    
% Second Hop
    [m_est02,c2] = ConvolutionalEncoder(L,m_est02(1:L));
    [m_est04,c4] = ConvolutionalEncoder(L,m_est04(1:L));
    
    amp = 1; % Amplification Gain
    h1 = 1/sqrt(2)*[randn(size(c)) + j*randn(size(c))];  % the channel gain of the second hop in the relay path
    n1 = 1/sqrt(2)*[randn(size(c)) + j*randn(size(c))];
        
    r_fading_hard1 = h1.*amp.*(1-2*c2) + 10^(-EN0_dB(i)/20)*n1;
    r_fading1 = h1.*amp.*(1-2*c4) + 10^(-EN0_dB(i)/20)*n1;
    
    r_fading_hard1 = real(r_fading_hard1./h1) < 0;
    r_fading1 = r_fading1 ./ h1;
    
% no relay   
    r_fading = h.*s + 10^(-EN0_dB(i)/20)*n;
    r_fading_hard = real(r_fading./h) < 0;
    
    % Decoding
    m_est1 = ViterbiDecoder(r_fading_hard,'hard');
%     m_est2 = ViterbiDecoder(r_fading_hard1,'hard');
    m_est2 = Viterbi_window200(r_fading_hard1,L,'hard');
    m_est4 = Viterbi_window200(r_fading1,L,'soft');%ViterbiDecoder(r_fading1,'soft');
    
    % counting the errors
    errViterbi_h(i) = size(find([m - m_est1(1:L)]),2);
    errViterbi_fading(i) = size(find([m - m_est2(1:L)]),2);
%     errViterbi_soft(i) = size(find([m - m_est3(1:L)]),2);
    errViterbi_fading_soft(i) = size(find([m - m_est4(1:L)]),2);
end

%% BER graphs
BER_Viterbi_fading = errViterbi_fading / L;
BER_Viterbi_f_s = errViterbi_fading_soft / L;
BER_Viterbi_h = errViterbi_h / L;
% BER_Viterbi_s = errViterbi_soft / L;

% Rayleigh Tehoretical BER
SNR = 10.^(EbN0_dB/10);
BER_theor_fading = 0.5.*(1-sqrt(SNR./(SNR+1)));

figure
semilogy(EbN0_dB,BER_theor_fading,'LineWidth',1.5);
hold on
semilogy(EbN0_dB,BER_Viterbi_fading,'LineWidth',1.5);
hold on
semilogy(EbN0_dB,BER_Viterbi_f_s,'LineWidth',1.5);
hold on
semilogy(EbN0_dB,BER_Viterbi_h,'LineWidth',1.5);
% hold on
% semilogy(EbN0_dB,BER_Viterbi_s,'LineWidth',1.5);

axis([0 15 10^-4 0.5])
grid on
legend('Theoretical, uncoded', 'Viterbi(hard), 2-hop relay', 'Viterbi(soft), 2-hop relay', 'Viterbi(hard), no relay');
xlabel('Eb/No, dB');
ylabel('Bit Error Rate');
title('BER with Viterbi decoding for BPSK in Rayleigh fading');

toc;
%% Convolutional Coding Encoder
% n = 2; K = 3; k=1
% L-bit message sequence
function [m,c] = ConvolutionalEncoder(L,m)
    if length(m) == 0
        m = rand(1,L) > 0.5;
    end
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

    survivorPath(:,1) = [1;0;1;0];
    r_12 = real(r(1:4));
    r_12 = [r_12;r_12;r_12;r_12]; % for comparing to 'state'

    if mode == 'hard'
        state = [0 0;0 1;1 0;1 1];
        BM_k = sum(r_12 ~= [[0 0;1 1;1 1;0 0], state],2); % branch metric
    elseif mode == 'soft'
        state = [1 1;1 -1;-1 1;-1 -1];
        BM_k = abs(r_12-[[1 1;-1 -1;-1 -1;1 1],state]);%.^2; % branch metric
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

%% Viterbi Decoder by windowing size of 200 in relays
function m_est = Viterbi_window200(r,L,mode)
    window = 200;    % size of the window
    depth = L + 2;  % n * (L + M -1) / 2
    outer = ceil(depth/window);     % # of windows
    if mode == 'hard'
        state = [0 0;0 1;1 0;1 1];
    elseif mode == 'soft'
        state = [1 1;1 -1;-1 1;-1 -1];
    end
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
            ls = ss;
            m_est(window*(out-2)+1:window*(out-1)) = traceback(ls, blockPath);
        end
        
        SM_k = zeros(4,1);  % Reset in every window
        survivorPath = [];
        blockPath =[];
        if (out==1)||(depth <= window)      % first window or short block
            window = min(depth/2, window);
            ss = 1;
            [survivorPath, SM_k] = initial_two_steps(ss,1,r,survivorPath, SM_k, state,mode);
            blockPath(:,1:2) = survivorPath;
            [survivorPath, SM_k] = further_steps(3, 2*window, r, survivorPath, SM_k, state,mode);
            blockPath(:,3:window) = survivorPath(:,1:window-2);
 
        elseif out == outer     % final window
            [survivorPath, SM_k] = initial_two_steps(ss,2*window*(out-1)+1,r,survivorPath, SM_k, state,mode);            
            blockPath(:,1:2) = survivorPath;
            
            [survivorPath, SM_k] = further_steps(window*(out-1)+3, depth, r, survivorPath, SM_k, state,mode);
            blockPath(:,3:depth-window*(out-1)) = survivorPath;
 
        else    % mid-process
            [survivorPath, SM_k] = initial_two_steps(ss,2*window*(out-1)+1,r,survivorPath, SM_k, state,mode);
            blockPath(:,1:2) = survivorPath;

            rig = min(depth, 2*window*out);
            [survivorPath, SM_k] = further_steps(window*(out-1)+3, rig, r, survivorPath, SM_k, state,mode);
            blockPath(:,3:window) = survivorPath(:,1:window-2);
        end
    end
    ls = 1;
    m_est(window*(out-1)+1:depth) = traceback(ls, blockPath); 
end

% Initialization for each window
function [survivorPath, SM_k] = initial_two_steps(ss,start,r,survivorPath, SM_k, state,mode)
    list = [0 0 1 1 1 0 0 1;1 1 0 0 0 1 1 0;1 1 0 0 0 1 1 0;0 0 1 1 1 0 0 1];
    if (ss == 1)||(ss == 2)     % found states at the end of window
        survivorPath(:,1)= ss * ones(4,1);
        r_12 = r(start:start+3);
        r_12 = [r_12;r_12;r_12;r_12];   % for comparing to 'state'
        if mode == 'hard'
            BM_k = sum(r_12 ~= [list(:,2*ss-1:2*ss), state],2); % branch metric
        elseif mode == 'soft'
            list = 1 - 2.*list;
            BM_k = abs(r_12-[list(:,2*ss-1:2*ss),state]);
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

    elseif (ss == 3)||(ss == 4)     % found states at the end of window
        survivorPath(:,1)= ss * ones(4,1);
        r_12 = r(start:start+3);
        r_12 = [r_12;r_12;r_12;r_12];   % for comparing to 'state'
        if mode == 'hard'
            BM_k = sum(r_12 ~= [list(:,2*ss-1:2*ss), state],2); % branch metric
        elseif mode == 'soft'
            list = 1 - 2.*list;
            BM_k = abs(r_12-[list(:,2*ss-1:2*ss),state]);
            BM_k = sum(BM_k,2);
        end
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
function [survivorPath, SM_k] = further_steps(lef, rig, r, survivorPath, SM_k, state,mode)
    for k = lef:rig     % window range: from lef-th element to rig-th element
        r_k = r(2*k-1:2*k);
        r_k = [r_k;r_k;r_k;r_k];    % for comparing to 'state'
        if mode == 'hard'
            BM_k = sum(r_k ~= state,2); % Branch Metric: Hamming Distance 
        elseif mode == 'soft'           
            BM_k = sum(abs(r_k - state),2); % Branch Metric: Hamming Distance 
        end
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
function part_est = traceback(ls, survivorPath)
    stateTable = [ 0   0   0   0; 0   0   0   0; 1   1   0   0; 0   0   1   1 ]; 
    currState = ls; 
    currState = currState(1);
    [x y] = size(survivorPath);
    part_est = zeros(1,y);
    for l = y:-1:1
        prevState = survivorPath(currState,l); 
        part_est(l) = stateTable(currState,prevState);
        currState = prevState;
    end
end