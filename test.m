% Viterbi Decoding Algorithms
clear all; clc;
% n = 2; K = 3; k=1
% L-bit message sequence
%% Convolutional Coding Encoder
L = 3;

    [m,c] = ConvolutionalEncoder(L);
%         s = 1 - 2*c; % 0 -> 1; 1 -> -1 % *sqrt(Eb)
%    sigma = 1; 
%     n = sigma/sqrt(2)*[randn(size(c)) + j*randn(size(c))]; % Gaussian noise, 0 mean and sigma^2 variance 
%     % Noise addition
%     y = s + n; % additive white gaussian noise
% 
%     % receiver - hard decision decoding
%     r = real(y)>0;
    m
    c
[survivorPath, m_est] = ViterbiDecoder(c)

function [m,c] = ConvolutionalEncoder(L)
    m = rand(1,L) > 0.5;
    g1 = [1 1 1];
    g2 = [1 0 1];
    upper = mod(conv(g1,m),2);
    lower = mod(conv(g2,m),2);
    c = [upper; lower]; 
    c = c(:).';     % codeword
end

function [survivorPath, m_est] = ViterbiDecoder(r)%window_size
    state = [0 0;0 1;1 0;1 1];
    survivorPath = zeros(4,length(r)/2);
    SM_k = zeros(4,1);  % State Metric
    
    survivorPath(:,1)=[1;0;1;0];
    r_12 = r(1:4);
    r_12 = [r_12;r_12;r_12;r_12]; % for comparing to 'state'
    BM_k = sum(r_12 ~= [[0 0;1 1;1 1;0 0], state],2); % branch metric
    % state 00
    SM_k(1) = BM_k(1);
%    survivorPath(1,2*k-1:2*k) = [1 1]*(idx-1);
    survivorPath(1,2) = 1;
    % state 01
    SM_k(2) = BM_k(3);
%    survivorPath(2,2*k-1:2*k) = [1 1*(idx-1)];
    survivorPath(2,2) = 3;
    % state 10
    SM_k(3) = BM_k(4);
%    survivorPath(3,2*k-1:2*k) = [0 1*(idx-1)];
    survivorPath(3,2) = 1;
    % state 11
    SM_k(4) = BM_k(2);
%    survivorPath(4,2*k-1:2*k) = [1 0*(idx-1)];
    survivorPath(4,2) = 3;
    
    for k = 3:length(r)/2
        r_k = r(2*k-1:2*k);
        r_k = [r_k;r_k;r_k;r_k]; % for comparing to 'state'
        BM_k = sum(r_k ~= state,2); % Hamming Distance: Branch Metric

        SM = SM_k; % SM: k-1 th step
        % state 00
        [SM_k(1), idx] = min([SM(1)+BM_k(1),SM(2)+BM_k(4)]);
%            survivorPath(1,2*k-1:2*k) = [1 1]*(idx-1);
        survivorPath(1,k) = idx;
        % state 01
        [SM_k(2), idx] = min([SM(3)+BM_k(3),SM(4)+BM_k(2)]);
%            survivorPath(2,2*k-1:2*k) = [1 1*(idx-1)];
        survivorPath(2,k) = idx;
        % state 10
        [SM_k(3), idx] = min([SM(1)+BM_k(4),SM(2)+BM_k(1)]);
%            survivorPath(3,2*k-1:2*k) = [0 1*(idx-1)];
        survivorPath(3,k) = idx;
        % state 11
        [SM_k(4), idx] = min([SM(3)+BM_k(2),SM(4)+BM_k(3)]);
%            survivorPath(4,2*k-1:2*k) = [1 0*(idx-1)];
        survivorPath(4,k) = idx;    
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