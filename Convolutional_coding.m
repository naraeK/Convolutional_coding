clear all; close all; clc
% n = 2; K = 3; k=1
% L-bit message sequence
%% Basic Encoder
m = [1 0 0 1 1];
g1 = [1 1 1];
g2 = [1 0 1];
upper = mod(conv(g1,m),2);
lower = mod(conv(g2,m),2);
c = []; % codeword
for i = 1:length(upper)
    c(2*i-1) = upper(i);
    c(2*i) = lower(i);
end

%% Built-in functions
trellis = poly2trellis(3,[7 5]);

M = 4;                 % Modulation order
k = log2(M);            % Bits per symbol
EbNoVec = (2:10)';       % Eb/No values (dB)
numSymPerFrame = 128;   % L; n_min % Number of QAM symbols per frame

berEstSoft = zeros(size(EbNoVec)); 
berEstHard = zeros(size(EbNoVec));

% Decode the data using the Viterbi algorithm
% window l <= m: tblen = 5 * K (the constraint length of the code)
l = 5*3; % tblen
rate = 1/2;

for n = 1:length(EbNoVec)
    % Convert Eb/No to SNR
    snrdB = EbNoVec(n) + 10*log10(k*rate);
    % Noise variance calculation for unity average signal power.
    noiseVar = 10.^(-snrdB/10);
    % Reset the error and bit counters
    [numErrsSoft,numErrsHard,numBits] = deal(0); % all value to 0
    
    while  numErrsSoft < 100 && numBits < 1e7
        % Generate binary data and convert to symbols
        message = randi([0 1],numSymPerFrame*k,1);
        
        % Convolutionally encode the data
        codeword = convenc(message,trellis);
        
        % QAM modulate
        txSig = qammod(codeword,M,'InputType','bit','UnitAveragePower',true);
        
        % Pass through AWGN channel
        rxSig = awgn(txSig,snrdB,'measured');
        
        % Demodulate the noisy signal using hard decision (bit) and
        % soft decision (approximate LLR) approaches.
        rxDataHard = qamdemod(rxSig,M,'OutputType','bit','UnitAveragePower',true);
        rxDataSoft = qamdemod(rxSig,M,'OutputType','approxllr', ...
            'UnitAveragePower',true,'NoiseVariance',noiseVar);
        
        % Viterbi decode the demodulated data
        dataHard = vitdec(rxDataHard,trellis,l,'cont','hard');
        dataSoft = vitdec(rxDataSoft,trellis,l,'cont','unquant');
        
        % Calculate the number of bit errors in the frame. Adjust for the
        % decoding delay, which is equal to the traceback depth.
        numErrsInFrameHard = biterr(message(1:end-l),dataHard(l+1:end));
        numErrsInFrameSoft = biterr(message(1:end-l),dataSoft(l+1:end));
        
        % Increment the error and bit counters
        numErrsHard = numErrsHard + numErrsInFrameHard;
        numErrsSoft = numErrsSoft + numErrsInFrameSoft;
        numBits = numBits + numSymPerFrame*k;

    end
    
    % Estimate the BER for both methods
    berEstSoft(n) = numErrsSoft/numBits;
    berEstHard(n) = numErrsHard/numBits;
end

semilogy(EbNoVec,[berEstSoft berEstHard],'-*')
hold on
semilogy(EbNoVec,berawgn(EbNoVec,'qam',M))
legend('Soft','Hard','Uncoded','location','best')
grid
xlabel('Eb/No (dB)')
ylabel('Bit Error Rate')