%msg = [1 1 0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1 0 0]; 
%[d1 m1 p1 in1]=vitdec(code(1:end/2),trellis,tblen,'cont','hard') 

m = randi([0 1],70,1);
% Generate random binary data
codeword = convenc(m,trellis);



decodedData = vitdec(codeword,trellis,l,'trunc','hard');
% bit errors
biterr(m,decodedData)


% five two-bit symbols using a rate 2/3 convolutional code
data = randi([0 1],10,1);
trellis1 = poly2trellis([5 4],[23 35 0; 0 5 13]);
code1 = convenc(data,poly2trellis([5 4],[23 35 0; 0 5 13]));
%Verify that the encoded output is 15 bits, 3/2 times the length of the input sequence, data.
length(code1)

trellis2 = struct('numInputSymbols',2,'numOutputSymbols',4,...
'numStates',4,'nextStates',[0 2;0 2;1 3;1 3],...
'outputs',[0 3;1 2;3 0;2 1]);
code2 = convenc(randi([0 1],10,1),trellis2);

[code3,fstate] = convenc(data(1:6),trellis1);
code4 = convenc(data(7:10),trellis1,fstate);
isequal(code1,[code3; code4])

% r = 1/2
trellis = poly2trellis(7,{'1 + x^3 + x^4 + x^5 + x^6', ...
    '1 + x + x^3 + x^4 + x^6'})
data = randi([0 1],70,1);
% Generate random binary data, convolutionally encode the data, and decode the data using the Viterbi algorithm.
codedData = convenc(data,trellis);
decodedData = vitdec(codedData,trellis,34,'trunc','hard');
% bit errors
biterr(data,decodedData)


clear; close all
rng default
M = 64;                 % Modulation order
k = log2(M);            % Bits per symbol
EbNoVec = (4:10)';       % Eb/No values (dB)
numSymPerFrame = 1000;   % Number of QAM symbols per frame

berEstSoft = zeros(size(EbNoVec)); 
berEstHard = zeros(size(EbNoVec));

trellis = poly2trellis(7,[171 133]);
tbl = 32;
rate = 1/2;

for n = 1:length(EbNoVec)
    % Convert Eb/No to SNR
    snrdB = EbNoVec(n) + 10*log10(k*rate);
    % Noise variance calculation for unity average signal power.
    noiseVar = 10.^(-snrdB/10);
    % Reset the error and bit counters
    [numErrsSoft,numErrsHard,numBits] = deal(0);
    
    while numErrsSoft < 100 && numBits < 1e7
        % Generate binary data and convert to symbols
        dataIn = randi([0 1],numSymPerFrame*k,1);
        
        % Convolutionally encode the data
        dataEnc = convenc(dataIn,trellis);
        
        % QAM modulate
        txSig = qammod(dataEnc,M,'InputType','bit','UnitAveragePower',true);
        
        % Pass through AWGN channel
        rxSig = awgn(txSig,snrdB,'measured');
        
        % Demodulate the noisy signal using hard decision (bit) and
        % soft decision (approximate LLR) approaches.
        rxDataHard = qamdemod(rxSig,M,'OutputType','bit','UnitAveragePower',true);
        rxDataSoft = qamdemod(rxSig,M,'OutputType','approxllr', ...
            'UnitAveragePower',true,'NoiseVariance',noiseVar);
        
        % Viterbi decode the demodulated data
        dataHard = vitdec(rxDataHard,trellis,tbl,'cont','hard');
        dataSoft = vitdec(rxDataSoft,trellis,tbl,'cont','unquant');
        
        % Calculate the number of bit errors in the frame. Adjust for the
        % decoding delay, which is equal to the traceback depth.
        numErrsInFrameHard = biterr(dataIn(1:end-tbl),dataHard(tbl+1:end));
        numErrsInFrameSoft = biterr(dataIn(1:end-tbl),dataSoft(tbl+1:end));
        
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