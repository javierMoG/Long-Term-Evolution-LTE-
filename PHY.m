%% PDSCH modo TM1

%% Introduccion
% Se calcula el rendimiento del canal para diversos valores de SNR por medio
% de la simulación de la transimision de subtramas. For each of the
% considered SNR points a populated resource grid is generated and OFDM
% modulated to create a transmit waveform. The generated waveform is passed
% through a noisy fading channel. The following operations are then
% performed by the receiver: channel estimation, equalization, demodulation
% and decoding. The throughput performance of the PDSCH is determined using
% the block CRC result at the output of the channel decoder.

%% Simulation Configuration
% The example is executed for a simulation length of 2 frames for a number
% of SNR points. A large number of |NFrames| should be used to produce
% meaningful throughput results. |SNRIn| can be an array of values or a
% scalar. 

NFrames = 2;                % Number of frames
SNRIn = [0 5 10 15];   % SNR range in dB



%% 
% For simplicity all TMs modeled in this example have a bandwidth of 50
% resource blocks with a full allocation and a code rate of 0.5. Not
% specifying an RMC number ensures that all downlink subframes are
% scheduled. If RMC is specified (e.g. 'R.0'), the subframe scheduling is
% as defined in TS 36.101 [ <#19 1> ] where subframe 5 is not scheduled in
% most cases.
% This example does not perform DCI
% format decoding, so the |DCIFormat| field is not strictly necessary.
% However, since the DCI format is closely linked to the TM, it is
% included for completeness.

simulationParameters = []; % clear simulationParameters   
simulationParameters.NDLRB = 50;        
simulationParameters.PDSCH.TargetCodeRate = 0.5;
simulationParameters.PDSCH.PRBSet = (0:49)';


fprintf('\nTM1 - Single antenna (port 0)\n');
simulationParameters.PDSCH.TxScheme = 'Port0';
simulationParameters.PDSCH.DCIFormat = 'Format1';
simulationParameters.CellRefP = 1;
simulationParameters.PDSCH.Modulation = {'16QAM'};

% Set other simulationParameters fields applying to all TMs
simulationParameters.TotSubframes = 1; % Generate one subframe at a time
simulationParameters.PDSCH.CSI = 'On'; % Soft bits are weighted by CSI

%%
% Call <matlab:doc('lteRMCDL') lteRMCDL> to generate the default eNodeB
% parameters not specified in |simulationParameters|. These will be
% required later to generate the waveform using <matlab:doc('lteRMCDLTool')
% lteRMCDLTool>.

enb = lteRMCDL(simulationParameters);

%%
% The output |enb| structure contains, amongst other fields, the transport
% block sizes and redundancy version sequence for each codeword subframe
% within a frame. These will be used later in the simulation.

rvSequence = enb.PDSCH.RVSeq;
trBlkSizes = enb.PDSCH.TrBlkSizes;

%%
% The number of codewords, |ncw|, is the number of entries in the
% |enb.PDSCH.Modulation| field.
ncw = length(string(enb.PDSCH.Modulation));



%%
% Next we print a summary of some of the more relevant simulation
% parameters. Check these values to make sure they are as expected. The
% code rate displayed can be useful to detect problems if manually
% specifying the transport block sizes. Typical values are 1/3, 1/2 and
% 3/4.

hDisplayENBParameterSummary(enb, 'TM1');

%% Propagation Channel Model Configuration
% The structure |channel| contains the channel model configuration
% parameters.

channel.Seed = 6;                    % Channel seed
channel.NRxAnts = 1;                 % 2 receive antennas
channel.DelayProfile = 'EPA';        % Delay profile
channel.DopplerFreq = 5;             % Doppler frequency
channel.MIMOCorrelation = 'Low';     % Multi-antenna correlation
channel.NTerms = 16;                 % Oscillators used in fading model
channel.ModelType = 'GMEDS';         % Rayleigh fading model type
channel.InitPhase = 'Random';        % Random initial phases
channel.NormalizePathGains = 'On';   % Normalize delay profile power  
channel.NormalizeTxAnts = 'On';      % Normalize for transmit antennas

% The sampling rate for the channel model is set using the value returned
% from <matlab:doc('lteOFDMInfo') lteOFDMInfo>.

ofdmInfo = lteOFDMInfo(enb);
channel.SamplingRate = ofdmInfo.SamplingRate;

%% Channel Estimator Configuration
% The variable |perfectChanEstimator| controls channel estimator behavior.
% Valid values are |true| or |false|. When set to |true| a perfect channel
% response is used as estimate, otherwise an imperfect estimation based on
% the values of received pilot signals is obtained.

% Perfect channel estimator flag
perfectChanEstimator = false;

%%
% If |perfectChanEstimator| is set to false a configuration structure |cec|
% is needed to parameterize the channel estimator. The channel changes
% slowly in time and frequency, therefore a large averaging window is used
% in order to average the noise out.

% Configure channel estimator
cec.PilotAverage = 'UserDefined';   % Type of pilot symbol averaging
cec.FreqWindow = 41;                % Frequency window size in REs
cec.TimeWindow = 27;                % Time window size in REs
cec.InterpType = 'Cubic';           % 2D interpolation type
cec.InterpWindow = 'Centered';      % Interpolation window type
cec.InterpWinSize = 1;              % Interpolation window size

%% Display Simulation Information
% The variable |displaySimulationInformation| controls the display of
% simulation information such as the HARQ process ID used for each
% subframe. In case of CRC error the value of the index to the RV sequence
% is also displayed.

displaySimulationInformation = true;

%% Processing Loop
% To determine the throughput at each SNR point, the PDSCH data is analyzed
% on a subframe by subframe basis using the following steps:
%
% * _Update Current HARQ Process._ The HARQ process either carries new
% transport data or a retransmission of previously sent transport data
% depending upon the Acknowledgment (ACK) or Negative Acknowledgment (NACK)
% based on CRC results. All this is handled by the HARQ scheduler,
% <matlab:edit('hHARQScheduling.m') hHARQScheduling.m>. The PDSCH data is
% updated based on the HARQ state.
%
% * _Create Transmit Waveform._ The data generated by the HARQ process is
% passed to <matlab:doc('lteRMCDLTool') lteRMCDLTool> which produces an
% OFDM modulated waveform, containing the physical channels and signals.
%
% * _Noisy Channel Modeling._ The waveform is passed through a fading
% channel and noise (AWGN) is added.
%
% * _Perform Synchronization and OFDM Demodulation._ The received symbols
% are offset to account for a combination of implementation delay and
% channel delay spread. The symbols are then OFDM demodulated.
%
% * _Perform Channel Estimation._ The channel response and noise levels are
% estimated. These estimates are used to decode the PDSCH.
%
% * _Decode the PDSCH._ The recovered PDSCH symbols for all transmit and
% receive antenna pairs, along with a noise estimate, are demodulated and
% descrambled by <matlab:doc('ltePDSCHDecode') ltePDSCHDecode> to obtain an
% estimate of the received codewords.
%
% * _Decode the Downlink Shared Channel (DL-SCH) and Store the Block CRC
% Error for a HARQ Process._ The vector of decoded soft bits is passed to
% <matlab:doc('lteDLSCHDecode') lteDLSCHDecode>; this decodes the codeword
% and returns the block CRC error used to determine the throughput of the
% system. The contents of the new soft buffer, |harqProc(harqID).decState|,
% is available at the output of this function to be used when decoding the
% next subframe.
%
% The number of transmit antennas P is obtained from the resource grid
% dimensions. 'dims' is M-by-N-by-P where M is the number of subcarriers, N
% is the number of symbols and P is the number of transmit antennas.
dims = lteDLResourceGridSize(enb);
P = dims(3);

% Initialize variables used in the simulation and analysis
% Array to store the maximum throughput for all SNR points
maxThroughput = zeros(length(SNRIn),1); 
% Array to store the simulation throughput for all SNR points
simThroughput = zeros(length(SNRIn),1);

% The temporary variables 'enb_init' and 'channel_init' are used to create
% the temporary variables 'enb' and 'channel' within the SNR loop to create
% independent simulation loops for the 'parfor' loop
enb_init = enb;
channel_init = channel;
legendString = ['Throughput: ' char(enb.PDSCH.TxScheme)];
allRvSeqPtrHistory = cell(1,numel(SNRIn));
nFFT = ofdmInfo.Nfft;
totErr=[];
totTot=[];
BER=[];

for snrIdx = 1:numel(SNRIn)
% To enable the use of parallel computing for increased speed comment out
% the 'for' statement above and uncomment the 'parfor' statement below.
% This needs the Parallel Computing Toolbox. If this is not installed
% 'parfor' will default to the normal 'for' statement. If 'parfor' is
% used it is recommended that the variable 'displaySimulationInformation'
% above is set to false, otherwise the simulation information displays for
% each SNR point will overlap.

    % Set the random number generator seed depending to the loop variable
    % to ensure independent random streams
    rng(snrIdx,'combRecursive');
    
    SNRdB = SNRIn(snrIdx);
    fprintf('\nSimulating at %g dB SNR for %d Frame(s)\n' ,SNRdB, NFrames);
    
    % Initialize variables used in the simulation and analysis
    offsets = 0;            % Initialize frame offset value
    offset = 0;             % Initialize frame offset value for radio frame
    blkCRC = [];            % Block CRC for all considered subframes
    bitTput = [];           % Number of successfully received bits per subframe
    txedTrBlkSizes = [];    % Number of transmitted bits per subframe
    enb = enb_init;         % Initialize RMC configuration
    channel = channel_init; % Initialize channel configuration
    %pmiIdx = 0;             % PMI index in delay queue
    
    % The variable harqPtrTable stores the history of the value of the
    % pointer to the RV sequence values for all the HARQ processes.
    % Pre-allocate with NaNs as some subframes do not have data
    rvSeqPtrHistory = NaN(ncw, NFrames*10);        
    
    % Initialize state of all HARQ processes
    harqProcesses = hNewHARQProcess(enb);
       
    % Initialize HARQ process IDs to 1 as the first non-zero transport
    % block will always be transmitted using the first HARQ process. This
    % will be updated with the full sequence output by lteRMCDLTool after
    % the first call to the function
    harqProcessSequence = 1;

    for subframeNo = 0:(NFrames*10-1)
        
        % Update subframe number
        enb.NSubframe = subframeNo;

        % Get HARQ process ID for the subframe from HARQ process sequence
        harqID = harqProcessSequence(mod(subframeNo, length(harqProcessSequence))+1);
                
        % If there is a transport block scheduled in the current subframe
        % (indicated by non-zero 'harqID'), perform transmission and
        % reception. Otherwise continue to the next subframe
        if harqID == 0
            continue;
        end
        
        % Update current HARQ process
        harqProcesses(harqID) = hHARQScheduling( ...
            harqProcesses(harqID), subframeNo, rvSequence);

        % Extract the current subframe transport block size(s)
        trBlk = trBlkSizes(:, mod(subframeNo, 10)+1).';

        % Display run time information
        if displaySimulationInformation
            disp(' ');
            disp(['Subframe: ' num2str(subframeNo)...
                            '. HARQ process ID: ' num2str(harqID)]);
        end
        
        % Update RV sequence pointer table
        rvSeqPtrHistory(:,subframeNo+1) = ...
                               harqProcesses(harqID).txConfig.RVIdx.';

        % Update the PDSCH transmission config with HARQ process state
        enb.PDSCH = harqProcesses(harqID).txConfig;      
        data = harqProcesses(harqID).data;

        % Create transmit waveform and get the HARQ scheduling ID sequence
        % from 'enbOut' structure output which also contains the waveform
        % configuration and OFDM modulation parameters
        [txWaveform,~,enbOut] = lteRMCDLTool(enb, data);
        
        % Add 25 sample padding. This is to cover the range of delays
        % expected from channel modeling (a combination of
        % implementation delay and channel delay spread)
        txWaveform =  [txWaveform; zeros(25, P)]; 
        
        % Get the HARQ ID sequence from 'enbOut' for HARQ processing
        harqProcessSequence = enbOut.PDSCH.HARQProcessSequence;

        % Initialize channel time for each subframe
        channel.InitTime = subframeNo/1000;

        % Pass data through channel model
        rxWaveform = lteFadingChannel(channel, txWaveform);

        % Calculate noise gain including compensation for downlink power
        % allocation
        SNR = 10^((SNRdB-enb.PDSCH.Rho)/20);

        % Normalize noise power to take account of sampling rate, which is
        % a function of the IFFT size used in OFDM modulation, and the 
        % number of antennas
        N0 = 1/(sqrt(2.0*enb.CellRefP*double(nFFT))*SNR);

        % Create additive white Gaussian noise
        noise = N0*complex(randn(size(rxWaveform)), ...
                            randn(size(rxWaveform)));

        % Add AWGN to the received time domain waveform        
        rxWaveform = rxWaveform + noise;

        % Once every frame, on subframe 0, calculate a new synchronization
        % offset
        if (mod(subframeNo,10) == 0)
            offset = lteDLFrameOffset(enb, rxWaveform);
            if (offset > 25)
                offset = offsets(end);
            end
            offsets = [offsets offset]; %#ok
        end
        
        % Synchronize the received waveform
        rxWaveform = rxWaveform(1+offset:end, :);

        % Perform OFDM demodulation on the received data to recreate the
        % resource grid
        rxSubframe = lteOFDMDemodulate(enb, rxWaveform);

        % Channel estimation
        if(perfectChanEstimator) 
            estChannelGrid = lteDLPerfectChannelEstimate(enb, channel, offset); %#ok
            noiseGrid = lteOFDMDemodulate(enb, noise(1+offset:end ,:));
            noiseEst = var(noiseGrid(:));
        else
            [estChannelGrid, noiseEst] = lteDLChannelEstimate( ...
                enb, enb.PDSCH, cec, rxSubframe);
        end

        % Get PDSCH indices
        pdschIndices = ltePDSCHIndices(enb, enb.PDSCH, enb.PDSCH.PRBSet);

        % Get PDSCH resource elements from the received subframe. Scale the
        % received subframe by the PDSCH power factor Rho. The PDSCH is
        % scaled by this amount, while the cell reference symbols used for
        % channel estimation (used in the PDSCH decoding stage) are not.
        [pdschRx, pdschHest] = lteExtractResources(pdschIndices, ...
            rxSubframe*(10^(-enb.PDSCH.Rho/20)), estChannelGrid);

        % Decode PDSCH
        dlschBits = ltePDSCHDecode(...
                             enb, enb.PDSCH, pdschRx, pdschHest, noiseEst);

        % Decode the DL-SCH
        [decbits, harqProcesses(harqID).blkerr,harqProcesses(harqID).decState] = ...
            lteDLSCHDecode(enb, enb.PDSCH, trBlk, dlschBits, ...
                           harqProcesses(harqID).decState);

        % Display block errors
        if displaySimulationInformation
            if any(harqProcesses(harqID).blkerr)
                disp(['Block error. RV index: ' num2str(harqProcesses(harqID).txConfig.RVIdx)...
                      ', CRC: ' num2str(harqProcesses(harqID).blkerr)])
            else
                disp(['No error. RV index: ' num2str(harqProcesses(harqID).txConfig.RVIdx)...
                      ', CRC: ' num2str(harqProcesses(harqID).blkerr)])
            end
        end

        % Store values to calculate throughput
        % Only for subframes with data
        if any(trBlk)
            blkCRC = [blkCRC harqProcesses(harqID).blkerr]; 
            bitTput = [bitTput trBlk.*(1- ...
                harqProcesses(harqID).blkerr)]; 
            txedTrBlkSizes = [txedTrBlkSizes trBlk];
        end

    end
    
    % Calculate maximum and simulated throughput
    maxThroughput(snrIdx) = sum(txedTrBlkSizes); % Max possible throughput
    simThroughput(snrIdx) = sum(bitTput,2);      % Simulated throughput
    
    % Display the results dynamically in the command window
    fprintf([['\nThroughput(Mbps) for ', num2str(NFrames) ' Frame(s) '],...
        '= %.4f\n'], 1e-6*simThroughput(snrIdx)/(NFrames*10e-3));
    fprintf(['Throughput(%%) for ', num2str(NFrames) ' Frame(s) = %.4f\n'],...
        simThroughput(snrIdx)*100/maxThroughput(snrIdx));
    
    allRvSeqPtrHistory{snrIdx} = rvSeqPtrHistory;
    totTot= [totTot sum(txedTrBlkSizes)];
    totErr= [totErr sum(txedTrBlkSizes)-sum(bitTput)];
end
BER= totErr./totTot;
%% Throughput Results
% The throughput results for the simulation are displayed in the MATLAB(R)
% command window after each SNR point is completed. They are also captured
% in |simThroughput| and |maxThroughput|. |simThroughput| is an array with
% the measured throughput in number of bits for all simulated SNR points.
% |maxThroughput| stores the maximum possible throughput in number of bits
% for each simulated SNR point.

% Plot throughput
figure 
plot(SNRIn, 1e-6*simThroughput/(NFrames*10e-3),'*-.');
xlabel('SNR (dB)');
ylabel('Throughput (Mbps)');
legend(legendString,'Location','NorthWest');
grid on;

%Plot BER
figure 
plot(SNRIn, BER,'*-.');
xlabel('SNR (dB)');
ylabel('BER');
legend(legendString,'Location','NorthWest');
grid on;
