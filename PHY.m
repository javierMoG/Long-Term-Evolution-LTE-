function [SNRIn, BER,trgh]=PHY(modulacion, nFr, espSnr,dispS)
%% PDSCH (Physical Downlink Shared Channel) modo TM1
% Utilizaremos la paqueteria LTEToolBox y funciones del ejemplo de Matlab 
% 'lte/PDSCHThroughputConformanceExample'
% Se calcula el rendimiento del canal y la Pb para diversos valores de SNR 
% por medio de la simulacion de la transimision de subtramas. 
% Para cada valor de SNR considerado, un malla fuente es generada y 
% modulada con OFDM para crear la onda de transmision. La onda generada se 
% pasa a traves de un canal de desvanecimiento ruidoso.
% En el receptor se realiza la estimacion del canal, la equalizacion, la
% demodulacion y la decodificacion. El rendimiento del canal del PSDCH se
% determina usando el resultado del bloque de CRC a la salida del
% decodificador del canal.

%% Transmisor
NFrames = nFr;%Numero de tramas generadas
SNRIn = [0:espSnr:18];% rango del SNR en dB

% Por simplicidad se utilizó un ancho de bando de 50 bloques fuente (9MHz)
% con asignacion completa y una tasa de codigo de 0.5.
% No se implementa el formato de decodificacion DCI 

simulationParameters = []; %Inicializamos los parametros de la simulacion   
simulationParameters.NDLRB = 50; %Ancho de banda (bloques fuente)        
simulationParameters.PDSCH.TargetCodeRate = 0.5; %Tasa de codigo
simulationParameters.PDSCH.PRBSet = (0:49)'; %Indices del PRB para cada 
%ranura de la trama


fprintf('\nTM1 - Single antenna (port 0)\n');
simulationParameters.PDSCH.TxScheme = 'Port0';
simulationParameters.PDSCH.DCIFormat = 'Format1'; %Se llena aunque no se 
%utilice para poder correr la simulacion
simulationParameters.CellRefP = 1; %numero de puertos de antena de señal
%específicos de celda
simulationParameters.PDSCH.Modulation = {modulacion};
simulationParameters.TotSubframes = 1; % Generar una subtrama a la vez
simulationParameters.PDSCH.CSI = 'On'; % Los bits suaves son pesados por 
% el CSI (Channel Status Information) 

enb = lteRMCDL(simulationParameters); %configuramos la antena de lte

rvSequence = enb.PDSCH.RVSeq; %Secuencia de version de redundancia 
trBlkSizes = enb.PDSCH.TrBlkSizes; %Guardamos el tamaño de los bloques de 
% transporte

ncw = length(string(enb.PDSCH.Modulation)); %numero de palabras de código
% (codewords)

%Desplegamos los parametros configurados para la simulacion
hDisplayENBParameterSummary(enb, 'TM1');

%% Configuracion del canal de propagacion

channel.Seed = 6;                    % Generador de numeros aleatorios
channel.NRxAnts = 1;                 % Una antena receptora
channel.DelayProfile = 'EPA';        % Perfil del retardo
channel.DopplerFreq = 5;             % Frecuencia Doppler
channel.MIMOCorrelation = 'Low';     % Correlacion de la antena y el UE
channel.NTerms = 16;                 % Numero de osiladores usados en el 
% fading 
channel.ModelType = 'GMEDS';         % Modelo de Rayleigh fading 
channel.InitPhase = 'Random';        % Fases iniciales aleatorias
channel.NormalizePathGains = 'On';   % Normalizamos la potencia del perfil 
% del retardo  
channel.NormalizeTxAnts = 'On';      % Normalizamos para antenas de 
% transmision

ofdmInfo = lteOFDMInfo(enb);         %Fijamos los parametros de la 
% modulacion OFDM
channel.SamplingRate = ofdmInfo.SamplingRate; %Fijamos la tasa de muestreo

%% Configuracion del Estimador del Canal

perfectChanEstimator = false; %Utilizamos un canal imperfecto
%El canal cambia lentamente en tiempo y frecuencia
%Parametros del estimador del canal
cec.PilotAverage = 'UserDefined';   % Establecemos el modo manual de 
% definicion de la ventana para leer los simbolos de referencia
cec.FreqWindow = 41;                % Tamaño de la ventana en frecuencia
cec.TimeWindow = 27;                % Tamaño de la ventana en tiempo
cec.InterpType = 'Cubic';           % Tipo de interpolacion usada con los 
%simbolos de referencia
cec.InterpWindow = 'Centered';      % Tipo de ventana de interpolacion 
cec.InterpWinSize = 1;              % Tamaño de la ventana de interpolacion
%Lo anterior se configura para reducir el efecto del ruido
%Para poder desplegar informacion de la transmision durante la simulacion
displaySimulationInformation = dispS;
%% Ciclo de procesamiento
% Para obtener el rendimiento y la Pb en cada valor de SNR, la informacion
% del PDSCH es analizada subtrama por subtrama 

% Arreglo para guardar el maximo rendimiento para cada valor de SNR 
maxThroughput = zeros(length(SNRIn),1); 
% Arreglo para guardar el redimiento de la simulacion para cada valor de
% SNR
simThroughput = zeros(length(SNRIn),1);

legendString = ['Throughput: ' char(enb.PDSCH.TxScheme)];
legendStringBer= [modulacion];
allRvSeqPtrHistory = cell(1,numel(SNRIn));
nFFT = ofdmInfo.Nfft; %Numero de puntos FTT 
totErr=[];
totTot=[];
BER=[];
trgh=[];
for snrIdx = 1:numel(SNRIn)
    % Fijamos el generador de numero aleatorios tal que dependa de la 
    % variable del ciclo para asegurarnos flujos aleatorios independientes
    rng(snrIdx,'combRecursive');
    
    SNRdB = SNRIn(snrIdx);
    fprintf('\nSimulando en %g dB SNR para %d Frame(s)\n' ,SNRdB, NFrames);
    
    offsets = 0;            % Inicializamos el valor del offset de la trama
    offset = 0;             % Inicializamo el valor del offset para las 
    % tramas de radio 
    blkCRC = [];            % Bloque CRC para todas las subtramas consideradas
    bitTput = [];           % Numero de bits recibidos con exito por subtrama
    txedTrBlkSizes = [];    % Numero de bits transmitidos por subtrama
     
    % Guardamos la historia de los valores del apuntador a los valores de
    % las secuencias RV para todos los procesos HARQ. Incializamos el
    % historial con NaNs porque algunas subtramas no tienen informacion
    rvSeqPtrHistory = NaN(ncw, NFrames*10);        
    
    % Inicializamos todos los procesos HARQ
    harqProcesses = hNewHARQProcess(enb);
    
    
    % Inicializamos las IDs de los procesos HARQ a 1 ya que el primer bloque 
    % de transporte distinto de cero siempre sera transmitido usando el
    % primer proceso HARQ
    harqProcessSequence = 1;

    for subframeNo = 0:(NFrames*10-1)
        
        % Actualizamos el numero de subtrama
        enb.NSubframe = subframeNo;

        % Obtenemos el ID del proceso HARQ para la subtrama de la secuencia
        % del proceso HARQ
        harqID = harqProcessSequence(mod(subframeNo, length(harqProcessSequence))+1);
                
        % Si hay un bloque de transporte no vacio en la subtrama actual
        % se realiza la transmision y recepcion
        if harqID == 0
            continue;
        end
        
        % Actualizacion del proceso HARQ actual
        harqProcesses(harqID) = hHARQScheduling( ...
            harqProcesses(harqID), subframeNo, rvSequence);

        % Extraemos el tamano del bloque de transporte
        trBlk = trBlkSizes(:, mod(subframeNo, 10)+1).';

        % Desplegamos la informacion de la simulacion en proceso
        if displaySimulationInformation
            disp(' ');
            disp(['Subframe: ' num2str(subframeNo)...
                            '. HARQ process ID: ' num2str(harqID)]);
        end
        
        % Actualizamos el apuntador de la tabla de la secuencia RV
        rvSeqPtrHistory(:,subframeNo+1) = ...
                               harqProcesses(harqID).txConfig.RVIdx.';
        % Actualizamos la configuracion de la transmision PDSCH con el
        % estado del proceso HARQ
        enb.PDSCH = harqProcesses(harqID).txConfig;      
        data = harqProcesses(harqID).data;

        % Creamos una forma de onda de transmision y obtenemos el ID de la
        % secuencua de programacion HARQ de la estructura de salida 'enbOut' 
        % que tambien contiene la configuracion de la forma de onda y los 
        % parametros de la modulacion OFDM

        [txWaveform,~,enbOut] = lteRMCDLTool(enb, data);
        % Agregamos 25 rellenos de muestra. Esto para cubrir el rango de
        % retardos esperados del modelado del canal
        txWaveform =  [txWaveform; zeros(25, 1)]; 
        
        % Obtenemos el ID de la secuencia HARQ de 'enbOut' para el
        % procesamiento HARQ
        harqProcessSequence = enbOut.PDSCH.HARQProcessSequence;

        % Inicializamos el tiempo de canal para cada subtrama
        channel.InitTime = subframeNo/1000;

        % Pasamos informacion a traves de modelo del canal
        rxWaveform = lteFadingChannel(channel, txWaveform);

        % Calculamos la ganancia del ruido incluyendo la compensacion para
        % por la asignacion de potencia del enlace descendente
        SNR = 10^((SNRdB-enb.PDSCH.Rho)/20);

        % Normalizamos la potencia del ruido para tomar en cuenta la tasa
        % de muestreo, la cual es una funcion de el tamano de IFFT usado en
        % la modulacion OFDM, y el numero de antenas
        N0 = 1/(sqrt(2.0*enb.CellRefP*double(nFFT))*SNR);

        % Creamos AWGN
        noise = N0*complex(randn(size(rxWaveform)), ...
                            randn(size(rxWaveform)));

        % Agregamos AWGN a la onda recibida                      
        rxWaveform = rxWaveform + noise;
        
        % Calculamos un nuevo valor del offset para una nueva trama
        if (mod(subframeNo,10) == 0)
            offset = lteDLFrameOffset(enb, rxWaveform);
            if (offset > 25)
                offset = offsets(end);
            end
            offsets = [offsets offset]; %#ok
        end
        
        %Sincronizamos la onda recibida
        rxWaveform = rxWaveform(1+offset:end, :);

        % Demodulacion OFDM en la informacion recibida para recrear la
        % malla fuente
        rxSubframe = lteOFDMDemodulate(enb, rxWaveform);

        % Estimacion del canal
        [estChannelGrid, noiseEst] = lteDLChannelEstimate( ...
                enb, enb.PDSCH, cec, rxSubframe);

        % Obtenemos los indices PDSCH
        pdschIndices = ltePDSCHIndices(enb, enb.PDSCH, enb.PDSCH.PRBSet);

        % Obtenemos los elementos fuente del PDSCH de la subtrama recibida
        [pdschRx, pdschHest] = lteExtractResources(pdschIndices, ...
            rxSubframe*(10^(-enb.PDSCH.Rho/20)), estChannelGrid);

        % Decodificamos el PDSCH
        dlschBits = ltePDSCHDecode(...
                             enb, enb.PDSCH, pdschRx, pdschHest, noiseEst);

        % Decodificamos el DL-SCH
        [decbits, harqProcesses(harqID).blkerr,harqProcesses(harqID).decState] = ...
            lteDLSCHDecode(enb, enb.PDSCH, trBlk, dlschBits, ...
                           harqProcesses(harqID).decState);

        % Desplegamos los errores de bloque
        if displaySimulationInformation
            if any(harqProcesses(harqID).blkerr)
                disp(['Block error. RV index: ' num2str(harqProcesses(harqID).txConfig.RVIdx)...
                      ', CRC: ' num2str(harqProcesses(harqID).blkerr)])
            else
                disp(['No error. RV index: ' num2str(harqProcesses(harqID).txConfig.RVIdx)...
                      ', CRC: ' num2str(harqProcesses(harqID).blkerr)])
            end
        end

        % Guardamos los valores para calcular el rendimiento
        % Solo tramas no vacias
        if any(trBlk)
            blkCRC = [blkCRC harqProcesses(harqID).blkerr]; 
            bitTput = [bitTput trBlk.*(1- ...
                harqProcesses(harqID).blkerr)]; 
            txedTrBlkSizes = [txedTrBlkSizes trBlk];
        end

    end
    
    % Calculamos el rendimiento maximo y el simulado 
    maxThroughput(snrIdx) = sum(txedTrBlkSizes); % Rendimiento maximo posible
    simThroughput(snrIdx) = sum(bitTput,2);      % Rendimiento simulado
    
    % Desplegamos los resultados
    fprintf([['\nThroughput(Mbps) for ', num2str(NFrames) ' Frame(s) '],...
        '= %.4f\n'], 1e-6*simThroughput(snrIdx)/(NFrames*10e-3));
    fprintf(['Throughput(%%) for ', num2str(NFrames) ' Frame(s) = %.4f\n'],...
        simThroughput(snrIdx)*100/maxThroughput(snrIdx));
    
    allRvSeqPtrHistory{snrIdx} = rvSeqPtrHistory;
    totTot= [totTot sum(txedTrBlkSizes)];
    totErr= [totErr sum(txedTrBlkSizes)-sum(bitTput)];
    trgh=[trgh 1e-6*simThroughput(snrIdx)*100/(NFrames*10e-3)];
end
BER= (totErr./totTot);

%% Resultados
% Grafica throughput

% Grafica BER
%figure(2)
end 