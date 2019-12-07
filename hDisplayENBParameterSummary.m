% hDisplayENBParameterSummary Imrpime el resumen de los parámetros 
%   hDisplayENBParameterSummary(ENB) muestra un resumen de los 
%	parámetros de la simulación con base en ENB

%   Copyright 2016 The MathWorks, Inc.

function hDisplayENBParameterSummary(enb, txMode)

% Número de codewords
ncw = length(enb.PDSCH.Modulation);

%Número de antenas transmisoras
[wave, ~] = lteRMCDLTool(enb,[]);
ntxants = size(wave,2);

% Imprime resumen
fprintf('\n-- Parameter summary: --------------------------------------------\n');
disp(['                      Duplexing mode: ' enb.DuplexMode]);
disp(['                   Transmission mode: ' txMode]);
disp(['                 Transmission scheme: ' enb.PDSCH.TxScheme]);
disp(['  Number of downlink resource blocks: ' num2str(enb.NDLRB)]);
disp([' Number of allocated resource blocks: ' num2str(length(enb.PDSCH.PRBSet))]);
disp(['Cell-specific reference signal ports: ' num2str(enb.CellRefP)]);
disp(['         Number of transmit antennas: ' num2str(ntxants)]);
disp(['                 Transmission layers: ' num2str(enb.PDSCH.NLayers)]);
disp(['                 Number of codewords: ' num2str(ncw)]);
disp(['               Modulation codeword 1: ' enb.PDSCH.Modulation{1}]);
disp(['    Transport block sizes codeword 1: ' sprintf('%8d',enb.PDSCH.TrBlkSizes(1,:))]);
disp(['                Code rate codeword 1: ' sprintf('%8.4g',enb.PDSCH.ActualCodeRate(1,:))]);

if size(enb.PDSCH.TrBlkSizes,1)==2
    disp(['               Modulation codeword 2: ' enb.PDSCH.Modulation{2}]);
    disp(['    Transport block sizes codeword 2: ' sprintf('%8d',enb.PDSCH.TrBlkSizes(2,:))]);
    disp(['                Code rate codeword 2: ' sprintf('%8.4g',enb.PDSCH.ActualCodeRate(2,:))]);
end
disp('------------------------------------------------------------------');

end
