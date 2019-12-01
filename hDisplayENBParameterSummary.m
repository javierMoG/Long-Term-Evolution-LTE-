% hDisplayENBParameterSummary Print enb parameter summary 
%   hDisplayENBParameterSummary(ENB) prints a summary of some of the more
%   relevant simulation parameters based on the ENB structure.

%   Copyright 2016 The MathWorks, Inc.

function hDisplayENBParameterSummary(enb, txMode)

% Number of codewords
ncw = length(enb.PDSCH.Modulation);

% Number of transmit antennas, get them from the waveform size. We do not
% use the generated waveform here, only needed to calculate ntxants
[wave, ~] = lteRMCDLTool(enb,[]);
ntxants = size(wave,2);

% Print summary
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
