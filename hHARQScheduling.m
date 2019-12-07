%hHARQScheduling Update HARQ process
%   HARQPROCESS = hHARQScheduling(HARQPROCESS,SF,RVSEQUENCE) actualiza el estado
%   del HARQ dado el número del sub frame y la sequencia de redundancia

%   Copyright 2009-2014 The MathWorks, Inc.
function harqProcess = hHARQScheduling(harqProcess,sf,rvSequence)
    
	blkpass = ~(harqProcess.blkerr);
    ncw = harqProcess.ncw;
    
    % Permite que una secuencia de RV se pueda compartir en dos codewords
    if(ncw==2 && size(rvSequence, 1)==1)
        rvSequence = repmat(rvSequence, 2, 1);
    end
    
    % Tamaño del bloque de transporte para el sub frame dado
    trblksize = harqProcess.txConfig.TrBlkSizes(:, mod(sf, 10)+1).';
    
    % Avanza al siguiente RV y genéra el nuevo dato
    harqProcess.txConfig.RVIdx = mod(harqProcess.txConfig.RVIdx + 1, ...
        size(rvSequence, 2)+1);
    blkpass(harqProcess.txConfig.RVIdx == 0) = true;
    if(ncw==2)
       if(length(trblksize)==1)
           trblksize = [trblksize trblksize];
       end
       tempData = [{randi([0 1], trblksize(1), 1)} ...
           {randi([0 1], trblksize(2), 1)}];
    else
       tempData = {randi([0 1], trblksize(blkpass(1)), 1)};
    end        
    
    % Transmite de nuevo si es necesario
    harqProcess.data(blkpass) = tempData(blkpass);
    
    % Limpia los buffers para la nueva transmision
    if(all(blkpass))
       harqProcess.decState = struct();
       [harqProcess.decState(1:ncw)]= deal(struct());
       [harqProcess.decState(1:ncw).CBSBuffers] = deal({});
       [harqProcess.decState(1:ncw).CBSCRC] = deal([]);
       [harqProcess.decState(1:ncw).BLKCRC] = deal([]);
    elseif(any(blkpass))
        harqProcess.decState(blkpass).CBSBuffers = {};
    end
    
    % Actualiza el RV para la siguiente transmisión
    harqProcess.txConfig.RVIdx(blkpass) = 1;
    
    % Actualiza el RV del lado del transmisor
    harqProcess.txConfig.RVSeq = rvSequence(1, ...
        harqProcess.txConfig.RVIdx(1)).';
    if(ncw==2)
        harqProcess.txConfig.RVSeq(2) = rvSequence(2, ...
            harqProcess.txConfig.RVIdx(2));
        harqProcess.txConfig.RVSeq = harqProcess.txConfig.RVSeq.';
    end
    
    % Actualiza el RV del decoder local
    harqProcess.txConfig.RV = harqProcess.txConfig.RVSeq;
    
end
