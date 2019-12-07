% hPlotRVSequence gráfica el historial de la secuencia RV
%   hPlotRVSequence (SNRIN,RVSEQPTRHISTORY,NFRAMES) grafica los contenidos
%   RVSEQPTRHISTORY para un numero de puntos SNR especificados en el
%   vector entrada SNRIN. NFRAMES especifica el numero de tramas a
%   considerar, RVSEQPTHISTORY es una celda de arreglos con una entrada por
%   punto SNR. La historia de los valores de la historia de la secuencia de
%   apuntadores RV se guarda en una matrix para el correspondiente punto de 
%   SNR. Esta matrix puede tener una o dos filas para el caso de una o dos
%   palabras código respectivamente. Cada columna de está matrix guarda la
%   correspondiente secuencia de valores de apuntadores RV para una
%   subtrama especifica.
%
%   La historia de cada punto SNR se imprime en una figura diferente.
%   Copyright 2016 The MathWorks, Inc.

function hPlotRVSequence(SNRIn,allRvSeqPtrHistory,NFrames)

% obtiene el numero de palabras codigo
ncw = size(allRvSeqPtrHistory{1},1);

for snrIdx = 1:numel(SNRIn)
    SNRdB = SNRIn(snrIdx);
    rvSeqPtrHistory = allRvSeqPtrHistory{snrIdx};
    % imprime la historia para la primera palabra codigo
    figure;
    line((0:(NFrames*10-1)),rvSeqPtrHistory(1,:),'Color','k',...
        'Marker','o','LineStyle','none')
    ax1 = gca; % ejes actuales
    ax1.YLabel.String = 'cw1  RV sequence index';
    ax1.XLabel.String = 'subframe';
    ax1.YLim = [min(rvSeqPtrHistory(:))-1 max(rvSeqPtrHistory(:))+1];
    ax1.YTick = [0 1 2 3 4 5];
    %imprime la historia para la segunda palabra codigo si existe
    if ncw == 2
        ax1_pos = ax1.Position; % posicion del primer eje
        ax2 = axes('Position',ax1_pos,...
            'YAxisLocation','right','Color','none','YColor','r');
        ax2.YLabel.String = 'cw2  RV sequence index';
        ax2.YLim = ax1.YLim;
        ax2.YTick = [0 1 2 3 4 5];
        line((0:(NFrames*10-1)),rvSeqPtrHistory(2,:),'Parent',ax2,...
            'Color','r','Marker','+','LineStyle',':')
    end
    title(['RV Sequence Indices: SNR ' num2str(SNRdB) ' dB'])
end

end

