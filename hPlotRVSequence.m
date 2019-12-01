% hPlotRVSequence Plots RV sequence history 
%   hPlotRVSequence(SNRIN,RVSEQPTRHISTORY,NFRAMES) plots the contents of
%   RVSEQPTRHISTORY for a number of SNR points specified in the input
%   vector SNRIN. NFRAMES specifies the number of frames to consider,
%   RVSEQPTRHISTORY is a cell array with one entry per SNR point. The
%   history of the values of the RV sequence pointer history is stored in a
%   matrix for the corresponding SNR point. This matrix can have one or two
%   rows for the one and two codeword case respectively. Each column of
%   this matrix stores the corresponding RV sequence pointer value for a
%   specific subframe.
%
%   The history for each SNR point is plotted in a different figure. 

%   Copyright 2016 The MathWorks, Inc.

function hPlotRVSequence(SNRIn,allRvSeqPtrHistory,NFrames)

% get number of codewords
ncw = size(allRvSeqPtrHistory{1},1);

for snrIdx = 1:numel(SNRIn)
    SNRdB = SNRIn(snrIdx);
    rvSeqPtrHistory = allRvSeqPtrHistory{snrIdx};
    % plot history for first codeword
    figure;
    line((0:(NFrames*10-1)),rvSeqPtrHistory(1,:),'Color','k',...
        'Marker','o','LineStyle','none')
    ax1 = gca; % current axes
    ax1.YLabel.String = 'cw1  RV sequence index';
    ax1.XLabel.String = 'subframe';
    ax1.YLim = [min(rvSeqPtrHistory(:))-1 max(rvSeqPtrHistory(:))+1];
    ax1.YTick = [0 1 2 3 4 5];
    % plot history for second codeword if it exists
    if ncw == 2
        ax1_pos = ax1.Position; % position of first axes
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

