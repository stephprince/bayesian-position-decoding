function posEstForBinAvg = plotEstActPosition_avg(positionEstimate, posActual, eventTimes, savedfiguresdir, dayindex, recIdx,PFbinsize)
% SP 9.27.18
% this function averages all of the estimated positions within 2 degree
% bins of the actual positions to visualize the accuracy of the decoder

if isempty(positionEstimate)
    posEstForBinAvg = nan;
    return
elseif isempty(positionEstimate{recIdx})
    posEstForBinAvg = nan;
    return
end

%% compile all the data
plotbinsize = 0.02;
posActualForBlock = []; posEstForBlock = [];
for eventIdx = 1:size(positionEstimate{recIdx},2)
    if ~isempty(positionEstimate{recIdx}{eventIdx})
        blockSize = eventTimes{eventIdx}(1):plotbinsize:eventTimes{eventIdx}(end);
        posEstForBlock = [posEstForBlock repmat(positionEstimate{recIdx}{eventIdx},1,size(blockSize,2))];
        posActualForBlock = [posActualForBlock; posActual{recIdx}{eventIdx}];
    end
end

%% plot the averaged estimated position vs actual position
if ~isempty(posActualForBlock)
    %get average for each 2 degree bin of actual position
    posbins = [0:PFbinsize:360];
    [poshist,posedges,poshistidx] = histcounts(posActualForBlock,posbins);
    for binidx = 1:length(poshist)
        posEstForBin = posEstForBlock(:,find(poshistidx == binidx));
        posEstForBinAvg(:,binidx) = nanmean(posEstForBin,2);
    end
    
    %plot the data
    figure; imagesc('XData',[0:PFbinsize:360-PFbinsize],'YData',[0:2:360-PFbinsize],'CData',posEstForBinAvg)
    xlim([0 360]); ylim([0 360])
    xlabel('Actual position'); ylabel('Estimated position')
    animaldir = [savedfiguresdir 'F' num2str(dayindex(1)) '_' num2str(dayindex(2)) '\'];
    filename = [animaldir 'decodedposvsestpos_avg_rec' num2str(recIdx)];
    saveas(gcf,filename,'png'); saveas(gcf,filename,'fig');
end
end