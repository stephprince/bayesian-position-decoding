function [posEstAvg] = getdecodingdurtheta(sessindex,behavioridx,trainingdata,testingdata,savedfiguresdir, processeddatadir, params)
% SP 2.28.19
%this function uses Bayesian pop decoding to estimate position during
%testing times from training data (specifically during theta)

disp(['Decoding position during theta for sess ' num2str(sessindex(1,2))])
posEstAvg = [];

%% get actual position data
spikedatadir = [processeddatadir 'F' num2str(sessindex(1,1)) '_' num2str(sessindex(1,2)) '\Clusters\Sorted\'];
virmendatadir = [processeddatadir 'F' num2str(sessindex(1,1)) '_' num2str(sessindex(1,2)) '\' num2str(sessindex(1,4)) '\'];
recs = behavioridx(ismember(behavioridx(:,2),sessindex(1,2)),3);
if ~isempty(recs)
    positionInfoRaw = getVirmenPositionInfo(virmendatadir, sessindex(1,1), sessindex(1,2), recs);
    positionInfo = smoothPosition(positionInfoRaw, recs);
end
    
%% decode estimated position
for recIdx = 1:numel(recs)
    %find matching cells and decode position for each block
    [posActual positionEstimate eventTimes] = calcPositionEstimate190131(trainingdata,testingdata,positionInfo,sessindex(1,:),recIdx,params.decodebinsize);
    
    % plot the results
    if exist('positionEstimate')
        %plot the data in 100s blocks to look at decoder on trial by trial basis
        plotEstActPosition_100sblocks(positionEstimate, positionInfo, eventTimes, sessindex(1,:), savedfiguresdir, recIdx, params.PFbinsize)
        
        %plot the averaged estimated position vs actual position
        posEstAvg{recIdx} = plotEstActPosition_avg(positionEstimate, posActual, eventTimes, savedfiguresdir, sessindex(1,:), recIdx,params.PFbinsize)
        close all
    else
        posEstAvg{recIdx} = nan;
    end
    clear posActual positionEstimate
end
    
%concatenate data across recording files
posEstAvgtemp = [];
for recIdx = 1:numel(recs)
    if ~isnan(posEstAvg{recIdx})
        posEstAvgtemp = cat(3, posEstAvgtemp, posEstAvg{recIdx});
    end
end
posEstAvg_day = nanmean(posEstAvgtemp,3);

%plot the average for that day
figure; imagesc('XData',[0:params.PFbinsize:360-params.PFbinsize],'YData',[0:params.PFbinsize:360-params.PFbinsize],'CData',posEstAvg_day)
xlim([0 360]); ylim([0 360])
xlabel('Actual position'); ylabel('Estimated position')
animaldir = [savedfiguresdir 'F' num2str(sessindex(1,1)) '_' num2str(sessindex(1,2)) '\'];
filename = [animaldir 'decodedposvsestpos_avg_wholeday'];
saveas(gcf,filename,'png'); saveas(gcf,filename,'fig');
end