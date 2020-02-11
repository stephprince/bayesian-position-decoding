
function out = calcspikestodecode_testingdata(processeddatadir, dayindex, behavioridx, testingtimes)
% SP 9.19.18
% this function gets spikes during 1ms blocks for decoding

cellcountthresh = 2; %minimum number of cells needed to try to decode the block

%% set data directories
recs = behavioridx(ismember(behavioridx(:,2),dayindex(2)),3); 
spikedatadir = [processeddatadir 'F' num2str(dayindex(1)) '_' num2str(dayindex(2)) '\Clusters\Sorted\'];
virmendatadir = [processeddatadir 'F' num2str(dayindex(1)) '_' num2str(dayindex(2)) '\' num2str(dayindex(3)) '\'];

% get virmen and cluster info
if ~isempty(recs)
    [clusterInfo recInfo] = getClusterInfo(spikedatadir, recs);
else
    out = [];
    return
end
positionInfoRaw = getVirmenPositionInfo(virmendatadir, dayindex(1), dayindex(2), recs);
positionInfo = smoothPosition(positionInfoRaw, recs);
speedThreshold = 1; 
positionInfo = findWhenMoving(recInfo, positionInfo, speedThreshold);

%% get all spike times during the testing blocks while moving
for recIdx = 1:size(recInfo,2)
    %get periods to decode
    isMovingIdx = find(positionInfo(recIdx).speedSmooth > speedThreshold);
    timesMoving = positionInfo(recIdx).timeSmooth(isMovingIdx);
    movetestIdx = isExcluded(timesMoving,testingtimes{recIdx});
    moveTestOverlap = timesMoving(logical(movetestIdx));
    stoppedMoveTestIdx = find(round(diff(moveTestOverlap),2) > 0.0200);
    startedMoveTestIdx = [1; stoppedMoveTestIdx+1];
    blocksToDecode{recIdx} = [moveTestOverlap(startedMoveTestIdx(1:end-1)) moveTestOverlap(stoppedMoveTestIdx)]; %the periods will be 3.98 instead of 4 sec long max if the moving start times is offset
end

%% get spikecounts within decoding blocks 
for recIdx = 1:size(recInfo,2)
    celldata = [];
    spikecounts = [];
    if ~isempty(blocksToDecode{recIdx})
        for clusterIdx = 1:length(clusterInfo.labels)
            %get spike times to include
            spikeTimes = recInfo(recIdx).spikeTimes(clusterIdx).TS;
            testingSpikesIdx = isExcluded(spikeTimes,blocksToDecode{recIdx});
            spikeTimesToDecode = spikeTimes(logical(testingSpikesIdx));
            spikebins = periodAssign(spikeTimesToDecode,blocksToDecode{recIdx});
            
            %get spikes that occur within decoding blocks
            if ~isempty(spikeTimesToDecode)
                validspikes = find(spikebins);
                spikeTimesToDecode = spikeTimesToDecode(validspikes);
                spikebins = spikebins(validspikes);
            end
            
            %load data into matrix
            tmpcelldata = [spikeTimesToDecode' spikebins];
            tmpcelldata(:,3) = clusterIdx;
            celldata = [celldata; tmpcelldata]; %append cluster data
            
            %get spikecounts
            spikecount = zeros(1,size(blocksToDecode{recIdx},1));
            for i = 1:length(spikebins)
                spikecount(spikebins(i)) = spikecount(spikebins(i))+1;
            end
            spikecounts = [spikecounts; spikecount]; %number of spikes per block per cell, each row is a cell and each column is a block
        end
        
        celldata = sortrows(celldata,1); 
        cellcounts = sum((spikecounts > 0));
        eventindex = find(cellcounts >= cellcountthresh); %blocks with at least 2 units firing
    
        for blockIdx = 1:length(eventindex)
            tmpind = find(celldata(:,2) == eventindex(blockIdx));
            out{recIdx}.blocktime(blockIdx,1:2) = blocksToDecode{recIdx}(eventindex(blockIdx),[1 2]); %block time
            out{recIdx}.blockdata(blockIdx).spiketimes = celldata(tmpind,1); %spike times
            out{recIdx}.blockdata(blockIdx).clusindex = celldata(tmpind,3); %cell info
        end
    end
end
end

