
function out = calcspikestodecode_testingdata_stabletimes(processeddatadir, dayindex, behavioridx, testingtimes)
% SP 9.19.18
% this function gets spikes during 1ms blocks for decoding

cellcountthresh = 1; %minimum number of cells needed to try to decode the block

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

%get stable firing rate times
tempindex = behavioridx(behavioridx(:,2) == dayindex(2),:);
mintimestable = 5; windowsize = 5;
stableFRInfo = getstableclustertimes_gauss(tempindex, spikedatadir, windowsize, mintimestable, 0);

%% get all spike times during the testing blocks while moving
for recIdx = 1:size(recInfo,2)
    %get periods when moving
    timesWhenMoving = positionInfo(recIdx).whenMoving;
    
    %get testing periods when moving
    movetestIdx = isExcluded(timesWhenMoving,testingtimes{recIdx});
    moveTestOverlap = timesWhenMoving(logical(movetestIdx));
    
    stoppedMoveTestIdx = find(round(diff(moveTestOverlap),2) > 0.0200);
    startedMoveTestIdx = [1; stoppedMoveTestIdx+1];
    blocksToDecode{recIdx} = [moveTestOverlap(startedMoveTestIdx(1:end-1)) moveTestOverlap(stoppedMoveTestIdx)]; 
end

%% get spikecounts within decoding blocks 
for recIdx = 1:size(recInfo,2)
    celldata = [];
    spikecounts = [];
    if ~isempty(blocksToDecode{recIdx})
        for clusterIdx = 1:length(clusterInfo.labels)
            %get spike times to include
            spikeTimes = recInfo(recIdx).spikeTimes(clusterIdx).TS;
          
            %find spiketimes when firing rate was stable
            stableTimes = stableFRInfo.stabletimes{clusterIdx}(recIdx,:);
            stableSpikesIdx = isExcluded(spikeTimes,stableTimes);
            spikeTimesWhenStable = spikeTimes(logical(stableSpikesIdx));
            
            %spike times during testing periods
            testingSpikesIdx = isExcluded(spikeTimesWhenStable,blocksToDecode{recIdx});
            spikeTimesToDecode = spikeTimesWhenStable(logical(testingSpikesIdx));
 
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

