function out = calcspikeratenormocc_trainingdata_stabletimes(processeddatadir, dayindex, behaviorindex, sessindex, trainingtimes,binsize)
% SP 9.12.18
% this function gets the spike rates x time normalized for occupancy 

% OUTPUT
% ratenormocc -  firing rate normalized for occupancy
% occ - smoothed time in bin
% pos - position
% binsize - length of spatial bin

%% set data dirs
recs = behaviorindex(ismember(behaviorindex(:,2),dayindex(1,2)),3);
spikedatadir = [processeddatadir 'F' num2str(dayindex(1)) '_' num2str(dayindex(2)) '\Clusters\Sorted\'];
virmendatadir = [processeddatadir 'F' num2str(dayindex(1)) '_' num2str(dayindex(2)) '\' num2str(dayindex(3)) '\'];

%% load spike times and clusters
if ~isempty(recs)
    [clusterInfo recInfo] = getClusterInfo(spikedatadir, recs);
else
    out = [];
    return
end

%correct for stable firing rate times
tempindex = behaviorindex(behaviorindex(:,2) == dayindex(2),:);
mintimestable = 5; windowsize = 5;
FRInfotemp = getstableclustertimes_gauss(sessindex, spikedatadir, windowsize, mintimestable, 0);
sess2keep = find(ismember(sessindex(:,3), tempindex(:,3)));
FRInfo = FRInfotemp;
FRInfo.stabletimes = cellfun(@(x) x(sess2keep,:), FRInfotemp.stabletimes,'UniformOutput',0);

%% load virmen data
%get raw position file and smooth
positionInfoRaw = getVirmenPositionInfo(virmendatadir, dayindex(1), dayindex(2), recs);
positionInfo = smoothPosition(positionInfoRaw, recs);
speedThreshold = 1; 

%% get spike counts during movement
for clusterIdx = 1:length(clusterInfo.labels)
    %get position corresponding to each spike
    spikePositions(clusterIdx).position = getSpikePositions_stableFR_1msblocks(clusterIdx, recInfo, positionInfo, FRInfo, speedThreshold, trainingtimes);
    
    %get spike counts for each position
    edges = 0:binsize:360;
    spikePositions(clusterIdx).spikeCount = histcounts(spikePositions(clusterIdx).position, edges);
    
    %smooth the spike counts
    spikePositions(clusterIdx).spikeCountSmooth = gaussSmooth(spikePositions(clusterIdx).spikeCount,2);
    disp(['Smoothing spike counts for cluster ' num2str(clusterIdx) ' out of ' num2str(length(clusterInfo.labels))])
end
    
%% get time in each position (when above speed threshold)
clear totalTimeInBin
positionInfo = findWhenMoving(recInfo, positionInfo, speedThreshold);
for clusterIdx = 1:length(clusterInfo.labels)
    stableTimes = FRInfo.stabletimes{clusterIdx};
    totalTimeInBin(clusterIdx,:) = getTimeSpentInPosition_stableFR_1msblocks(positionInfo, stableTimes, edges, recs, trainingtimes);
    
    %smooth the occupancy
    occupancySmooth(clusterIdx,:) = gaussSmooth(totalTimeInBin(clusterIdx,:),2);
    
    % if low occupancy should exclude by replacing firing rate in low occ areas to nan
    lowocc = find(occupancySmooth(clusterIdx,:)<binsize*0.1);
    spikePositions(clusterIdx).spikeCountSmooth(lowocc) = nan;
    
    disp(['Smoothing ooccupancy for cluster ' num2str(clusterIdx) ' out of ' num2str(length(clusterInfo.labels))])
end

%% get spike rate normalized for occupancy

for clusterIdx = 1:length(clusterInfo.labels)
    spikePositionHistNorm =  spikePositions(clusterIdx).spikeCountSmooth./occupancySmooth(clusterIdx,:);
    out.ratenormocc(clusterIdx,:) = spikePositionHistNorm;
    out.rate(clusterIdx,:) = spikePositions(clusterIdx).spikeCountSmooth ;
    out.occ(clusterIdx,:) = occupancySmooth(clusterIdx,:);
    out.stabletimes{clusterIdx} = FRInfo.stabletimes{clusterIdx};
    out.clusterID = clusterInfo.labels;
    out.pos = edges;
    out.binsize = binsize;
end

end