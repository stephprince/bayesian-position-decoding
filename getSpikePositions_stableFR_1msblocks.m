function spikePositions = getSpikePositions_stableFR_1msblocks(clusterIdx, recInfo, positionInfo, stableFRInfo, speedThreshold,trainingtimes)
% this function gets the corresponding virmen position for every spike time
% when the animal is moving above a certain speed threshold
% SP 12.19.17

%% initialize variable
spikePositions = [];

%% get spike positions for each recording
for recIdx = 1:size(recInfo,2)
    %find spiketimes when firing rate was stable
    spikeTimes = recInfo(recIdx).spikeTimes(clusterIdx).TS;
    stableTimes = stableFRInfo.stabletimes{clusterIdx}(recIdx,:);
    stableSpikesIdx = isExcluded(spikeTimes,stableTimes);
    spikeTimesWhenStable = spikeTimes(logical(stableSpikesIdx));
    
    %find spiketimes during the 1ms training blocks
    trainingSpikesIdx = isExcluded(spikeTimesWhenStable,trainingtimes{recIdx});
    spikeTimesWhenTraining = spikeTimesWhenStable(logical(trainingSpikesIdx));

    %find corresponding virmen time for every spike time
    positionTimes = positionInfo(recIdx).timeSmooth;
    spikeTimeVirmenTimeIdx = lookup2(spikeTimesWhenTraining,positionTimes);
    spikeTimesInVirmen = positionTimes(spikeTimeVirmenTimeIdx);
    
    %find spiketimes when animal was above speed threshold
    isMovingIdx = find(positionInfo(recIdx).speedSmooth > speedThreshold);
    timesMoving = positionInfo(recIdx).timeSmooth(isMovingIdx);
    spikeTimesMovingIdx = ismember(spikeTimesInVirmen, timesMoving);
    
    %get spikes and position only when above speed threshold
    spikesWhenMovingIdx = spikeTimeVirmenTimeIdx(spikeTimesMovingIdx);
    spikesWhenMovingPosition = positionInfo(recIdx).thetaSmooth(spikesWhenMovingIdx);
    
    %compile for all recordings
    spikePositions = [spikePositions; spikesWhenMovingPosition];
end

end