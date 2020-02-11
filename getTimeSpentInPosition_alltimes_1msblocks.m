function totalTimeInBin = getTimeSpentInPosition_alltimes_1msblocks(positionInfo, edges, recs,trainingtimes)
%this function gets the time spent in certain position bins of the virmen
%track for each recording. Use for normalizing place field firing
%SP 12.19.17

%% initialize variables
timeInBinPerRec = zeros(length(recs), length(edges)-1);
step = 0.02; %resampling frequency for position data

%% make total time in bin distribution
for recIdx = 1:size(recs,1)
    %get relevant virmen data
    thetaBreaks = positionInfo(recIdx).thetaBreaksSmooth;
    time = positionInfo(recIdx).timeSmooth;
    thetaSmooth = positionInfo(recIdx).thetaSmooth;
    timesWhenMoving = positionInfo(recIdx).whenMoving;
    
    %get times during 1ms training blocks
    timesWhenTrainingIdx = isExcluded(timesWhenMoving,trainingtimes{recIdx});
    timesWhenMovingAndTraining = timesWhenMoving(logical(timesWhenTrainingIdx));
    
    %go through theta and compile time in bin
    for binIdx = 1:length(edges)-1 
        for trialIdx = 1:length(thetaBreaks)-1
            thetaforTrial = thetaSmooth(thetaBreaks(trialIdx)+1:thetaBreaks(trialIdx+1));
            timesforTrial = time(thetaBreaks(trialIdx)+1:thetaBreaks(trialIdx+1));
            
            thetaforBinIdx = (thetaforTrial < edges(binIdx+1))  & (thetaforTrial >= edges(binIdx));
            thetaforBin = thetaforTrial(thetaforBinIdx);
            timesforBin = timesforTrial(thetaforBinIdx);
            
            %get timesforBin when speed greater than speed threshold
            timesforBinWhenMoving = ismember(timesforBin, timesWhenMovingAndTraining);
            
            %get time in bin for each trial
            samplesInBin = sum(timesforBinWhenMoving);
            timeInBinPerTrial(trialIdx) = samplesInBin*step;
        end
        %get time in bin for each recording
        timeInBinPerRec(recIdx, binIdx) = sum(timeInBinPerTrial);
    end
end

%% get time in bin for all recordings
totalTimeInBin = sum(timeInBinPerRec,1);
end