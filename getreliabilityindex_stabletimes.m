function out = getreliabilityindex_stabletimes(spikeratenormocc, PFbinsize, dayindex, behavioridx, sessindex, reliabilitythreshold, processeddatadir)
% SP 9.27.18
% this function gets the reliability index (ie. the prediction quality of
% the place field model, as a fraction of the explained variance)
% see Saleem et al. 2013 or 2018 for explanation and example

%% get directory info
disp(['Calculating reliability index - F ' num2str(dayindex(1)) '_' num2str(dayindex(2))])
spikedatadir = [processeddatadir 'F' num2str(dayindex(1)) '_' num2str(dayindex(2)) '\Clusters\Sorted\'];
virmendatadir = [processeddatadir 'F' num2str(dayindex(1)) '_' num2str(dayindex(2)) '\' num2str(dayindex(3)) '\'];
recs = behavioridx(ismember(behavioridx(:,2),dayindex(2)),3);
if ~isempty(recs)
    positionInfoRaw = getVirmenPositionInfo(virmendatadir, dayindex(1), dayindex(2), recs);
    positionInfo = smoothPosition(positionInfoRaw, recs);
end

%% split up data set into training and testing data
[trainingtimes testingtimes] = splittraintestdata(spikedatadir, recs, '8020');

%% get testing and training data
trainingdata = calcspikeratenormocc_trainingdata_stabletimes(processeddatadir, dayindex, behavioridx, sessindex, trainingtimes,PFbinsize);
testingdata = calcspikestodecode_testingdata(processeddatadir, dayindex, behavioridx, testingtimes); %bc you can't figure out if testing block is within a stable time until you loop through the clusters
%the stable time part is below when you calculate the reliability index

%% get the reliability index
for clusIdx = 1:size(spikeratenormocc.rate,1)
    counter = 1;
    for recIdx = 1:numel(recs)
        for eventIdx = 1:length(testingdata{recIdx}.blockdata)
            %get position to figure out what the firing rate is expected to be
            startevent = testingdata{recIdx}.blocktime(eventIdx,1);
            endevent = testingdata{recIdx}.blocktime(eventIdx,2);
            startpos = find(positionInfo(recIdx).timeSmooth == startevent);
            endpos =  find(positionInfo(recIdx).timeSmooth == endevent);
            posForEvent = positionInfo(recIdx).thetaSmooth(startpos:endpos);
            
            %check if event is within the cell's stable times
            cellstabletimes = trainingdata.stabletimes{clusIdx}(recIdx,:);
            if sum(isExcluded([startevent endevent],cellstabletimes)) == 2 %both start and end time in the stable time periods 
                %find the mean, actual, and expected firing rates
                posidx = lookup2(posForEvent, trainingdata.pos(1:end-1));
                meantrainingdata = nanmean(trainingdata.ratenormocc(clusIdx,:));
                expectedFR = nanmean(trainingdata.ratenormocc(clusIdx,posidx));
                spikeIdx = find(testingdata{recIdx}.blockdata(eventIdx).clusindex == clusIdx);
                actualFR = length(testingdata{recIdx}.blockdata(eventIdx).spiketimes(spikeIdx))/(endevent-startevent);
                
                %get prediction error and variance
                predictionerror(counter) = (actualFR - expectedFR)^2;
                variance(counter) = (actualFR - meantrainingdata)^2;
                counter = counter+1;
            end
        end
    end
    
    % get reliability index
    if exist('predictionerror')
        out(clusIdx) = 1 - (nansum(predictionerror)/nansum(variance));
        clear predictionerror variance
    end
end

end