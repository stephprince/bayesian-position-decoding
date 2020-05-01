function rippleseqout = getRippleSequences(trainingdata,testingdata,mincells, allindex,dayindex, processeddatadir)
% SP 9.27.18
% this function goes through all of the ripple events and decodes replay/preplay sequences

%set parameters
minspikes = 2; %minimum number of spikes in order to attempt to decode block
binsize = 0.015; %default temporal bin
recs = allindex(find(ismember(allindex(:,2),dayindex(1,2))),3); 
spikedatadir = [processeddatadir 'F' num2str(dayindex(1)) '_' num2str(dayindex(2)) '\Clusters\Sorted\'];
if ~isempty(recs)
    [clusterInfo recInfo] = getClusterInfo(spikedatadir, recs);
else
    rippleseqout = [];
    return
end

%% get estimated and actual position for each ripple event
for recIdx = 1:size(recInfo,2)
    if isstruct(testingdata{recIdx})
        for eventIdx = 1:length(testingdata{recIdx}.ripdata)
            %initialize variables
            activespiketimes = [];
            activerates = [];
            activecount = 0;
            
            %get start and end times to get time bins
            startevent = testingdata{recIdx}.riptime(eventIdx,1);
            endevent = testingdata{recIdx}.riptime(eventIdx,2);
            timebins{eventIdx} = startevent:binsize:endevent;
            
            %pick out matching units from the training data and the decoding data
            cellsdurtest = testingdata{recIdx}.ripdata(eventIdx).clusindex;
            cellsdurtrain = find(~isnan(trainingdata.rate(:,1)));
            matches = cellsdurtrain(ismember(cellsdurtrain,cellsdurtest));
            
            %decode position for each match
            for trainingcell = 1:length(matches)
                %pick out matching units that are stable during these event times
                cellstabletimes = trainingdata.stabletimesall{matches(trainingcell)}(recIdx,:);
                if sum(isExcluded([timebins{eventIdx}(1) timebins{eventIdx}(end)],cellstabletimes)) == 2 %both start and end time in the stable time periods
                    if matches(trainingcell) > 0 %if there is a cell in both the training and decoding data during this epoch
                        cellIdx = find(testingdata{recIdx}.ripdata(eventIdx).clusindex == matches(trainingcell));
                        tmpspiketimes = testingdata{recIdx}.ripdata(eventIdx).spiketimes(cellIdx);
                        
                        %save info for cells active during this epoc
                        if ~isempty(tmpspiketimes)
                            activecount = activecount+1;
                            activespiketimes{activecount} = tmpspiketimes;
                            activerates = [activerates; trainingdata.ratenormocc(matches(trainingcell),:)];
                        end
                    end
                end
            end
            activerates = activerates*binsize;
            
            %decode the data
            if ~isempty(activespiketimes)
                if activecount >= mincells
                %if (max(sum(cellfun(@(x) length(x), activespiketimes))) >= minspikes) old way of using min number of spikes, now at least 5 cells must be active during the ripple
                    %activespiketimes is the data to test, activerates is the training data
                    decodingOutput{recIdx}{eventIdx} = decodeSequencesDurRipples(activespiketimes, activerates, timebins{eventIdx});
                else
                    decodingOutput{recIdx}{eventIdx} = [];
                end
            else
                decodingOutput{recIdx}{eventIdx} = [];
            end
            
            %calc replay stats
            disp(['Calculating replay stats for rec ' num2str(recIdx) ' ripple ' num2str(eventIdx)])
            if ~isempty(decodingOutput{recIdx}{eventIdx})
                replayStats{recIdx}{eventIdx} = calcReplaySeqStats(decodingOutput{recIdx}{eventIdx});
            else
                replayStats{recIdx}{eventIdx} = [];
            end
            
            %calc replay stats by combining areas of the track with the same visual cues
            if ~isempty(decodingOutput{recIdx}{eventIdx})
                [replayStatsCues{recIdx}{eventIdx} decodingOutputCues{recIdx}{eventIdx}] = calcReplayStatsCues(decodingOutput{recIdx}{eventIdx});
            else
                replayStatsCues{recIdx}{eventIdx} = [];
                decodingOutputCues{recIdx}{eventIdx} = [];
            end
        end
    else
        replayStats{recIdx} = nan;
        decodingOutput{recIdx} = nan;
    end
end

%% concatenate data across recordings
allDecodingOutput = []; allReplayStats = []; allReplayStatsCues = []; allDecodingOutputCues = [];
eventCounter = 1;
for recIdx = 1:length(recs)
    if iscell(decodingOutput{recIdx})
        for eventIdx = 1:length(decodingOutput{recIdx})
            allDecodingOutput{eventCounter} = decodingOutput{recIdx}{eventIdx};
            allReplayStats{eventCounter} = replayStats{recIdx}{eventIdx};
            allReplayStatsCues{eventCounter} = replayStatsCues{recIdx}{eventIdx};
            allDecodingOutputCues{eventCounter} = decodingOutputCues{recIdx}{eventIdx};
            eventCounter = eventCounter + 1;
        end
    end
end

%% combine data for output
rippleseqout.decodingOutput = allDecodingOutput;
rippleseqout.replayStats = allReplayStats;
rippleseqout.replayStatsCues = allReplayStatsCues;
rippleseqout.decodingOutputCues = allDecodingOutputCues;
end