function [posActual positionEstimate timebins] = calcPositionEstimate(trainingdata,testingdata,positionInfo,dayindex,recIdx,binsize)
% SP 9.27.18
% this function goes through all of the testing blocks, finds cells that
% are place cells and also active in the testing blocks, and then decodes p
% osition using Bayesian pop decoding methods

%set parameters
minspikes = 2; %minimum number of spikes in order to attempt to decode block

%% get estimated and actual position for each event/block
for eventIdx = 1:length(testingdata{dayindex(1)}{dayindex(2)}{recIdx}.blockdata)
    %get actual position data for each rec and event index
    startevent = testingdata{dayindex(1)}{dayindex(2)}{recIdx}.blocktime(eventIdx,1);
    endevent = testingdata{dayindex(1)}{dayindex(2)}{recIdx}.blocktime(eventIdx,2);
    startpos = find(positionInfo(recIdx).timeSmooth == startevent);
    endpos =  find(positionInfo(recIdx).timeSmooth == endevent);
    posActual{recIdx}{eventIdx} = positionInfo(recIdx).thetaSmooth(startpos:endpos);
    
    %initialize variables
    disp(['Decoding block number ' num2str(eventIdx)])
    activespiketimes = [];
    activerates = [];
    activecount = 0;
            
    %pick out matching units from the training data and the decoding data
    cellsdurtest = testingdata{dayindex(1)}{dayindex(2)}{recIdx}.blockdata(eventIdx).clusindex;
    cellsdurtrain = find(~isnan(trainingdata{dayindex(1)}{dayindex(2)}.rate(:,1)));
    matches = cellsdurtrain(ismember(cellsdurtrain,cellsdurtest));
    %if binsize < 1
    if (endevent-startevent) < binsize
        continue
    else
        timebins{eventIdx} = startevent:binsize:endevent;
    end
    %else
    %    timebins{eventIdx} = [startevent endevent];
    %end
            
    %decode position for each match
    for trainingcell = 1:length(matches)
        %pick out matching units that are stable during these event times
        cellstabletimes = trainingdata{dayindex(1)}{dayindex(2)}.stabletimes{matches(trainingcell)}(recIdx,:);
        if sum(isExcluded(timebins{eventIdx},cellstabletimes)) == 2 %both start and end time in the stable time periods
            if matches(trainingcell) > 0 %if there is a cell in both the training and decoding data during this epoch
                cellIdx = find(testingdata{dayindex(1)}{dayindex(2)}{recIdx}.blockdata(eventIdx).clusindex == matches(trainingcell));
                tmpspiketimes = testingdata{dayindex(1)}{dayindex(2)}{recIdx}.blockdata(eventIdx).spiketimes(cellIdx);
                
                %save info for cells active during this epoch
                if ~isempty(tmpspiketimes)
                    activecount = activecount+1;
                    activespiketimes{activecount} = tmpspiketimes;
                    activerates = [activerates; trainingdata{dayindex(1)}{dayindex(2)}.ratenormocc(matches(trainingcell),:)];
                end
            end
        end
    end
    activerates = activerates*binsize;
            
    %decode the data
    if ~isempty(activespiketimes)
        if max(cellfun(@(x) length(x), activespiketimes)) >= minspikes
            %activespiketimes is the data to test, activerates is the training data
            decodingOutput{recIdx}{eventIdx} = decodePositionDurTheta(activespiketimes, activerates, timebins{eventIdx});
            positionEstimate{recIdx}{eventIdx} = sum(decodingOutput{recIdx}{eventIdx}.spatialprob,2);
        else
            decodingOutput{recIdx}{eventIdx} = [];
            positionEstimate{recIdx}{eventIdx} = [];
        end
    end
end

if ~exist('positionEstimate')
    positionEstimate = [];
end
end