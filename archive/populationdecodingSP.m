%% population decoding
% SP 9.12.18
% this script uses Bayesian population decoding to decode location using
% place cell spike trains

%% get index
animals = [1,7,8 9];

%set data directories and animals
savedfiguresdir = '\\neuro-cloud\labs\singer\Steph\Figures\decoding_std0peak0_reliability\';
processeddatadir = '\\neuro-cloud\labs\singer\ProcessedData\VR_AnnularTrack\';
spreadsheetdir = '\\neuro-cloud\labs\singer\Steph\Code\spreadsheets\VRAnnularTrackSpreadsheet.xlsx';
celltypedir = '\\neuro-cloud\labs\singer\Steph\Figures\celltype\celltypeProps\';

%all rec index
[allindex iden] = getdefaultindex5XFAD(processeddatadir, animals);
allindex2 = allindex(ismember(allindex(:,1),animals),[1:2 5]);
dayindex = unique(allindex2,'rows');

%behavior index (only recordings when running through virtual environment)
behavioridx = getBehaviorRecIndex(animals, spreadsheetdir);
allbehavioridx = allindex(ismember(allindex(:,[1:3 5]),behavioridx,'rows'),:);

stdthreshold = 0;
peakratethreshold = 0;
%% load place field info
for i = 1:size(dayindex) 
    if ~exist(savedfiguresdir); mkdir(savedfiguresdir); end;
    datadir = [savedfiguresdir 'F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '\'];
    filename = [datadir 'placefields_F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '.mat'];
    
        % get linearized spike rates x time normalized for occupancy
        spikeratenormocc{dayindex(i,1)}{dayindex(i,2)} = calcspikeratenormocc(processeddatadir, dayindex(i,:), behavioridx);
        
        if ~isempty(spikeratenormocc{dayindex(i,1)}{dayindex(i,2)})
            %exclude interneurons
            filename = [celltypedir 'cellTypeProps_' iden num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '.mat'];
            load(filename);
            cellTypeInfo = allrecProps{dayindex(i,1)}{dayindex(i,2)}.cellTypeInfo;
            fnames = fieldnames(spikeratenormocc{dayindex(i,1)}{dayindex(i,2)});
            for fieldIdx = 1:3
                temp = spikeratenormocc{dayindex(i,1)}{dayindex(i,2)}.(fnames{fieldIdx});
                temp(cellTypeInfo.INidx,:) = nan; %replace IN fields with nans
                spikeratenormocc{dayindex(i,1)}{dayindex(i,2)}.(fnames{fieldIdx}) = temp;
            end
            spikeratenormocc{dayindex(i,1)}{dayindex(i,2)}.exclINIdx = cellTypeInfo.INidx;
            
            % use peak rate threshold to exclude units
            peakrates = max(spikeratenormocc{dayindex(i,1)}{dayindex(i,2)}.ratenormocc,[],2);
            excl = find(peakrates < peakratethreshold); incl = find(peakrates > peakratethreshold);
            peakrateunits = spikeratenormocc{dayindex(i,1)}{dayindex(i,2)};
            for fieldIdx = 1:3
                temp = spikeratenormocc{dayindex(i,1)}{dayindex(i,2)}.(fnames{fieldIdx});
                temp(excl,:) = nan; %replace place cells with low peak FR with nans
                peakrateunits.(fnames{fieldIdx}) = temp;
            end
            peakrateunits.exclPeakIdx = excl;
            peakrateunits.inclPeakIdx = incl;
            
            % use peak rate + 2 std above the mean of peak to exclude units
            sd1 = nanstd(spikeratenormocc{dayindex(i,1)}{dayindex(i,2)}.ratenormocc,[],2);
            mean = nanmean(spikeratenormocc{dayindex(i,1)}{dayindex(i,2)}.ratenormocc,2);
            excl1 = find(peakrates < mean+(stdthreshold*sd1)); 
            excl2 = find(peakrates < peakratethreshold); 
            excl = unique([excl1; excl2]); 
            stdabovemeanunits = spikeratenormocc{dayindex(i,1)}{dayindex(i,2)};
            for fieldIdx = 1:3
                temp = spikeratenormocc{dayindex(i,1)}{dayindex(i,2)}.(fnames{fieldIdx});
                temp(excl,:) = nan; %replace place cells with low peak FR with nans
                stdabovemeanunits.(fnames{fieldIdx}) = temp;
            end
            stdabovemeanunits.exclStdIdx = excl;
            stdabovemeanunits.stdthreshold = stdthreshold;
            stdabovemeanunits.peakratethreshold = peakratethreshold;
            
            % use reliability threshold to exclude variable units (Saleem et al. 2013)
            
            % let's do this later bc I don't want to go through it rn
            reliableunits = [];
        else
            reliableunits = [];
            peakrateunits = [];
            stdabovemeanunits = [];
        end
        
        %plot all placefields to see span over track
        if ~isempty(stdabovemeanunits)
            for cellIdx = 1:size(stdabovemeanunits.ratenormocc,1)
                if ~isnan(stdabovemeanunits.ratenormocc(cellIdx,1))
                    maxind(cellIdx) = find(stdabovemeanunits.ratenormocc(cellIdx,:) == max(stdabovemeanunits.ratenormocc(cellIdx,:)));
                end
            end
            if ~exist('maxind'); maxind = []; end
            [sorted sortedind] = sort(maxind);
            sortedind(sorted == 0) = [];
            for cellIdx = 1:length(sortedind)
                stdabovemeanunits.sortedratenormocc(cellIdx,:) = stdabovemeanunits.ratenormocc(sortedind(cellIdx),:);
                stdabovemeanunits.sortedratenormocc_normfr(cellIdx,:) = stdabovemeanunits.ratenormocc(sortedind(cellIdx),:)./max(stdabovemeanunits.ratenormocc(sortedind(cellIdx),:));
            end
            if isfield(stdabovemeanunits,'sortedratenormocc_normfr')
                figure; hold on;
                imagesc(stdabovemeanunits.sortedratenormocc_normfr)
                title(['Place cells - std > ' num2str(stdthreshold) ' peakrate > ' num2str(peakratethreshold) ' - F' num2str(dayindex(i,1)) num2str(dayindex(i,2))])
                animaldir = [savedfiguresdir 'F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '\'];
                if ~exist(animaldir); mkdir(animaldir); end;
                filename = [animaldir 'placefields_sorted_std' num2str(stdthreshold) '_peakrate'];
                saveas(gcf,filename,'png'); saveas(gcf,filename,'fig');
            end
            
            % output - spike rates x time for place cell units
            placefields{dayindex(i,1)}{dayindex(i,2)}.reliabilityunits = reliableunits;
            placefields{dayindex(i,1)}{dayindex(i,2)}.peakrateunits = peakrateunits;
            placefields{dayindex(i,1)}{dayindex(i,2)}.peakstdunits = stdabovemeanunits;
        end
        clear maxind sortedind stdabovemeanunits
        
        filename = [datadir 'placefields_F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '.mat'];
        if ~exist(datadir); mkdir(datadir); end;
        save(filename, 'placefields');
    %end
end

%% decode position during theta
% use pop decoding on alternating 1ms blocks to decode theta position
minspikes = 2; %minimum number of spikes to attempt to decode a 1ms block

%split data up into alternating 1ms blocks
for i = 1:size(dayindex,1)
    spikedatadir = [processeddatadir 'F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '\Clusters\Sorted\'];
    recs = behavioridx(ismember(behavioridx(:,2),dayindex(i,2)),3);
    if ~isempty(recs)
        [clusterInfo recInfo] = getClusterInfo(spikedatadir, recs);
        for recIdx = 1:numel(recs)
            %%% THIS IS FOR 1 MS ALTERNATING BLOCKS
            %starttimes = 0:2:round(recInfo(recIdx).duration);
            %endtimes = 1:2:round(recInfo(recIdx).duration);
            
            %%% THIS IS FOR 80/20 4 MS TRAIN AND 1 TEST 
            tempstarttimes = 0:5:round(recInfo(recIdx).duration);
            tempendtimes = tempstarttimes-1;
            starttimes = tempstarttimes(1:end-1);
            endtimes = tempendtimes(2:end);
            
            if length(starttimes) > length(endtimes); starttimes = starttimes(1:end-1); end;
            trainingtimes{dayindex(i,1)}{dayindex(i,2)}{recIdx} = [starttimes', endtimes'];
            testingtimes{dayindex(i,1)}{dayindex(i,2)}{recIdx} = [endtimes(1:end-1)',starttimes(2:end)' ];
        end
    else
        trainingtimes{dayindex(i,1)}{dayindex(i,2)}{1} = [];
        testingtimes{dayindex(i,1)}{dayindex(i,2)}{1} = [];
    end
    clear clusterInfo recInfo
end

% get training data generated from 1ms blocks - linearized spike rates x spatial bins
for i = 1:size(dayindex,1)
    datadir = [savedfiguresdir 'F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '\'];
    filename = [datadir 'trainingdata_1msalternate_F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '.mat'];
    
    %if exist(filename)
    %    load(filename) 
    %else
        % get linearized spike rates x time normalized for occupancy
        placefields_1msblocks{dayindex(i,1)}{dayindex(i,2)} = calcspikeratenormocc_1msblocks(processeddatadir, dayindex(i,:), behavioridx, trainingtimes{dayindex(i,1)}{dayindex(i,2)});
        
        if ~isempty(placefields_1msblocks{dayindex(i,1)}{dayindex(i,2)})
            %exclude interneurons and non place fields
            fnames = fieldnames(placefields{dayindex(i,1)}{dayindex(i,2)}.peakrateunits);
            clear temp
            for fieldIdx = 1:3
                temp = placefields_1msblocks{dayindex(i,1)}{dayindex(i,2)}.(fnames{fieldIdx});
                temp(placefields{dayindex(i,1)}{dayindex(i,2)}.peakstdunits.exclINIdx,:) = nan; %replace IN fields with nans
                temp(placefields{dayindex(i,1)}{dayindex(i,2)}.peakstdunits.exclStdIdx,:) = nan; %replace non-place fields with nans
                placefields_1msblocks{dayindex(i,1)}{dayindex(i,2)}.(fnames{fieldIdx}) = temp;
            end
        
            %plot all placefields to see span over track
            for cellIdx = 1:size(placefields_1msblocks{dayindex(i,1)}{dayindex(i,2)}.ratenormocc,1)
                if ~isnan(placefields_1msblocks{dayindex(i,1)}{dayindex(i,2)}.ratenormocc(cellIdx,1))
                    maxind(cellIdx) = find(placefields_1msblocks{dayindex(i,1)}{dayindex(i,2)}.ratenormocc(cellIdx,:) == max(placefields_1msblocks{dayindex(i,1)}{dayindex(i,2)}.ratenormocc(cellIdx,:)));
                end
            end
            if ~exist('maxind'); maxind = []; end
            [sorted sortedind] = sort(maxind);
            sortedind(sorted == 0) = [];
            for cellIdx = 1:length(sortedind)
                placefields_1msblocks{dayindex(i,1)}{dayindex(i,2)}.sortedratenormocc(cellIdx,:) = placefields_1msblocks{dayindex(i,1)}{dayindex(i,2)}.ratenormocc(sortedind(cellIdx),:);
                placefields_1msblocks{dayindex(i,1)}{dayindex(i,2)}.sortedratenormocc_normfr(cellIdx,:) = placefields_1msblocks{dayindex(i,1)}{dayindex(i,2)}.ratenormocc(sortedind(cellIdx),:)./max(placefields_1msblocks{dayindex(i,1)}{dayindex(i,2)}.ratenormocc(sortedind(cellIdx),:));
            end
            if isfield(placefields_1msblocks{dayindex(i,1)}{dayindex(i,2)},'sortedratenormocc_normfr')
                figure; hold on;
                imagesc(placefields_1msblocks{dayindex(i,1)}{dayindex(i,2)}.sortedratenormocc_normfr)
                title(['Place cells training data - std > ' num2str(stdthreshold) ' peakrate > ' num2str(peakratethreshold) ' - F' num2str(dayindex(i,1)) num2str(dayindex(i,2))])
                animaldir = [savedfiguresdir 'F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '\'];
                if ~exist(animaldir); mkdir(animaldir); end;
                filename = [animaldir 'placefields_sorted_traindata_std' num2str(stdthreshold) '_peakrate'];
                saveas(gcf,filename,'png'); saveas(gcf,filename,'fig');
            end
            clear maxind sortedind
        end
        
        filename = [datadir 'trainingdata_1msalternate_F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '.mat'];
        if ~exist(datadir); mkdir(datadir); end;
        save(filename, 'placefields_1msblocks');
    %end 
end

% get decoding data generated from 1ms blocks - spike counts x temporal bins 
for i = 1:size(dayindex,1)
    datadir = [savedfiguresdir 'F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '\'];
    filename = [datadir 'testinggdata_1msalternate_F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '.mat'];
    
    %if exist(filename)
    %    load(filename) 
    %else
        % get spikecounts over temporal bins to be decoded
        testingdata_1msblocks{dayindex(i,1)}{dayindex(i,2)} = calcspikestodecode_1msblocks(processeddatadir, dayindex(i,:), behavioridx, testingtimes{dayindex(i,1)}{dayindex(i,2)});
        
        if ~exist(datadir); mkdir(datadir); end;
        save(filename, 'testingdata_1msblocks');
   % end 
end

%get the reliability of the units using the training and test data
for i = 1:size(dayindex,1)
    if ~isempty(placefields_1msblocks{dayindex(i,1)}{dayindex(i,2)})
        %get actual position data
        spikedatadir = [processeddatadir 'F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '\Clusters\Sorted\'];
        virmendatadir = [processeddatadir 'F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '\' num2str(dayindex(i,3)) '\'];
        recs = behavioridx(ismember(behavioridx(:,2),dayindex(i,2)),3);
        if ~isempty(recs)
            positionInfoRaw = getVirmenPositionInfo(virmendatadir, dayindex(i,1), dayindex(i,2), recs);
            positionInfo = smoothPosition(positionInfoRaw, recs);
        end
        for clusIdx = 1:size(placefields_1msblocks{dayindex(i,1)}{dayindex(i,2)}.rate,1)
            counter = 1;
            clusIdx 
            for recIdx = 1:numel(recs)
                for eventIdx = 1:length(testingdata_1msblocks{dayindex(i,1)}{dayindex(i,2)}{recIdx}.blockdata)
                    %see if cluster is active within this event
                    if ismember(clusIdx, testingdata_1msblocks{dayindex(i,1)}{dayindex(i,2)}{recIdx}.blockdata(eventIdx).clusindex) && ~isnan(placefields_1msblocks{dayindex(i,1)}{dayindex(i,2)}.ratenormocc(clusIdx,1))
                        %get position to figure out what the firing rate is expected to be
                        startevent = testingdata_1msblocks{dayindex(i,1)}{dayindex(i,2)}{recIdx}.blocktime(eventIdx,1);
                        endevent = testingdata_1msblocks{dayindex(i,1)}{dayindex(i,2)}{recIdx}.blocktime(eventIdx,2);
                        startpos = find(positionInfo(recIdx).timeSmooth == startevent);
                        endpos =  find(positionInfo(recIdx).timeSmooth == endevent);
                        posForEvent = positionInfo(recIdx).thetaSmooth(startpos:endpos);
                        
                        % get linearized spike rates x time normalized for occupancy
                        posidx = lookup2(posForEvent, placefields_1msblocks{dayindex(i,1)}{dayindex(i,2)}.pos(1:end-1));
                        meantrainingdata = nanmean(placefields_1msblocks{dayindex(i,1)}{dayindex(i,2)}.ratenormocc(clusIdx,:));
                        for j = 1:length(posidx)
                            t = posidx(j);
                            expectedFR = placefields_1msblocks{dayindex(i,1)}{dayindex(i,2)}.ratenormocc(clusIdx,t);
                            spikeIdx = find(testingdata_1msblocks{dayindex(i,1)}{dayindex(i,2)}{recIdx}.blockdata(eventIdx).clusindex == clusIdx);
                            actualFR = length(testingdata_1msblocks{dayindex(i,1)}{dayindex(i,2)}{recIdx}.blockdata(eventIdx).spiketimes(spikeIdx))/(endevent-startevent);
                            
                            predictionerror(counter) = (actualFR - expectedFR)^2;
                            variance(counter) = (actualFR - meantrainingdata)^2;
                            counter = counter+1;
                        end
                    end
                end  
            end
            if exist('predictionerror')
                % get reliability index
                reliability{dayindex(i,1)}{dayindex(i,2)}(clusIdx) = 1 - (nansum(predictionerror)/nansum(variance));
                clear predictionerror variance
            end
        end
    end
end

%eliminate place fields using a reliability threshold
reliabilitythreshold = 0.01;
for i = 1:size(dayindex,1)       
    %get units below threshold
    excl1 = find(reliability{dayindex(i,1)}{dayindex(i,2)} < reliabilitythreshold);
    excl2 = find(isnan(reliability{dayindex(i,1)}{dayindex(i,2)}));
    exclIdx = unique([excl1 excl2]);
    
    %replace unit firing rate info with a nan
    if ~isempty(placefields_1msblocks{dayindex(i,1)}{dayindex(i,2)})
        fnames = fieldnames(placefields_1msblocks{dayindex(i,1)}{dayindex(i,2)});
        clear temp
        for fieldIdx = 1:3
            temp = placefields_1msblocks{dayindex(i,1)}{dayindex(i,2)}.(fnames{fieldIdx});
            temp(exclIdx,:) = nan; %replace IN fields with nans
            placefields_1msblocks{dayindex(i,1)}{dayindex(i,2)}.(fnames{fieldIdx}) = temp;
        end
        
        %plot place fields with reliability idx
        for cellIdx = 1:size(placefields_1msblocks{dayindex(i,1)}{dayindex(i,2)}.ratenormocc,1)
            if ~isnan(placefields_1msblocks{dayindex(i,1)}{dayindex(i,2)}.ratenormocc(cellIdx,1))
                maxind(cellIdx) = find(placefields_1msblocks{dayindex(i,1)}{dayindex(i,2)}.ratenormocc(cellIdx,:) == max(placefields_1msblocks{dayindex(i,1)}{dayindex(i,2)}.ratenormocc(cellIdx,:)));
            end
        end
        if ~exist('maxind'); maxind = []; end
        [sorted sortedind] = sort(maxind);
        sortedind(sorted == 0) = [];
        for cellIdx = 1:length(sortedind)
            placefields_1msblocks{dayindex(i,1)}{dayindex(i,2)}.sortedratenormocc_reliable(cellIdx,:) = placefields_1msblocks{dayindex(i,1)}{dayindex(i,2)}.ratenormocc(sortedind(cellIdx),:);
            placefields_1msblocks{dayindex(i,1)}{dayindex(i,2)}.sortedratenormocc_reliable_normfr(cellIdx,:) = placefields_1msblocks{dayindex(i,1)}{dayindex(i,2)}.ratenormocc(sortedind(cellIdx),:)./max(placefields_1msblocks{dayindex(i,1)}{dayindex(i,2)}.ratenormocc(sortedind(cellIdx),:));
        end
        if isfield(placefields_1msblocks{dayindex(i,1)}{dayindex(i,2)},'sortedratenormocc_reliable_normfr')
            figure; hold on;
            imagesc(placefields_1msblocks{dayindex(i,1)}{dayindex(i,2)}.sortedratenormocc_reliable_normfr)
            title(['Place cells training data - reliability > ' num2str(reliabilitythreshold) ' - F' num2str(dayindex(i,1)) num2str(dayindex(i,2))])
            animaldir = [savedfiguresdir 'F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '\'];
            filename = [animaldir 'placefields_sorted_traindata_reliable'];
            saveas(gcf,filename,'png'); saveas(gcf,filename,'fig');
        end
        clear maxind sortedind
        
        filename = [datadir 'trainingdata_1msalternate_reliability_F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '.mat'];
        if ~exist(datadir); mkdir(datadir); end;
        save(filename, 'placefields_1msblocks');
    end
end

%use training and decoding data to decode position 
binsize = 1; %temporal binsize
for i = 1:size(dayindex,1)
    %get actual position data
    spikedatadir = [processeddatadir 'F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '\Clusters\Sorted\'];
    virmendatadir = [processeddatadir 'F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '\' num2str(dayindex(i,3)) '\'];
    recs = behavioridx(ismember(behavioridx(:,2),dayindex(i,2)),3);
    if ~isempty(recs)
        positionInfoRaw = getVirmenPositionInfo(virmendatadir, dayindex(i,1), dayindex(i,2), recs);
        positionInfo = smoothPosition(positionInfoRaw, recs);
    end
    
    %decode estimated position 
    for recIdx = 1:numel(recs)         
        for eventIdx = 1:length(testingdata_1msblocks{dayindex(i,1)}{dayindex(i,2)}{recIdx}.blockdata)
            %initialize variables'
            disp(['Decoding block number ' num2str(eventIdx)])
            trainingdata = [];
            spikedata = [];
            decodedata = [];
            activespiketimes = [];
            activerates = [];
            posvect = placefields_1msblocks{dayindex(i,1)}{dayindex(i,2)}.pos;
            
            %pick out matching units from the training data and the decoding data
            cellsdurtest = testingdata_1msblocks{dayindex(i,1)}{dayindex(i,2)}{recIdx}.blockdata(eventIdx).clusindex;
            cellsdurtrain = find(~isnan(placefields_1msblocks{dayindex(i,1)}{dayindex(i,2)}.rate(:,1)));
            matches = cellsdurtrain(ismember(cellsdurtrain,cellsdurtest));
            startevent = testingdata_1msblocks{dayindex(i,1)}{dayindex(i,2)}{recIdx}.blocktime(eventIdx,1);
            endevent = testingdata_1msblocks{dayindex(i,1)}{dayindex(i,2)}{recIdx}.blocktime(eventIdx,2);
            if binsize < 1 %%|| (endevent-startevent > 1)
                timebins{eventIdx} = startevent:binsize:endevent;
            else
                timebins{eventIdx} = [startevent endevent];
            end
                       
            eventcellsactive = [];
            activecount = 0;
            for trainingcell = 1:length(matches)
                %pick out matching units that are stable during these event times
                %cellstabletimes = placefields_1msblocks{dayindex(i,1)}{dayindex(i,2)}.stabletimes{matches(trainingcell)}(recIdx,:);
                %if sum(isExcluded(timebins{eventIdx},cellstabletimes)) == 2 %both start and end time in the stable time periods
                    if matches(trainingcell) > 0 %if there is a cell in both the training and decoding data during this epoch
                        trainingdata = [trainingdata; placefields_1msblocks{dayindex(i,1)}{dayindex(i,2)}.ratenormocc(matches(trainingcell),:)];
                        cellIdx = find(testingdata_1msblocks{dayindex(i,1)}{dayindex(i,2)}{recIdx}.blockdata(eventIdx).clusindex == matches(trainingcell));
                        tmpspiketimes = testingdata_1msblocks{dayindex(i,1)}{dayindex(i,2)}{recIdx}.blockdata(eventIdx).spiketimes(cellIdx);
                        
                        %save info for cells active during this epoch
                        if ~isempty(tmpspiketimes)
                            activecount = activecount+1;
                            activespiketimes{activecount} = tmpspiketimes;
                            activerates = [activerates; placefields_1msblocks{dayindex(i,1)}{dayindex(i,2)}.ratenormocc(matches(trainingcell),:)];
                        end
                    end
                %end
            end
            trainingdata = trainingdata*binsize; %transform rates to expected number of spikes, is this right?
            activerates = activerates*binsize;
            
            %decode the data
            if ~isempty(activespiketimes)
                if max(cellfun(@(x) length(x), activespiketimes)) >= minspikes %DOUBLE CHECK THIS, JUST SAYING AT LEAST 2 SPIKES IN THAT BLOCK NOT NECESSARILY 2 UNITS
                    %active spike times is the data to test during the 1ms block
                    %active rates is the training data used to decode position
                    decodingOutput{recIdx}{eventIdx} = decodePositionDurTheta(activespiketimes, activerates, timebins{eventIdx});
                    positionEstimate{recIdx}{eventIdx} = sum(decodingOutput{recIdx}{eventIdx}.spatialprob,2);
                else
                    decodingOutput{recIdx}{eventIdx} = [];
                    positionEstimate{recIdx}{eventIdx} = [];
                end
            end
            
            
            %get actual position data for each rec and event index
            startpos = find(positionInfo(recIdx).timeSmooth == startevent);
            endpos =  find(positionInfo(recIdx).timeSmooth == endevent);
            posActual{recIdx}{eventIdx} = positionInfo(recIdx).thetaSmooth(startpos:endpos);
        end
        
        if exist('positionEstimate{recIdx}') 
            %plot the data for all of the event indices
            plotbinsize = 0.02;
            totaltime = 0:plotbinsize:100-plotbinsize;
            posEstOverTime = nan(180,length(totaltime));
            counter = 0;
            for eventIdx = 1:size(positionEstimate{recIdx},2)%length(testingdata_1msblocks{dayindex(i,1)}{dayindex(i,2)}{recIdx}.blockdata)
                endofrec = 0;
                if totaltime(end) > positionInfo(recIdx).timeSmooth(end)
                    totaltime(end) = positionInfo(recIdx).timeSmooth(end);
                    endofrec = 1;
                end
                if (totaltime(end) > timebins{eventIdx}(end)) && endofrec ~=1
                    if ~isempty(positionEstimate{recIdx}{eventIdx})
                        blockSize = timebins{eventIdx}(1):plotbinsize:timebins{eventIdx}(2)-plotbinsize;
                        if ~isempty(blockSize)
                            posEstForBlock = repmat(positionEstimate{recIdx}{eventIdx},1,size(blockSize,2));
                            %posEstForBlock = repmat(positionEstimate{recIdx}{eventIdx},1,size(timebins{eventIdx},2));
                            timeIdx = find(ismember(round(totaltime,2),round(timebins{eventIdx},2)));
                            %posEstOverTime(:,timeIdx) = posEstForBlock;
                            posEstOverTime(:,timeIdx(1):timeIdx(2)-1) = posEstForBlock;
                        end
                    end
                elseif ((totaltime(end) < timebins{eventIdx}(end)) && (sum(sum(posEstOverTime)) ~= 0)) && endofrec ~=1
                    counter = counter+1;
                    %posEstOverTime(posEstOverTime == 0) = nan;
                    %plot the estimated position
                    figure; hold on;
                    subplot(2,1,1)
                    %imAlpha = ones(size(posEstOverTime));
                    %imAlpha(isnan(posEstOverTime)) = 0;
                    %imagesc(totaltime,posvect(1:end-1),posEstOverTime,'AlphaData',imAlpha);
                    imagesc(totaltime,posvect(1:end-1),posEstOverTime);
                    set(gca,'color',[1 1 1])
                    set(gca,'YDir','normal')
                    oldclims = caxis;
                    caxis([(oldclims(1)+(oldclims(2)-oldclims(1))/10) oldclims(2)]);
                    
                    %colormap(flipud(gray));
                    %colorbar
                    title('Decoded Position')
                    ylabel('Degrees')
                    
                    %plot the actual position
                    subplot(2,1,2)
                    startpos = find(round(positionInfo(recIdx).timeSmooth,2) == round(totaltime(1),2));
                    endpos = find(round(positionInfo(recIdx).timeSmooth,2) == round(totaltime(end),2));
                    %startpos = find(round(positionInfo(recIdx).timeSmooth) == round(totaltime(1)));
                    %endpos = find(round(positionInfo(recIdx).timeSmooth) == round(totaltime(end)));
                    posActualOverTime = positionInfo(recIdx).thetaSmooth(startpos:endpos);
                    %posActualOverTime_downsample = resample(posActualOverTime,length(totaltime),length(posActualOverTime));
                    plot(totaltime,posActualOverTime,'k','LineWidth',2)
                    title('Actual Position')
                    ylabel('Degrees')
                    xlabel('Time (s)')
                    set(gcf,'units','normalize','outerpos',[0 0 1 1])
                    suptitle(['Actual vs. Decoded Position for 100s periods - F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) ' - ' num2str(counter)])
                    animaldir = [savedfiguresdir 'F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '\'];
                    filename = [animaldir 'decodedpos_100sblocks_' num2str(counter)];
                    saveas(gcf,filename,'png'); saveas(gcf,filename,'fig');
                    
                    %reset variables
                    totaltime = totaltime + 100;
                    posEstOverTime = nan(180,length(totaltime));
                elseif endofrec == 1
                    totaltime = totaltime+100;
                    posEstOverTime = nan(180,length(totaltime));                
                end
            end
            
            %plot the estimated position vs actual position in small chunks
            posActualForBlock = []; posEstForBlock = [];
            for eventIdx = 1:size(positionEstimate{recIdx},2)
                if ~isempty(positionEstimate{recIdx}{eventIdx})
                    blockSize = timebins{eventIdx}(1):plotbinsize:timebins{eventIdx}(2);
                    posEstForBlock = [posEstForBlock repmat(positionEstimate{recIdx}{eventIdx},1,size(blockSize,2))];
                    %posEstForBlock = [posEstForBlock repmat(positionEstimate{recIdx}{eventIdx},1,size(timebins{eventIdx},2))];
                    posActualForBlock = [posActualForBlock; posActual{recIdx}{eventIdx}];
                    if size(posEstForBlock,2) ~= size(posActualForBlock,1)
                        pause
                    end
                end
            end
            if ~isempty(posActualForBlock)
                [sortedPosActual sortedPosActualInd] = sort(posActualForBlock);
                for locIdx = 1:length(sortedPosActualInd)
                    sortedPosEst(:,locIdx) = posEstForBlock(:,sortedPosActualInd(locIdx));
                end
                figure; imagesc('XData',sortedPosActual,'YData',[0:2:360-2],'CData',sortedPosEst,[0 0.5])
                xlim([0 360]); ylim([0 360])
                xlabel('Actual position'); ylabel('Estimated position')
                filename = [animaldir 'decodedposvsestpos_allsorted'];
                saveas(gcf,filename,'png'); saveas(gcf,filename,'fig');
            end
                
            %plot the averaged estimated position vs actual position
            posActualForBlock = []; posEstForBlock = [];
            for eventIdx = 1:size(positionEstimate{recIdx},2)
                if ~isempty(positionEstimate{recIdx}{eventIdx})
                    blockSize = timebins{eventIdx}(1):plotbinsize:timebins{eventIdx}(2);
                    posEstForBlock = [posEstForBlock repmat(positionEstimate{recIdx}{eventIdx},1,size(blockSize,2))];
                    %posEstForBlock = [posEstForBlock repmat(positionEstimate{recIdx}{eventIdx},1,size(timebins{eventIdx},2))];
                    posActualForBlock = [posActualForBlock; posActual{recIdx}{eventIdx}];
                end
            end
            if ~isempty(posActualForBlock)
                posbins = [0:2:360];
                [poshist,posedges,poshistidx] = histcounts(posActualForBlock,posbins);
                for binidx = 1:length(poshist)
                    posEstForBin = posEstForBlock(:,find(poshistidx == binidx));
                    posEstForBinAvg(:,binidx) = nanmean(posEstForBin,2);
                end
                
                figure; imagesc('XData',[0:2:360-2],'YData',[0:2:360-2],'CData',posEstForBinAvg)
                xlim([0 360]); ylim([0 360])
                xlabel('Actual position'); ylabel('Estimated position')
                filename = [animaldir 'decodedposvsestpos_avg'];
                saveas(gcf,filename,'png'); saveas(gcf,filename,'fig');
            end
            close all
        end
    end
    clear positionEstimate decodingOutput posActual
end 

%% decode sequences during ripples
% get spike counts during ripples for all place cells

% find matching cells with training and decoding data

% use training data to get p(spikes|location) = p(spikes|location)p(location)/p(spikes)

% use test data to decode position

% calculate significance of replay events using regression

%% pairwise analysis of spike timing during theta vs. during ripples

% get peak spike counts during ripples for all place cells

% get peak spike counts during theta for all place cells

% find matching cells with training and decoding data

% compare timing between ripple peaks and theta peaks for all cell pairs within a ripple



