%% population decoding
% SP 9.12.18
% this script uses Bayesian population decoding to decode location using
% place cell spike trains

%% get index
animals = [1,7,8 9];

%set data directories and animals
savedfiguresdir = '\\neuro-cloud\labs\singer\Steph\Figures\decoding_reliability_v2\';
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

%set threshold flags
stdthreshold = 0;
peakratethreshold = 0;
reliabilitythreshold = 0.1;

%% load place field info
for i = 1:size(dayindex) 
    if ~exist(savedfiguresdir); mkdir(savedfiguresdir); end;
    datadir = [savedfiguresdir 'F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '\'];
    filename = [datadir 'placefields_F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '.mat'];
    
    % get linearized spike rates x time normalized for occupancy
    spikeratenormocc{dayindex(i,1)}{dayindex(i,2)} = calcspikeratenormocc(processeddatadir, dayindex(i,:), behavioridx);
        
    %apply exclusion criteria
    placefields{dayindex(i,1)}{dayindex(i,2)} = applyplacefieldcriteria(spikeratenormocc{dayindex(i,1)}{dayindex(i,2)}, dayindex(i,:), behavioridx, processeddatadir, ...
        'pyr', celltypedir, 'std', stdthreshold, 'peak', peakratethreshold, 'reliability_alltimes', reliabilitythreshold)
    placefields_stabletimes{dayindex(i,1)}{dayindex(i,2)} = applyplacefieldcriteria(spikeratenormocc{dayindex(i,1)}{dayindex(i,2)}, dayindex(i,:), behavioridx, processeddatadir, ...
        'pyr', celltypedir, 'std', stdthreshold, 'peak', peakratethreshold, 'reliability_stabletimes', reliabilitythreshold);
    
    %plot placefields to see span over track
    plotplacefields_sortedpeakfr(placefields{dayindex(i,1)}{dayindex(i,2)},dayindex(i,:),savedfiguresdir);

    %save the data
    filename = [datadir 'placefields_F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '.mat'];
    if ~exist(datadir); mkdir(datadir); end;
    save(filename, 'placefields');
end

%% get training data set (spike rate x space)
for i = 1:size(dayindex,1)
    % get times of training and testing split
    spikedatadir = [processeddatadir 'F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '\Clusters\Sorted\'];
    recs = behavioridx(ismember(behavioridx(:,2),dayindex(i,2)),3);
    [trainingtimes{dayindex(i,1)}{dayindex(i,2)} testingtimes{dayindex(i,1)}{dayindex(i,2)}] = splittraintestdata(spikedatadir, recs, '5050');
    
    %generate training data place fields
    trainingdata{dayindex(i,1)}{dayindex(i,2)} = calcspikeratenormocc_trainingdata(processeddatadir, dayindex(i,:), behavioridx, trainingtimes{dayindex(i,1)}{dayindex(i,2)});
    
    %apply exclusion criteria from orignal place fields
    excl = placefields{dayindex(i,1)}{dayindex(i,2)}.exclallcriteriaidx;
    fnames = fieldnames( trainingdata{dayindex(i,1)}{dayindex(i,2)});
    for fieldIdx = 1:3
        temp = trainingdata{dayindex(i,1)}{dayindex(i,2)}.(fnames{fieldIdx});
        temp(excl,:) = nan; %replace excl indices with nans
        trainingdata_exclcriteria{dayindex(i,1)}{dayindex(i,2)}.(fnames{fieldIdx}) = temp;
    end
    
    %plot training data set place fields
    savedfiguresdir2 = [savedfiguresdir 'placefields_trainingdata\'];
    if ~exist(savedfiguresdir2); mkdir(savedfiguresdir2); end
    plotplacefields_sortedpeakfr(trainingdata_exclcriteria{dayindex(i,1)}{dayindex(i,2)}, savedfiguresdir2);
    
    %save the data
    datadir = [savedfiguresdir 'F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '\'];
    filename = [datadir 'trainingdata_F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '.mat'];
    if ~exist(datadir); mkdir(datadir); end;
    save(filename, 'trainingdata_exclcriteria');
end

%% get testing/decoding data set (spike counts x temporal bins)
for i = 1:size(dayindex,1)
    % get spikecounts over temporal bins to be decoded
    testingdata{dayindex(i,1)}{dayindex(i,2)} = calcspikestodecode_testingdata(processeddatadir, dayindex(i,:), behavioridx, testingtimes{dayindex(i,1)}{dayindex(i,2)});
    
    %save the data
    datadir = [savedfiguresdir 'F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '\'];
    filename = [datadir 'testinggdata_F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '.mat'];
    if ~exist(datadir); mkdir(datadir); end;
    save(filename, 'testingdata');
end

%% decode position during theta
minspikes = 2; %minimum number of spikes in order to attempt to decode block
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
        %find matching cells and decode position for each block
        [posActual{recIdx} positionEstimate{recIdx} eventTimes] = calcPositionEstimate(trainingdata_exclcriteria,testingdata,positionInfo);        
        
        % plot the results
        if exist('positionEstimate{recIdx}') 
            %plot the data in 100s blocks to look at decoder on trial by trial basis
            plotEstActPosition_100sblocks(positionEstimate, positionInfo, eventTimes, dayindex, savedfiguresdir)
                
            %plot the averaged estimated position vs actual position
            plotEstActPosition_avg(positionEstimate, posActual, eventTimes, savedfiguresdir)
            close all
        end
    end
    clear positionEstimate posActual
end 