%% population decoding
% SP 9.12.18
% this script uses Bayesian population decoding to decode location using
% place cell spike trains

%% get index
animals = [1,7,8,9,10,11,12,13,14];

%set data directories and animals
savedfiguresdir = '\\neuro-cloud\labs\singer\Steph\Figures\decoding_spatialinfo_5050_bins4_peakrate1_maxavgrate10_190123\';
if ~exist(savedfiguresdir); mkdir(savedfiguresdir); end;
processeddatadir = '\\neuro-cloud\labs\singer\ProcessedData\VR_AnnularTrack\';
spreadsheetdir = '\\neuro-cloud\labs\singer\Steph\Code\spreadsheets\VRAnnularTrackSpreadsheet.xlsx';
celltypedir = '\\neuro-cloud\labs\singer\Steph\Figures\celltype_newstruct\celltypeProps\';

%all rec index
[allindex iden] = getdefaultindex5XFAD(processeddatadir, animals);
allindex2 = allindex(ismember(allindex(:,1),animals),[1:2 5]);
allrecidx = allindex(ismember(allindex(:,1),animals),:);
dayindex = unique(allindex2,'rows');

%behavior index (only recordings when running through virtual environment)
behavioridx = getBehaviorRecIndex(animals, spreadsheetdir);
allbehavioridx = allindex(ismember(allindex(:,[1:3 5]),behavioridx,'rows'),:);

%set threshold flags
stdthreshold = 0;
peakratethreshold = 1;
maxmeanratethreshold = 9;
reliabilitythreshold = 0.01; %0.01 is what it was before
spatialinfothreshold = 95;
traintesttype = '5050large';
decodebinsize = 2; %sec, don't think it's actually doing this
PFbinsize = 4; %degrees

%% load place field info
for i = 1:size(dayindex) 
    sessindex = allrecidx(ismember(allrecidx(:,2), dayindex(i,2)),[1:3,5]);
    datadir = [savedfiguresdir 'F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '\'];
    if ~exist(datadir); mkdir(datadir); end;
    filename = [datadir 'placefields_F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '.mat'];
    
    % get linearized spike rates x time normalized for occupancy
    filename = [datadir 'spikeratenormocc_F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '.mat'];
    if exist(filename)
        load(filename)
    else
        spikeratenormocc{dayindex(i,1)}{dayindex(i,2)} = calcspikeratenormocc(processeddatadir, dayindex(i,:), behavioridx, sessindex,PFbinsize);
        save(filename,'spikeratenormocc')
    end
        
    %apply exclusion criteria
    placefields{dayindex(i,1)}{dayindex(i,2)} = applyplacefieldcriteria(spikeratenormocc{dayindex(i,1)}{dayindex(i,2)}, dayindex(i,:), behavioridx, sessindex, processeddatadir,PFbinsize, datadir,...
        'pyr', celltypedir, 'std', stdthreshold, 'peak', peakratethreshold, 'maxmean',maxmeanratethreshold,'spatialinfo',spatialinfothreshold);
    
    %plot placefields to see span over track
    plotplacefields_sortedpeakfr(placefields{dayindex(i,1)}{dayindex(i,2)},dayindex(i,:),savedfiguresdir);
    
    %save the data
    filename = [datadir 'placefields_F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '.mat'];
    if ~exist(datadir); mkdir(datadir); end;
    save(filename, 'placefields');
end

%% plot place fields on circular track
% for i = 1:size(dayindex) 
%     sessindex = allrecidx(ismember(allrecidx(:,2), dayindex(i,2)),[1:3,5]);
%     datadir = [savedfiguresdir 'F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '\'];
%     if ~exist(datadir); mkdir(datadir); end;
%     filename = [datadir 'placefields_F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '.mat'];
%     
%     % get linearized spike rates x time normalized for occupancy
%     filename = [datadir 'spikeratenormocc_F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '.mat'];
%     if exist(filename)
%         load(filename)
%     else
%         spikeratenormocc{dayindex(i,1)}{dayindex(i,2)} = calcspikeratenormocc(processeddatadir, dayindex(i,:), behavioridx, sessindex);
%         save(filename,'spikeratenormocc')
%     end
%         
%     %apply exclusion criteria
%     placefields{dayindex(i,1)}{dayindex(i,2)} = applyplacefieldcriteria(spikeratenormocc{dayindex(i,1)}{dayindex(i,2)}, dayindex(i,:), behavioridx, sessindex, processeddatadir, datadir, ...
%         'pyr', celltypedir, 'std', stdthreshold, 'peak', peakratethreshold, 'reliability_stabletimes', reliabilitythreshold);
%     
%     %plot placefields to see span over track
%     plotplacefields_circular(placefields{dayindex(i,1)}{dayindex(i,2)},dayindex(i,:),savedfiguresdir, PFbinsize);
% end

%% get training data set (spike rate x space)
for i = 1:size(dayindex,1)
    sessindex = allrecidx(ismember(allrecidx(:,2), dayindex(i,2)),[1:3,5]);
    % get times of training and testing split
    spikedatadir = [processeddatadir 'F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '\Clusters\Sorted\'];
    recs = behavioridx(ismember(behavioridx(:,2),dayindex(i,2)),3);
    [trainingtimes{dayindex(i,1)}{dayindex(i,2)} testingtimes{dayindex(i,1)}{dayindex(i,2)}] = splittraintestdata(spikedatadir, recs, traintesttype);
    
    %generate training data place fields
    datadir = [savedfiguresdir 'F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '\'];
    if ~exist(datadir); mkdir(datadir); end
    filename = [datadir 'rawtrainingdata_F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '.mat'];
    if exist(filename)
        load(filename)
    else
        trainingdata{dayindex(i,1)}{dayindex(i,2)} = calcspikeratenormocc_trainingdata_stabletimes(processeddatadir, dayindex(i,:), behavioridx, sessindex, trainingtimes{dayindex(i,1)}{dayindex(i,2)},PFbinsize);
        save(filename,'trainingdata')
    end
    
    %apply exclusion criteria from orignal place fields
    if ~isempty(placefields{dayindex(i,1)}{dayindex(i,2)})
        excl = placefields{dayindex(i,1)}{dayindex(i,2)}.exclallcriteriaidx;
        fnames = fieldnames(trainingdata{dayindex(i,1)}{dayindex(i,2)});
        for fieldIdx = 1:3
            temp = trainingdata{dayindex(i,1)}{dayindex(i,2)}.(fnames{fieldIdx});
            temp(excl,:) = nan; %replace excl indices with nans
            trainingdata_exclcriteria{dayindex(i,1)}{dayindex(i,2)}.(fnames{fieldIdx}) = temp;
        end
        trainingdata_exclcriteria{dayindex(i,1)}{dayindex(i,2)}.(fnames{4}) = trainingdata{dayindex(i,1)}{dayindex(i,2)}.(fnames{4});
    else
        trainingdata_exclcriteria{dayindex(i,1)}{dayindex(i,2)} = [];
    end
    
    %plot training data set place fields
    savedfiguresdir2 = [savedfiguresdir 'placefields_trainingdata\'];
    if ~exist(savedfiguresdir2); mkdir(savedfiguresdir2); end
    plotplacefields_sortedpeakfr(trainingdata_exclcriteria{dayindex(i,1)}{dayindex(i,2)}, dayindex(i,:), savedfiguresdir2);
    
    %save the data
    filename = [datadir 'trainingdata_F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '.mat'];
    if ~exist(datadir); mkdir(datadir); end;
    save(filename, 'trainingdata_exclcriteria');
end

%% get testing/decoding data set (spike counts x temporal bins)
for i = 1:size(dayindex,1)
    % get spikecounts over temporal bins to be decoded
    datadir = [savedfiguresdir 'F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '\'];
    if ~exist(datadir); mkdir(datadir); end;
    filename = [datadir 'testinggdata_F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '.mat'];
    if exist(filename)
        load(filename)
    else
        testingdata{dayindex(i,1)}{dayindex(i,2)} = calcspikestodecode_testingdata(processeddatadir, dayindex(i,:), behavioridx, testingtimes{dayindex(i,1)}{dayindex(i,2)});
        save(filename,'testingdata')
    end
end

%% decode position during theta
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
        [posActual positionEstimate eventTimes] = calcPositionEstimate(trainingdata_exclcriteria,testingdata,positionInfo,dayindex(i,:),recIdx,decodebinsize);        
        
        % plot the results
        if exist('positionEstimate')
            %plot the data in 100s blocks to look at decoder on trial by trial basis
            plotEstActPosition_100sblocks(positionEstimate, positionInfo, eventTimes, dayindex(i,:), savedfiguresdir, recIdx, PFbinsize)
                
            %plot the averaged estimated position vs actual position
            posEstAvg{dayindex(i,1)}{dayindex(i,2)}{recIdx} = plotEstActPosition_avg(positionEstimate, posActual, eventTimes, savedfiguresdir, dayindex(i,:), recIdx,PFbinsize)
            close all
        else
            posEstAvg{dayindex(i,1)}{dayindex(i,2)}{recIdx} = nan;
        end
    end
    
    %concatenate data across recording files
    posEstAvgtemp = [];
    for recIdx = 1:numel(recs)
        if ~isnan(posEstAvg{dayindex(i,1)}{dayindex(i,2)}{recIdx})
            posEstAvgtemp = cat(3, posEstAvgtemp, posEstAvg{dayindex(i,1)}{dayindex(i,2)}{recIdx});
        end
    end
    posEstAvg_day = nanmean(posEstAvgtemp,3);
    
    %plot the average for that day
    figure; imagesc('XData',[0:PFbinsize:360-PFbinsize],'YData',[0:PFbinsize:360-PFbinsize],'CData',posEstAvg_day)
    xlim([0 360]); ylim([0 360])
    xlabel('Actual position'); ylabel('Estimated position')
    animaldir = [savedfiguresdir 'F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '\'];
    filename = [animaldir 'decodedposvsestpos_avg_wholeday'];
    saveas(gcf,filename,'png'); saveas(gcf,filename,'fig');
    
    clear positionEstimate posActual
end 

%% plot decoded positions for all animals
posEstAvgtemp = [];
for i = 1:size(dayindex,1)
    for recIdx = 1:numel(posEstAvg{dayindex(i,1)}{dayindex(i,2)})
        if ~isnan(posEstAvg{dayindex(i,1)}{dayindex(i,2)}{recIdx})
            posEstAvgtemp = cat(3, posEstAvgtemp, posEstAvg{dayindex(i,1)}{dayindex(i,2)}{recIdx});
        end
    end
end
posEstAvg_all = nanmean(posEstAvgtemp,3);
%plot the average for that day
figure; imagesc('XData',[0:PFbinsize:360-PFbinsize],'YData',[0:PFbinsize:360-PFbinsize],'CData',posEstAvg_all)
xlim([0 360]); ylim([0 360])
xlabel('Actual position'); ylabel('Estimated position')
colorbar
filename = [savedfiguresdir 'ALLdecodedposvsestpos_avg'];
saveas(gcf,filename,'png'); saveas(gcf,filename,'fig');
    