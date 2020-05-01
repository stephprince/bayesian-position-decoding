%% population decoding
% SP 9.12.18
% this script uses Bayesian population decoding to decode location using
% place cell spike trains

%% get index
animals = [1,7,8,9,10,11,12,13,14];

%set data directories and animals
savedfiguresdir = '\\neuro-cloud\labs\singer\Steph\Figures\decoding_5050_peakrate1_maxavgrate10_190228\';
if ~exist(savedfiguresdir); mkdir(savedfiguresdir); end;
processeddatadir = '\\neuro-cloud\labs\singer\ProcessedData\VR_AnnularTrack\';
spreadsheetdir = '\\neuro-cloud\labs\singer\Steph\Code\spreadsheets\VRAnnularTrackSpreadsheetNoNovel.xlsx';
celltypedir = '\\neuro-cloud\labs\singer\Steph\Figures\celltype_newstruct\celltypeProps\';

%all rec index
[allindex iden] = getdefaultindex5XFAD(processeddatadir, animals);
allindex2 = allindex(ismember(allindex(:,1),animals),[1:2 5]);
allrecidx = allindex(ismember(allindex(:,1),animals),:);
dayindex = unique(allindex2,'rows');

%behavior index (only recordings when running through virtual environment)
behavioridx = getBehaviorRecIndex(animals, spreadsheetdir);
allbehavioridx = allindex(ismember(allindex(:,[1:3 5]),behavioridx,'rows'),:);

%set params
params.stdthreshold = 0;
params.peakratethreshold = 1;
params.maxmeanratethreshold = 10;
params.reliabilitythreshold = 0; %0.01 is what it was before
params.spatialinfothreshold = 0;
params.traintesttype = '5050';
params.decodebinsize = 0.5; %sec
params.PFbinsize = 2; %degrees

%% load place field info
for i = 1:size(dayindex) 
    sessindex = allrecidx(ismember(allrecidx(:,2), dayindex(i,2)),[1:3,5]);
    placefields{dayindex(i,1)}{dayindex(i,2)} = getplacecells(sessindex, behavioridx, celltypedir, savedfiguresdir, params);
end

%% get training data set (spike rate x space)
for i = 1:size(dayindex,1)
    sessindex = allrecidx(ismember(allrecidx(:,2), dayindex(i,2)),[1:3,5]);
    [trainingtimes{dayindex(i,1)}{dayindex(i,2)} testingtimes{dayindex(i,1)}{dayindex(i,2)} trainingdata_exclcriteria{dayindex(i,1)}{dayindex(i,2)}]...
        = gettrainingdata(sessindex, dayindex(i,:), behavioridx, placefields{dayindex(i,1)}{dayindex(i,2)}, savedfiguresdir, processeddatadir, params);
end

%% get testing/decoding data set (spike counts x temporal bins)
for i = 1:size(dayindex,1)
    sessindex = allrecidx(ismember(allrecidx(:,2), dayindex(i,2)),[1:3,5]);
    testingdata{dayindex(i,1)}{dayindex(i,2)}...
        = gettestingdata(sessindex,dayindex(i,:),behavioridx,testingtimes{dayindex(i,1)}{dayindex(i,2)},savedfiguresdir,processeddatadir);
end

%% decode position during theta
for i = 1:size(dayindex,1)
   sessindex = allrecidx(ismember(allrecidx(:,2), dayindex(i,2)),[1:3,5]);
   [posEstAvg{dayindex(i,1)}{dayindex(i,2)}] = getdecodingdurtheta(sessindex,behavioridx,trainingdata_exclcriteria{dayindex(i,1)}{dayindex(i,2)},testingdata{dayindex(i,1)}{dayindex(i,2)},savedfiguresdir, processeddatadir, params);
end 

%% plot decoded positions for all animals
plotdecodingdurtheta(dayindex,posEstAvg,trainingdata_exclcriteria,savedfiguresdir,params);

%% plot decoding errors for individual animals
figure(100); clf; hold on;
for i = 1:size(dayindex,1)
    %concatenate all recording files for a single session
    posEstAvgtemp = [];
    for recIdx = 1:numel(posEstAvg{dayindex(i,1)}{dayindex(i,2)})
        if ~isnan(posEstAvg{dayindex(i,1)}{dayindex(i,2)}{recIdx})
            posEstAvgtemp = cat(3, posEstAvgtemp, posEstAvg{dayindex(i,1)}{dayindex(i,2)}{recIdx});
        end
    end
    posEstAvg_all = nanmean(posEstAvgtemp,3);
    
    if ~isempty(posEstAvg_all)
        %calculate error
        actualpos = [1:params.PFbinsize:360];
        [posEstPeaki posEstPeakj] = find(posEstAvg_all == max(posEstAvg_all));
        estpos = posEstPeaki'*params.PFbinsize;
        
        %get RMSerror for the whole session
        for binIdx = 1:length(actualpos)
            residsquared(binIdx) = (estpos(binIdx)-actualpos(binIdx))^2;
        end
        rmsError(i) = sqrt(sum(residsquared)/length(actualpos));
        
        %get cumulative distribution function of errors (each bin)
        edges = 0:10:360;
        errorhist = histcounts(estpos-actualpos,edges);
        errorhistnorm(i,:) = errorhist/sum(errorhist);
        plot(edges(1:end-1),cumsum(errorhistnorm(i,:)));
        xlabel('Error (degrees)'); ylabel('Cumulative fraction');
        title('Decoding error')
    end
end
filename = [savedfiguresdir 'decodingerror_all'];
saveas(gcf,filename,'png'); saveas(gcf,filename,'fig');
