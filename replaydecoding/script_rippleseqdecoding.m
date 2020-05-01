%% rippleseq decoding
% SP 2.4.19
% this script uses Bayesian population decoding to decode ripple seqs using
% place cell spike trains

%% get index
animals = [1,7,8,9,10,11,12,13,14];
dates = [];

%set data directories and animals
savedfiguresdir = '\\neuro-cloud\labs\singer\Steph\Figures\rippleseqdecoding_minrate1_peakrate10_mincell3_nonthetanew_190226\';
if ~exist(savedfiguresdir); mkdir(savedfiguresdir); end;
processeddatadir = '\\neuro-cloud\labs\singer\ProcessedData\VR_AnnularTrack\';
spreadsheetdir = '\\neuro-cloud\labs\singer\Steph\Code\spreadsheets\VRAnnularTrackSpreadsheetNoNovel.xlsx';
celltypedir = '\\neuro-cloud\labs\singer\Steph\Figures\celltype_newstruct\celltypeProps\';

%all rec index
[allindex iden] = getdefaultindex5XFAD(processeddatadir, animals);
if ~isempty(dates) 
    allindex = allindex(ismember(allindex(:,2),dates),:); %exclude by date
end
allindex2 = allindex(ismember(allindex(:,1),animals),[1:2 5]);
allrecidx = allindex(ismember(allindex(:,1),animals),:);
dayindex = unique(allindex2,'rows');

%behavior index (only recordings when running through virtual environment)
behavioridx = getBehaviorRecIndex(animals, spreadsheetdir);
allbehavioridx = allindex(ismember(allindex(:,[1:3 5]),behavioridx,'rows'),:);

%set threshold flags
stdthreshold = 0;
peakratethreshold = 1;
maxmeanratethreshold = 10;
reliabilitythreshold = 0; %0.01 is what it was before
spatialinfothreshold = 0;
PFbinsize = 2; %degrees
mincells = 3;

%% load place field info
for i = 1:size(dayindex,1) 
    sessindex = allrecidx(ismember(allrecidx(:,2), dayindex(i,2)),[1:3,5]);
    session = dayindex(i,:);
    datadir = [savedfiguresdir 'F' num2str(session(1)) '_' num2str(session(2)) '\'];
    if ~exist(datadir,'dir'); mkdir(datadir); end;
    filename = [datadir 'placefields_F' num2str(session(1)) '_' num2str(session(2)) '.mat'];
    
    % get linearized spike rates x time normalized for occupancy
    filename = [datadir 'spikeratenormocc_F' num2str(session(1)) '_' num2str(session(2)) '.mat'];
    if exist(filename,'file')
        load(filename);
    else
        spikeratenormocc{session(1)}{session(2)} = calcspikeratenormocc(processeddatadir, session, behavioridx, sessindex,PFbinsize);
        save(filename,'spikeratenormocc')
    end
        
    %apply exclusion criteria
    disp(['Applying place field criteria for F' num2str(session(1)) num2str(session(2))])
    placefields{dayindex(i,1)}{dayindex(i,2)} = applyplacefieldcriteria(spikeratenormocc{session(1)}{session(2)}, session, behavioridx, sessindex, processeddatadir,PFbinsize, datadir,...
        'pyr', celltypedir, 'std', stdthreshold, 'peak', peakratethreshold, 'maxmean',maxmeanratethreshold,'spatialinfo',spatialinfothreshold);
    
    %plot placefields to see span over track
    plotplacefields_sortedpeakfr(placefields{dayindex(i,1)}{dayindex(i,2)},session,savedfiguresdir);
  
    %save the data
    filename = [datadir 'placefields_F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '.mat'];
    if ~exist(datadir); mkdir(datadir); end;
    save(filename, 'placefields');
    close all;
end

%% get ripple event data (spike counts x temporal bins)
for i = 1:size(dayindex,1)
    sessindex = allrecidx(ismember(allrecidx(:,2), dayindex(i,2)),[1:3,5]);
    datadir = [savedfiguresdir 'F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '\'];
    if ~exist(datadir); mkdir(datadir); end;
    filename = [datadir 'rippleevents_F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '.mat'];
    
    %if exist(filename)
    %    load(filename)
    %    rippleeventstemp{dayindex(i,1)}{dayindex(i,2)} = rippleevents{dayindex(i,1)}{dayindex(i,2)};
    %else
        %spikes during ripples
        disp(['Getting spikes during ripple events for F' num2str(dayindex(i,1)) num2str(dayindex(i,2))])
        rippleevents{dayindex(i,1)}{dayindex(i,2)} = calcspikesdurrippleevents(processeddatadir, mincells, dayindex(i,:),allindex);
        save(filename,'rippleevents')
        rippleeventstemp{dayindex(i,1)}{dayindex(i,2)} = rippleevents{dayindex(i,1)}{dayindex(i,2)};
    %end   
end
clear rippleevents thetaevents
rippleevents = rippleeventstemp;

%% get ripple sequences 
for i = 1:size(dayindex,1)
    %get replay fidelity struct
    datadir = [savedfiguresdir 'F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '\'];
    filename = [datadir 'rippleseqdata_' num2str(dayindex(i,1)) num2str(dayindex(i,2)) '.mat'];
%     if exist(filename)
%         load(filename);
%         rippleseqdatatemp{dayindex(i,1)}{dayindex(i,2)} = rippleseqdata{dayindex(i,1)}{dayindex(i,2)};
%     else
        disp(['Identifying ripple sequences for F' num2str(dayindex(i,1)) num2str(dayindex(i,2))])
        if ~isempty(placefields{dayindex(i,1)}{dayindex(i,2)})
            rippleseqdata{dayindex(i,1)}{dayindex(i,2)} = getRippleSequences(placefields{dayindex(i,1)}{dayindex(i,2)},rippleevents{dayindex(i,1)}{dayindex(i,2)},mincells,allindex,dayindex(i,:),processeddatadir);
        else
            rippleseqdata{dayindex(i,1)}{dayindex(i,2)} = [];
        end
        save(filename,'rippleseqdata');
        rippleseqdatatemp{dayindex(i,1)}{dayindex(i,2)} = rippleseqdata{dayindex(i,1)}{dayindex(i,2)};
%     end
end
clear rippleseqdata
rippleseqdata = rippleseqdatatemp;

%% concatenate data across recordings and plot ripple events
for i = 1:size(dayindex,1)
    datadir = [savedfiguresdir 'F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '\'];
    
    %get significant ripple events
    sigReplayIdx = [];  sigReplayCuesIdx = [];
    rvalthreshold = 0;
    if ~isempty(rippleseqdata{dayindex(i,1)}{dayindex(i,2)})
        for eventIdx = 1:length(rippleseqdata{dayindex(i,1)}{dayindex(i,2)}.replayStats)
            %for normal replay
            if ~isempty(rippleseqdata{dayindex(i,1)}{dayindex(i,2)}.replayStats{eventIdx}) ...
                    && ((rippleseqdata{dayindex(i,1)}{dayindex(i,2)}.replayStats{eventIdx}.pval < 0.05)...
                    || rippleseqdata{dayindex(i,1)}{dayindex(i,2)}.replayStats{eventIdx}.rvalprctile > rvalthreshold)
                sigReplayIdx = [sigReplayIdx eventIdx];
            end
            
            %for cues only replay
            if ~isempty(rippleseqdata{dayindex(i,1)}{dayindex(i,2)}.replayStatsCues{eventIdx}) ...
                    && ((rippleseqdata{dayindex(i,1)}{dayindex(i,2)}.replayStatsCues{eventIdx}.pval < 0.05)...
                    || rippleseqdata{dayindex(i,1)}{dayindex(i,2)}.replayStatsCues{eventIdx}.rvalprctile > rvalthreshold)
                sigReplayCuesIdx = [sigReplayCuesIdx eventIdx];
            end
        end
    end
    
    %plot PDF of significant ripple events
%     for eventIdx = 1:length(sigReplayIdx)
%         pdf = rippleseqdata{dayindex(i,1)}{dayindex(i,2)}.decodingOutput{sigReplayIdx(eventIdx)}.spatialprob;
%         cmap = [1 0 0; 1 0.5 0; 0 1 0; 0 0.5 0.5; 0 0 1; 0.25 0 1; 0.5 0 1];
%         %cmap = jet(size(pdf,2)*1); cmap = flipud(cmap);
%         colorset = cmap([1:size(pdf,2)],:);
%         figure; hold on; 
%         for m = 1:size(pdf,2)
%             plot(pdf(:,m),'Color',colorset(m,:),'LineWidth',2);
%             legend;
%         end
%     end
    
    %plot heat maps for significant ripple events
    for eventIdx = 1:length(sigReplayIdx)
        figure('units','normalized','outerposition',[0 0 1 1]); hold on;
        pdf = rippleseqdata{dayindex(i,1)}{dayindex(i,2)}.decodingOutput{sigReplayIdx(eventIdx)}.spatialprob;
        totalspikes = sum(rippleseqdata{dayindex(i,1)}{dayindex(i,2)}.decodingOutput{sigReplayIdx(eventIdx)}.totalspikes);
        totalcells = sum(rippleseqdata{dayindex(i,1)}{dayindex(i,2)}.decodingOutput{sigReplayIdx(eventIdx)}.totalcellsactive);
        rprctile = rippleseqdata{dayindex(i,1)}{dayindex(i,2)}.replayStats{sigReplayIdx(eventIdx)}.rvalprctile;
        imagesc([0:0.015:0.015*size(pdf,2)],[0:2:360],pdf,[0 0.03])
        colormap(hot); colorbar;
        xlim([0 0.015*size(pdf,2)]); ylim([0 360]);
        xlabel('Time from start of ripple (s)'); ylabel('Decoded position')
        title(['Replay seq F' num2str(dayindex(i,1)) num2str(dayindex(i,2)) ' ripple ' num2str(eventIdx) ' rprctile: ' num2str(rprctile)  'spikes: ' num2str(totalspikes) 'cells: ' num2str(totalcells)])
        datadir = [savedfiguresdir 'F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '\'];
        filename = [datadir 'rippleseqplots_' num2str(dayindex(i,1)) num2str(dayindex(i,2)) num2str(eventIdx)];
        saveas(gcf,filename,'png'); saveas(gcf,filename,'fig');
    end
    close all;
    
    %plot heat maps for significant ripple events
    for eventIdx = 1:length(sigReplayCuesIdx)
        figure('units','normalized','outerposition',[0 0 1 1]); hold on;
        pdf = rippleseqdata{dayindex(i,1)}{dayindex(i,2)}.decodingOutputCues{sigReplayCuesIdx(eventIdx)}.spatialprob;
        totalspikes = sum(rippleseqdata{dayindex(i,1)}{dayindex(i,2)}.decodingOutputCues{sigReplayCuesIdx(eventIdx)}.totalspikes);
        totalcells = sum(rippleseqdata{dayindex(i,1)}{dayindex(i,2)}.decodingOutputCues{sigReplayCuesIdx(eventIdx)}.totalcellsactive);
        rprctile = rippleseqdata{dayindex(i,1)}{dayindex(i,2)}.replayStatsCues{sigReplayCuesIdx(eventIdx)}.rvalprctile;
        imagesc([0:0.015:0.015*size(pdf,2)],[0:2:180],pdf,[0 0.03])
        colormap(hot); colorbar;
        xlim([0 0.015*size(pdf,2)]); ylim([0 180]);
        xlabel('Time from start of ripple (s)'); ylabel('Decoded cue position')
        title(['Replay seq cues only F' num2str(dayindex(i,1)) num2str(dayindex(i,2)) ' ripple ' num2str(eventIdx) ' rprctile: ' num2str(rprctile)  'spikes: ' num2str(totalspikes) 'cells: ' num2str(totalcells)])
        datadir = [savedfiguresdir 'F' num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '\'];
        filename = [datadir 'rippleseqcuesonlyplots_' num2str(dayindex(i,1)) num2str(dayindex(i,2)) num2str(eventIdx)];
        saveas(gcf,filename,'png'); saveas(gcf,filename,'fig');
    end
    close all;
    
    filename = [datadir 'rippleseqdata_' num2str(dayindex(i,1)) num2str(dayindex(i,2)) num2str(eventIdx)];
    save(filename,'rippleseqdata');
end

%% compare slopes and r values between 5XFAD and WT groups
group = [1,7,10,11,13 ; 8,9,12,14,nan]; clr = 'gk';

for g = [2 1]
    inclsess = find(ismember(dayindex(:,1), group(g,:)));
    rvalgroup{g} = []; pvalgroup{g} = []; rvalprctilegroup{g} = []; slopegroup{g} = [];
    for i = inclsess'
        sessindex = dayindex(i,:);
        %loop through events to concatenate metrics
        if ~isempty(rippleseqdata{sessindex(1)}{sessindex(2)})
            for eventIdx = 1:length(rippleseqdata{sessindex(1)}{sessindex(2)}.replayStats)
                if ~isempty(rippleseqdata{sessindex(1)}{sessindex(2)}.replayStats{eventIdx})
                    rvalgroup{g} = [rvalgroup{g} rippleseqdata{sessindex(1)}{sessindex(2)}.replayStats{eventIdx}.rval];
                    pvalgroup{g} = [pvalgroup{g} rippleseqdata{sessindex(1)}{sessindex(2)}.replayStats{eventIdx}.pval];
                    rvalprctilegroup{g} = [rvalprctilegroup{g} rippleseqdata{sessindex(1)}{sessindex(2)}.replayStats{eventIdx}.rvalprctile];
                    slopegroup{g} = [slopegroup{g} rippleseqdata{sessindex(1)}{sessindex(2)}.replayStats{eventIdx}.slope];
                end
            end
        end
    end
end

%plot 5XFAD vs WT rval percentiles
figure; hold on;
edges = 0:10:100;
for g = [2 1]
    rvalprctilehist = histcounts(rvalprctilegroup{g},edges);
    rvalprctilehistnorm = rvalprctilehist/sum(rvalprctilehist);
    plot(edges(2:end),rvalprctilehistnorm,clr(g));
    xlabel('r-value percentile')
    ylabel('fraction of decoded events')
end
[h p] = kstest2(rvalprctilegroup{1},rvalprctilegroup{2});
title(['Ripple seq rval percentiles 5XFAD vs. WT kstest2 - ' num2str(p)])
filename = [savedfiguresdir 'distrvalpercentilesnorm_5XFADvsWT'];
saveas(gcf,filename,'fig'); saveas(gcf,filename,'png')

%plot 5XFAD vs WT rval percentiles
figure; hold on;
edges = 0:10:100;
for g = [2 1]
    rvalprctilehist = histcounts(rvalprctilegroup{g},edges);
    plot(edges(2:end),rvalprctilehist,clr(g));
    xlabel('r-value percentile')
    ylabel('number of decoded events')
end
[h p] = kstest2(rvalprctilegroup{1},rvalprctilegroup{2});
title(['Ripple seq rval percentiles 5XFAD vs. WT kstest2 - ' num2str(p)])
filename = [savedfiguresdir 'distrvalpercentiles_5XFADvsWT'];
saveas(gcf,filename,'fig'); saveas(gcf,filename,'png')

%plot 5XFAD vs WT pvals
figure; hold on;
edges = 0:0.05:1;
for g = [2 1]
    pvalhist = histcounts(pvalgroup{g},edges);
    pvalhistnorm = pvalhist/sum(pvalhist);
    plot(edges(1:end-1),pvalhistnorm,clr(g));
    xlabel('p-value')
    ylabel('fraction of decoded events')
end
[h p] = kstest2(pvalgroup{1},pvalgroup{2});
title(['Ripple seq pvals 5XFAD vs. WT kstest2 - ' num2str(p)])
filename = [savedfiguresdir 'distpvalpercentilesnorm_5XFADvsWT'];
saveas(gcf,filename,'fig'); saveas(gcf,filename,'png')

%plot 5XFAD vs WT rvals
figure; hold on;
edges = 0:0.01:0.5;
for g = [2 1]
    rvalgroup{g}(rvalgroup{g} == nan) = [];
    rvalhist = histcounts(rvalgroup{g},edges);
    rvalhistnorm = rvalhist/sum(rvalhist);
    plot(edges(2:end),cumsum(rvalhistnorm),clr(g));
    xlabel('r-value')
    ylabel('fraction of decoded events')
end
[h p] = kstest2(rvalgroup{1},rvalgroup{2});
title(['Ripple seq rvals 5XFAD vs. WT kstest2 - ' num2str(p)])
filename = [savedfiguresdir 'distrvalsnorm_5XFADvsWT'];
saveas(gcf,filename,'fig'); saveas(gcf,filename,'png')

%plot 5XFAD vs WT pval slopes
figure; hold on;
edges = -60:5:130;
for g = [2 1]
    slopehist = histcounts(slopegroup{g},edges);
    slopehistnorm = slopehist/sum(slopehist);
    plot(edges(2:end),slopehistnorm,clr(g));
    xlabel('slope')
    ylabel('fraction of decoded events')
end
[h p] = kstest2(slopegroup{1},slopegroup{2});
title(['Ripple seq slopes 5XFAD vs. WT kstest2 - ' num2str(p)])
filename = [savedfiguresdir 'distslopesnorm_5XFADvsWT'];
saveas(gcf,filename,'fig'); saveas(gcf,filename,'png')

%plot 5XFAD vs WT slopes cumulative fraction
figure; hold on;
edges = -60:5:130;
for g = [2 1]
    slopehist = histcounts(slopegroup{g},edges);
    slopehistnorm = slopehist/sum(slopehist);
    plot(edges(2:end),cumsum(slopehistnorm),clr(g));
    xlabel('slope')
    ylabel('fraction of decoded events')
end
[h p] = kstest2(slopegroup{1},slopegroup{2});
title(['Ripple seq slopes 5XFAD vs. WT kstest2 - ' num2str(p)])
filename = [savedfiguresdir 'distslopescumfract_5XFADvsWT'];
saveas(gcf,filename,'fig'); saveas(gcf,filename,'png')