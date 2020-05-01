function out = calcspikesdurrippleevents(processeddatadir, cellcountthresh, dayindex, allindex)
%SP 11.15.18
% this function gets all the spike info during ripples

%% set data directories
spikedatadir = [processeddatadir 'F' num2str(dayindex(1)) '_' num2str(dayindex(2)) '\Clusters\Sorted\'];
recs = allindex(ismember(allindex(:,1:2), dayindex(1:2), 'rows'),3);
recDepth = dayindex(3);
fileindex = [repmat(dayindex(1:2),length(recs),1) recs];
bestRippleChan = getBestRippleChan(fileindex, processeddatadir, 'f', recDepth);
rippledatadir = [processeddatadir, 'f',num2str(dayindex(1)), '_', num2str(dayindex(2)) , '/' num2str(recDepth(1)) '/' num2str(bestRippleChan) '/'];

% get virmen and cluster info
if ~isempty(recs)
    [clusterInfo recInfo] = getClusterInfoSP(spikedatadir, recs);
else
    out = [];
    return
end

%% get ripple event times durin nontheta periods
ripples = loaddatastruct2(rippledatadir, dayindex, 'ripples', recs);
disp(['loading ripples index ', num2str(dayindex)]);
nonthetas = loaddatastruct2(rippledatadir, dayindex, 'nonthetas', recs);
disp(['loading nonthetas index ', num2str(dayindex)]);

%% get all spike times during ripple events
disp(['Getting spikes during ripple events for F' num2str(dayindex(1)) num2str(dayindex(2))]);
for recIdx = 1:size(recInfo,2)
    celldata = [];
    spikecounts = [];
    
    if ~isempty(ripples{dayindex(1)}{dayindex(2)}{recs(recIdx)}.startind)
        nonThetaPeriods = [nonthetas{dayindex(1)}{dayindex(2)}{recs(recIdx)}.startind nonthetas{dayindex(1)}{dayindex(2)}{recs(recIdx)}.endind];
        inclrips = isExcluded(ripples{dayindex(1)}{dayindex(2)}{recs(recIdx)}.midind, nonThetaPeriods); %ripples during non-thetas
        
        %get ripple start and end times
        %rippleStarts = ripples{dayindex(1)}{dayindex(2)}{recs(recIdx)}.starttime(logical(inclrips)) + 1; %times are in seconds
        %rippleEnds = ripples{dayindex(1)}{dayindex(2)}{recs(recIdx)}.endtime(logical(inclrips)) + 1;
        rippleMids = ripples{dayindex(1)}{dayindex(2)}{recs(recIdx)}.midtime(logical(inclrips)) + 1;
        rippletimes = [rippleMids-0.5 rippleMids+0.5];

        %loop through the different clusters
        for clusterIdx = 1:length(clusterInfo.labels)
            %get spike times to include
            spikeTimes = recInfo(recIdx).spikeTimes(clusterIdx).TS;
            ripSpikesIdx = isExcluded(spikeTimes,rippletimes);
            spikeTimesDurRips = spikeTimes(logical(ripSpikesIdx));
            spikebins = periodAssign(spikeTimesDurRips,rippletimes);
            
            %get spikes that occur within decoding blocks
            if ~isempty(spikeTimesDurRips)
                validspikes = find(spikebins);
                spikeTimesDurRips = spikeTimesDurRips(validspikes);
                spikebins = spikebins(validspikes);
            end
            
            %load data into matrix
            tmpcelldata = [spikeTimesDurRips' spikebins];
            tmpcelldata(:,3) = clusterIdx;
            celldata = [celldata; tmpcelldata]; %append cluster data
            
            %get spikecounts
            spikecount = zeros(1,size(rippletimes,1));
            for i = 1:length(spikebins)
                spikecount(spikebins(i)) = spikecount(spikebins(i))+1;
            end
            spikecounts = [spikecounts; spikecount]; %number of spikes per block per cell, each row is a cell and each column is a block
        end
        
        celldata = sortrows(celldata,1); 
        cellcounts = sum((spikecounts > 0));
        eventindex = find(cellcounts >= cellcountthresh); %blocks with at least X units firing
    
        for blockIdx = 1:length(eventindex)
            tmpind = find(celldata(:,2) == eventindex(blockIdx));
            out{recIdx}.riptime(blockIdx,1:2) = rippletimes(eventindex(blockIdx),[1 2]); %block time
            out{recIdx}.ripdata(blockIdx).spiketimes = celldata(tmpind,1); %spike times
            out{recIdx}.ripdata(blockIdx).clusindex = celldata(tmpind,3); %cell info
        end
        
        if isempty(eventindex)
            out{recIdx} = nan;
        end
    else
        out{recIdx} = nan;
    end
end
end