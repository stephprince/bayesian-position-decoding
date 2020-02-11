function out = getpopulationevents2(indices, excludeperiods, spikes, linpos, pos, ripples, cellinfo, task, cellcountthresh, ripVeqn, ripwelldist)
%
%Calculates the binned spike counts of all neurons in the index list.
%
%spikes - the 'spikes' cell array for the day you are analyzing
%linpos - the output of LINEARDAYPROCESS for the day you are analyzing.
%indices - [day epoch tetrode cell]
%timebin- the length of each temporal bin (default 0.01 sec)
%excludeperiods - [start end] times for each exlcude period
%ripVeqn, velocity eqn, eg '<=1' in cm
%ripwelldist distance from center well, suggest 20
%
%In the output, the spikecounts field is n by t, where n is the number of
%cells and t is the number of timebins.

index = [indices(1,1) indices(1,2)];

%get allowable ripple times
[y ind] = min(linpos{index(1,1)}{index(1,2)}.statematrix.linearDistanceToWells,[],2); %y is distance to nearest well, ind is nearest well
linv = calclinvelocity(linpos, index(1,:));
wellinfo = [abs(linv.velocity) y ind]; %linear velocity & dist to nearest welln & closest well
times = linpos{index(1,1)}{index(1,2)}.statematrix.time;
riptimeinclude = zeros(size(times));
eval(['[ripindex j] = find(wellinfo(:,1)', ripVeqn, ...
    ' & wellinfo(:,2) <= ripwelldist & wellinfo(:,3) == 1 );']);
riptimeinclude(ripindex) = 1;
ripexcludeperiods = getExcludePeriods(times, riptimeinclude);

posdata = pos{index(1)}{index(2)}.data;
%statematrix = linpos{index(1)}{index(2)}.statematrix;
riptimes = getripples(index, ripples, cellinfo, 'cellfilter', '(isequal($area, ''CA1''))','excludeperiods', ripexcludeperiods,'minstd',3);

[traj correct passnum] = getrippassinfo(index, linpos, task, riptimes);
[ripfutpast] = getripcodefutpast( indices, linpos, spikes, riptimes, traj);
%ripfutpast, -1 to select past, 1 to select future, 0 for both


out.index = [];
out.eventtraj = [];
out.eventcorrect = []; %ripple correct or incorrect
out.eventfutpast = []; % ripple agree with future (1) or past (-1)
out.eventpassnum = [];
out.eventdist = [];
out.eventtime = [];
out.preeventcount = [];
out.eventimmobiletime = [];
out.eventdata = [];
out.peak = 0;

spikecounts = [];
celldata = [];

if ~isempty(riptimes)
    %go through each cell and calculate the binned spike counts
    for cellcount = 1:size(indices,1) %for each cell
        
        index = indices(cellcount,:);
        if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)}.data)
            spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);
        else
            spiketimes = [];
        end
        spiketimes = spiketimes(find(~isExcluded(spiketimes, excludeperiods))); %get included spike times for this cell
        spikebins = periodAssign(spiketimes, riptimes(:,[1 2])); %assings each included spike to a ripple bin, index of ripple each spike belongs to
        
        %spikebins = lookup(spiketimes,riptimes(:,2));
        if ~isempty(spiketimes)
            validspikes = find(spikebins);
            spiketimes = spiketimes(validspikes); %spikes that occur during a ripple
            spikebins = spikebins(validspikes);
            %             spikedeviation = abs(spiketimes - riptimes(spikebins,2));
            %             validspikes = find(spikedeviation <= (window/2));
            %             spiketimes = spiketimes(validspikes);
            %             spikebins = spikebins(validspikes);
        end
        tmpcelldata = [spiketimes spikebins]; %spike time and ripple index
        tmpcelldata(:,3) = cellcount; %cell number
        celldata = [celldata; tmpcelldata]; %append to info from other cells
        spikecount = zeros(1,size(riptimes,1));
        for i = 1:length(spikebins)
            spikecount(spikebins(i)) = spikecount(spikebins(i))+1;
        end
        
        spikecounts = [spikecounts; spikecount]; %number of spikes per ripple for each cell.  each row is cell, each column is aripple/
        out.index = [out.index; index];
    end
    celldata = sortrows(celldata,1); %sort all spikes by time
    
    
    
    %newtimebinsInd = find(~isExcluded(timebins, excludeperiods));
    %newtimes = timebins(newtimebinsInd);
    
    cellcounts = sum((spikecounts > 0)); %each ripple with spikes from at least 1 cell
    %cellcounts = cellcounts(newtimebinsInd);
    eventindex = find(cellcounts >= cellcountthresh); %ripples with at least cellcountthresh cells active
    
    timeindex = lookup(riptimes(:,2),posdata(:,1)); %look up pos times that correspond to ripple events
    
    tmpvel = posdata(:,5);
    %tmpvel = tmpvel<(max(tmpvel)*.05);
    tmpvel = tmpvel<(2);
    
    
    resetpoints = find(diff(tmpvel) < 0)+1;
    immobiletime = cumsum(tmpvel);
    for i = 1:length(resetpoints)
        immobiletime(resetpoints(i):end) = immobiletime(resetpoints(i):end)-immobiletime(resetpoints(i)-1);
    end
    immobiletime = immobiletime/30;
    immobile = immobiletime(timeindex);
    
    
    %timeindex = lookup(riptimes(:,2),statematrix.time);
    %traj = statematrix.traj(timeindex);
    %dist = statematrix.lindist(timeindex);
    
    for event = 1:length(eventindex) % for each included ripple event
        out.eventtime(event,1:2) = riptimes(eventindex(event),[1 2]); %start and end time
        out.eventimmobiletime(event,1) = immobile(eventindex(event)); %time been immobile
        out.eventtraj(event) = traj(eventindex(event)); %trajectory from statematrix
        out.eventcorrect(event) = correct(eventindex(event)); %ripple correct or incorrect
        out.eventpassnum(event) = passnum(eventindex(event)); %ripple correct or incorrect
        out.eventfutpast(event) = ripfutpast(eventindex(event)); % ripple agree with future (1) or past (-1)
        tmpind = find(celldata(:,2) == eventindex(event));
        out.eventdata(event).spiketimes = celldata(tmpind,1); %spike times
        out.eventdata(event).cellindex = celldata(tmpind,3); %cell info
        
    end
else
    warning('No ripples found')
end



function out = periodAssign(times, periods)
% out = periodAssign(times, periods)
% TIMES is a vector of times
% PERIODS is an N by 2 list of start and end times
% Returns the index of the period that each time falls into.  If a time
% does not fall into one of the periods, a zero is returned.
% This function assumes that the periods are not overlapping.
%

if ~isempty(periods)
    oneborder = [(periods(:,1)-.0000001);periods(:,2)+.0000001];
    oneborder(:,2) = 0;
    insideborder = [(periods(:,1)+.0000001) (1:size(periods,1))'; (periods(:,2)-.0000001) (1:size(periods,1))'];
    sortedMatrix = [[-inf 0]; sortrows([oneborder;insideborder],1); [inf 0]];
else
    sortedMatrix = [[-inf 0]; [inf 0]];
end
out = sortedMatrix(lookup(times,sortedMatrix(:,1)),2);
