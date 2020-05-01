function out = calcepochreplaystats(trainingfilter, decodefilter)

%allow updates to propagate to main workspace
global decodefilter

trajmapping = [1 1 2 2];
binsize = .015; %default temporal bin
minspikes = 2;
out = [];
for animalnum = 1:length(decodefilter)
    for g = 1:length(decodefilter(animalnum).output)
        trainingindex = g;
        for epochnum = 1:length(decodefilter(animalnum).output{g})
            
            %get values to determine what rates past CP on each traj
            ratestrj1 = (trainingfilter(animalnum).output{trainingindex}(epochnum).traj == 1 | trainingfilter(animalnum).output{trainingindex}(epochnum).traj == 2) & trainingfilter(animalnum).output{trainingindex}(epochnum).dist > 80;
            ratestrj3 = (trainingfilter(animalnum).output{trainingindex}(epochnum).traj == 3 | trainingfilter(animalnum).output{trainingindex}(epochnum).traj == 4) & trainingfilter(animalnum).output{trainingindex}(epochnum).dist > 80;
            
            for eventindex = 1:length(decodefilter(animalnum).output{g}(epochnum).eventdata)
                disp(eventindex)
                trainingdata = [];
                spikedata = [];
                decodedata = [];
                indexlist = [];
                activespiketimes = [];
                activerates = [];
                
                %pick out all the matching cells from the training data and the
                %decoding data
                %traindata contains linear rates, and is n by x, where n is the
                %number of cells and x is the number of spatial bins
                %spikedata contains spikecounts, and is n by t, where t is the
                %number of temporal bins in the data to be decoded.
                matches = rowfind(trainingfilter(animalnum).output{trainingindex}(epochnum).index(:,[1 3 4]),decodefilter(animalnum).output{g}(epochnum).index(:,[1 3 4])); %find the matching cell indices
                startevent = decodefilter(animalnum).output{g}(epochnum).eventtime(eventindex,1);
                endevent = decodefilter(animalnum).output{g}(epochnum).eventtime(eventindex,2);
                if ((endevent-startevent) < 2)
                    timebins = startevent:binsize:endevent;
                    eventcellsactive = [];
                    activecount = 0;
                    for trainingcell = 1:length(matches)%for each cell in the training data
                        if (matches(trainingcell) > 0) %we have a match
                            indexlist = [indexlist; trainingfilter(animalnum).output{trainingindex}(epochnum).index(trainingcell,:)];
                            trainingdata = [trainingdata; trainingfilter(animalnum).output{trainingindex}(epochnum).rates(trainingcell,:)];
                            tmpspiketimes = decodefilter(animalnum).output{g}(epochnum).eventdata(eventindex).spiketimes(find(decodefilter(animalnum).output{g}(epochnum).eventdata(eventindex).cellindex == matches(trainingcell)));
                            %save all the info for the active cells, eg
                           % cells active in this ripple
                            if ~isempty(tmpspiketimes)
                                activecount = activecount+1;
                                activespiketimes{activecount} = tmpspiketimes;
                                activerates = [activerates; trainingfilter(animalnum).output{trainingindex}(epochnum).rates(trainingcell,:)];
                            end
                        end
                    end
                    trainingdata = trainingdata*binsize; %transform rates to expected number of spikes
                    activerates = activerates*binsize;
                    
                    %determine if codes for future or past. pos = more
                    %future, neg = more past
                    if ~isempty(activerates)
                        if decodefilter(animalnum).output{g}(epochnum).eventtraj(eventindex) == 1
                            ratiorates = (nansum(nansum(activerates(:,ratestrj1))) - nansum(nansum(activerates(:,ratestrj3))))/(nansum(nansum(activerates(:,ratestrj1)))+nansum(nansum(activerates(:,ratestrj3))));
                        elseif decodefilter(animalnum).output{g}(epochnum).eventtraj(eventindex) == 3
                            ratiorates = (nansum(nansum(activerates(:,ratestrj3))) - nansum(nansum(activerates(:,ratestrj1))))/(nansum(nansum(activerates(:,ratestrj1)))+nansum(nansum(activerates(:,ratestrj3))));
                        else
                            ratiorates = NaN;
                        end
                        activpastCP = any(any(activerates(:,(ratestrj1 | ratestrj3))./binsize>3));
                        
                        %determine most distant point with >25% of peak activerate
                        maxdist1 = max( trainingfilter(animalnum).output{trainingindex}(epochnum).dist(:,(sum(activerates,1)>max(sum(activerates,1))*.25)) );
                        
                        %find midpoint of pdf for traj with peak
                        tempactiverates = sum(activerates,1);
                        tempactiverates(isnan(tempactiverates)) = 0;
                        [m ind] = max(tempactiverates);
                        trj = trainingfilter(animalnum).output{trainingindex}(epochnum).traj(ind);
                        a = cumsum(tempactiverates(trainingfilter(animalnum).output{trainingindex}(epochnum).traj==trj));
                        ind2 = find(a>=a(end)/2, 1, 'first');
                        dists = trainingfilter(animalnum).output{trainingindex}(epochnum).dist(trainingfilter(animalnum).output{trainingindex}(epochnum).traj==trj);
                        maxdist2 =  dists(ind2);
                        
                        %find peak for each cell and then most distant max
                        [m ind3] = max(activerates, [], 2);
                        maxdist3 = trainingfilter(animalnum).output{trainingindex}(epochnum).dist(:,max(ind3));
                        
                        %determine peak of pdf
                        [m ind4] = max(sum(activerates,1));
                        maxdist4 = trainingfilter(animalnum).output{trainingindex}(epochnum).dist(:, ind4);
                       
                        maxdist = [maxdist1 maxdist2 maxdist3 maxdist4];
                       
                    else
                        ratiorates = [];
                        activpastCP = [];
                        maxdist = [];
                    end
                    
                    if (length(activespiketimes) >= minspikes)
                        %out = [out; calcReplayStats(activespiketimes,activerates,timebins,trainingfilter(animalnum).output{trainingindex}(epochnum).dist)];
                        
                        % rest box stats
                        %rstats = calcReplayStats(activespiketimes,activerates,timebins,trainingfilter(animalnum).output{trainingindex}(epochnum).dist);
                        %out = [out; [rstats decodefilter(animalnum).output{g}(epochnum).eventimmobiletime(eventindex) length(activespiketimes) ]];  %decodefilter(animalnum).output{g}(epochnum).std(eventindex)
                        
                        
                        %out = [out; [decodefilter(animalnum).output{g}(epochnum).eventimmobiletime(eventindex) length(activespiketimes) ]];  %decodefilter(animalnum).output{g}(epochnum).std(eventindex)
                        %out = [out;index(1)];
                        
                        % the following calculates the significance of the replay events
                        % and, with the addition of the global line at the top of the file,
                        % saves the stats for each event in the decode filter
                        rstats = calcReplayStats(activespiketimes,activerates,timebins, ...
                            trainingfilter(animalnum).output{trainingindex}(epochnum).dist);
                        %out = [out; [rstats decodefilter(animalnum).output{g}(epochnum).eventimmobiletime(eventindex) length(activespiketimes) decodefilter(animalnum).output{g}(epochnum).eventloc(eventindex)]];
                        
                        index = [animalnum decodefilter(animalnum).epochs{g}(epochnum,:) g]; %animal day epoch group
                        appenddata = [rstats decodefilter(animalnum).output{g}(epochnum).eventimmobiletime(eventindex)...
                            length(activespiketimes) decodefilter(animalnum).output{g}(epochnum).eventcorrect(eventindex)...
                            ratiorates activpastCP decodefilter(animalnum).output{g}(epochnum).eventpassnum(eventindex) maxdist ];
                        out = [out; index appenddata];
                        %[1:an 2:day 3:epoch 4:group(exposure) 5:slope 6:R^2 7:p
                        %8:immobiletime 9:numcellsactive? 10:in/correct 11:fut/past
                        %12:activpastCP 13:passnum 14-16:maxdist 17:pdfpkdist]
                        
                    else
                        rstats = NaN;
                    end
                    decodefilter(animalnum).output{g}(epochnum).eventdata(eventindex).replaystat = rstats; %slope  R^2  p
                    decodefilter(animalnum).output{g}(epochnum).eventdata(eventindex).eventpopfutpast = [ratiorates activpastCP];
                    decodefilter(animalnum).output{g}(epochnum).eventdata(eventindex).maxdists = maxdist;
                end
            end
        end
    end
end