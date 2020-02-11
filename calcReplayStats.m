function out = calcReplayStats(activespiketimes,expSpikeCounts,timebins,distvector)

totalsamples = 10000;
spikedata = [];
cellsactive = [];
exponentmatrix = exp(-expSpikeCounts);
%histogram the spikes for each active cell
for i = 1:length(activespiketimes)
    spikebins = lookup(activespiketimes{i},timebins);
    spikecount = zeros(1,length(timebins));
    for j = 1:length(spikebins)
        spikecount(spikebins(j)) = spikecount(spikebins(j))+1;
    end
    spikedata = [spikedata; spikecount];
    cellsactive = [cellsactive; (spikecount > 0)];
end

%the decoded data contains spatial probabilities, and is x by t
decodedata = zeros(size(expSpikeCounts,2),size(spikedata,2));
naninds = find(isnan(expSpikeCounts(1,:)));
for t = 1:size(spikedata,2) %time  
    Tspikecount = repmat(spikedata(:,t),1,size(expSpikeCounts,2));
    %calculate P(spikecount|x) for this timebin across all cells and all x 
    spatialprob = prod(((expSpikeCounts.^Tspikecount)./factorial(Tspikecount)).*exp(-expSpikeCounts),1)'; 
    spatialprob(naninds) = 0;
    spatialprob = spatialprob/sum(spatialprob);  %normalize across space
    %spatialprob(find(spatialprob < max(spatialprob/2))) = 0;
    decodedata(:,t) = spatialprob;  
end

totalspikecounts = sum(spikedata,1);
totalcellsactive = sum(cellsactive,1);
%totalcellsactive = smoothvect(totalcellsactive,gaussian(1,8));
%totalcellsactive = smoothvect(totalspikecounts,gaussian(1,8));
%disp(totalcellsactive)
nonzerobins = find(totalspikecounts > 0);
%nonzerobins = find(totalcellsactive > (length(activespiketimes)/20));
%nonzerobins = find(totalcellsactive > 1);


rvalues = [];
slopes = [];
for rloop = 1:10
    % sample from the total samples distribution weighted by total spike counts
    tBinPicks = distsample(totalsamples,totalspikecounts);
    regressdata = [];
    for i = 1:length(nonzerobins)
        if (totalspikecounts(nonzerobins(i)) > 0)
            tmpnumsamples = sum(tBinPicks == nonzerobins(i));
            distpicks = distvector(distsample(tmpnumsamples,decodedata(:,nonzerobins(i))))';
            %distpicks(:,2) = timebins(nonzerobins(i));
            distpicks(:,2) = i;
            regressdata = [regressdata; distpicks];
        end
    end
    regressdata(:,3) = 1;
    [b,bint,r,rint,stats] = regress(regressdata(:,1),regressdata(:,2:3));
    % get the r^2 value
    rvalues = [rvalues; stats(1)];
    slopes = [slopes; b(1)];
end

%figure
%plot(regressdata(:,2),regressdata(:,1),'.');

scrambleddata = [];
nonzerobins2 = nonzerobins;
for iteration = 1:500
    %the decoded data contains spatial probabilities, and is x by t
    
    %permindex = randperm(size(expSpikeCounts,1)); 
    permindex = 1:size(expSpikeCounts,1);
    
    nonzerobins2 = nonzerobins2(randperm(length(nonzerobins2)));
    
    tmpexpSpikeCounts = expSpikeCounts(permindex,:);
    tmpexponentmatrix = exponentmatrix(permindex,:);
    decodedata2 = zeros(size(expSpikeCounts,2),length(nonzerobins2));
    
    
    for t = 1:length(nonzerobins) %time
        Tspikecount = repmat(spikedata(:,nonzerobins2(t)),1,size(expSpikeCounts,2));
        factorialmatrix = repmat(factorial(spikedata(:,nonzerobins2(t))),1,size(expSpikeCounts,2));
        
        %calculate P(spikecount|x) for this timebin across all cells and all x
        spatialprob = prod(((tmpexpSpikeCounts.^Tspikecount)./factorialmatrix).*tmpexponentmatrix,1)';
        spatialprob(naninds) = 0;
        spatialprob = spatialprob/sum(spatialprob);  %normalize across space
        %spatialprob(find(spatialprob < max(spatialprob/2))) = 0;
        decodedata2(:,t) = spatialprob;
    end
    
    
    regressdata = [];
    tBinPicks = distsample(totalsamples,totalspikecounts);
    for i = 1:length(nonzerobins2)
        if (totalspikecounts(nonzerobins(i)) > 0)
            tmpnumsamples = sum(tBinPicks == nonzerobins(i));
            distpicks = distvector(distsample(tmpnumsamples,decodedata2(:,i)))';
            %distpicks(:,2) = timebins(nonzerobins(i));
            distpicks(:,2) = i;
            regressdata = [regressdata; distpicks];
        end
    end
    regressdata(:,3) = 1;
    [b,bint,r,rint,tmpstats] = regress(regressdata(:,1),regressdata(:,2:3));
    scrambleddata = [scrambleddata; tmpstats(1)];
    
    if (iteration == 100) %early break for cells that are clearly not significant
        tmpresult = [];
        for i = 1:length(rvalues)
            tmpresult(i) = sum(rvalues(i) < scrambleddata)/length(scrambleddata);
        end
        if (mean(tmpresult) > .2)
            break
        else %do more of the initial regressions
            for rloop = 1:490
                tBinPicks = distsample(totalsamples,totalspikecounts);
                regressdata = [];
                for i = 1:length(nonzerobins)
                    if (totalspikecounts(nonzerobins(i)) > 0)
                        tmpnumsamples = sum(tBinPicks == nonzerobins(i));
                        distpicks = distvector(distsample(tmpnumsamples,decodedata(:,nonzerobins(i))))';
                        %distpicks(:,2) = timebins(nonzerobins(i));
                        distpicks(:,2) = i;
                        regressdata = [regressdata; distpicks];
                    end
                end
                regressdata(:,3) = 1;
                [b,bint,r,rint,stats] = regress(regressdata(:,1),regressdata(:,2:3));
                rvalues = [rvalues; stats(1)];
                slopes = [slopes; b(1)];
            end
        end
    end
end

outvector = [];
for i = 1:length(rvalues)
    outvector(i) = sum(rvalues(i) < scrambleddata)/length(scrambleddata); %proportion of r squared values that are less than the scrambled data, eg p value
end
out = [mean(slopes) mean(rvalues) mean(outvector)]; %slope R^2 p

