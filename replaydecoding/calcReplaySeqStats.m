function statsOutput = calcReplayStats(decodingOutput)
% SP 2.5.19 this function calcs replay stats using the
% decodingSequencesDurRipples output

%% get slope and r-values for decoded data
totalsamples = 10000;
nonzerobins = find(decodingOutput.totalspikes > 0);
rvalues = [];
slopes = [];
for rloop = 1:10
    tBinPicks = distsample(totalsamples,decodingOutput.totalspikes); % sample from the total samples distribution weighted by total spike counts
    regressdata = [];
    for i = 1:length(nonzerobins)
        if (decodingOutput.totalspikes(nonzerobins(i)) > 0)
            tmpnumsamples = sum(tBinPicks == nonzerobins(i)); %proportional sampling from nonzero bins based on number of spikes in bin? 
            distpicks = distsample(tmpnumsamples,decodingOutput.spatialprob(:,nonzerobins(i))); %samples from PDF
            distpicks(:,2) = i; %would this be better as distpicks(:,2) = timebins(nonzerobins(i))
            regressdata = [regressdata; distpicks]; %assigns sampled locations to the bin number
        end
    end
    regressdata(:,3) = 1;
    [b,bint,r,rint,stats] = regress(regressdata(:,1),regressdata(:,2:3));
    
    % get the r^2 value
    rvalues = [rvalues; stats(1)];
    slopes = [slopes; b(1)];
end

%% get slope and r-values for shuffled data
% note - this was done a different way in Annabelle's old code that I'm going to not do right now
scrambleddata = []; scrambledslopes = [];
nonzerobinsshuffled = nonzerobins;
for iteration = 1:10000
    %get PDFs with shuffled bins
    nonzerobinsshuffled = nonzerobinsshuffled(randperm(length(nonzerobinsshuffled)));
   % decodingOutputShuffled.spatialprob = zeros(size(decodingOutput.spatialprob,1),length(nonzerobinsshuffled));
    for binIdx = 1:length(nonzerobins) %loop through new bins to assign spatial prob from random permutation of original decoding
        decodingOutputShuffled.totalspikes(binIdx) = decodingOutput.totalspikes(nonzerobinsshuffled(binIdx));
        decodingOutputShuffled.spatialprob(:,binIdx) = decodingOutput.spatialprob(:,nonzerobinsshuffled(binIdx));
    end
    
    %get slope and r-values for shuffled data
    tBinPicks = distsample(totalsamples,decodingOutput.totalspikes); % sample from the total samples distribution weighted by total spike counts
    regressdata = [];
    for i = 1:length(nonzerobins)
        if (decodingOutput.totalspikes(nonzerobins(i)) > 0)
            tmpnumsamples = sum(tBinPicks == nonzerobins(i)); %proportional sampling from nonzero bins based on number of spikes in bin? 
            distpicks = distsample(tmpnumsamples,decodingOutputShuffled.spatialprob(:,i)); %samples from PDF
            distpicks(:,2) = i; %would this be better as distpicks(:,2) = timebins(nonzerobins(i))
            regressdata = [regressdata; distpicks]; %assigns sampled locations to the bin number
        end
    end
    regressdata(:,3) = 1;
    [b,bint,r,rint,tmpstats] = regress(regressdata(:,1),regressdata(:,2:3));
    
    % get the r^2 value
    scrambleddata = [scrambleddata; tmpstats(1)];
    scrambledslopes = [scrambledslopes; b(1)];
end

%% compare normal and shuffled data to get significance
%method 1 - get a p-value
outvector = [];
for i = 1:length(rvalues)
    outvector(i) = sum(rvalues(i) < scrambleddata)/length(scrambleddata); %proportion of r squared values that are less than the scrambled data, eg p value
end
statsOutput.slope = mean(slopes);
statsOutput.rval = mean(rvalues);
statsOutput.pval = mean(outvector); %slope R^2 p

%method 2 - check if r^2 is in greater than 95th percentile
outvector = [];
for i = 1:length(rvalues)
    % significance check
    rvaluethreshold = prctile(scrambleddata,95);
    signif(i) = logical(mean(rvalues) > rvaluethreshold);
    
    %get percentile
    nless = sum(scrambleddata < mean(rvalues));
    nequal = sum(scrambleddata == mean(rvalues));
    rvalprctile(i) = 100*(nless + 0.5*nequal)/length(scrambleddata);
end
statsOutput.rvalprctile = mean(rvalprctile);
statsOutput.signif = mean(signif);

end
