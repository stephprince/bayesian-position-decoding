function out = decodePositionDurTheta(testingdata,trainingdata,timebins)
% SP 9.18.18
% this function usins Bayesian population decoding to decode expected
% position from place fields and compare to the actual position

%% initialize variables
spikedata = [];
cellsactive = [];

%% get histogram of spikes for each active cell
for spikeIdx = 1:length(testingdata)
    spikebins = lookup2(testingdata{spikeIdx},timebins);
    spikecount = zeros(1,length(timebins));
    for j = 1:length(spikebins)
        spikecount(spikebins(j)) = spikecount(spikebins(j))+1;
    end
    spikedata = [spikedata; spikecount];
    cellsactive = [cellsactive; (spikecount > 0)]; 
end

%% decode spatial probabilities from test spike counts
decodedata = zeros(size(trainingdata,2),size(spikedata,2));
naninds = find(isnan(trainingdata(1,:)));
for timeIdx = 1:size(spikedata,2) %loop through time
    Tspikecount = repmat(spikedata(:,timeIdx),1,size(trainingdata,2));
    
    %get p(spikecount|x) for this timebin across all cells and all x
    spatialprob = prod(((trainingdata.^Tspikecount)./factorial(Tspikecount)).*exp(-trainingdata),1)'; 
    spatialprob(naninds) = 0;
    
    %normalize across space to make the probabilities add up to 1
    spatialprob = spatialprob/sum(spatialprob);
    decodedata(:,timeIdx) = spatialprob;
end

totalspikecounts = sum(spikedata,1);
totalcellsactive = sum(cellsactive,1);
nonzerobins = find(totalspikecounts > 0);

out.spatialprob = decodedata;
out.totalspikes = totalspikecounts;
out.totalcellsactive = totalcellsactive;

end