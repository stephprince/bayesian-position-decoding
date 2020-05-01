function [trainingtimes testingtimes] = splittraintestdata(spikedatadir, recs, splittype)
% SP 9.27.18
% this function splits up time series data into blocks for training and
% testing amodel

%%%%% NOTE: THIS IS NOT RANDOMIZED, DATA IS BROKEN UP INTO BLOCKS %%%%%

% OPTIONS:
% splittype - '8020', '5050', '5050largebins'

%% split up data set into blocks
if ~isempty(recs)
    [clusterInfo recInfo] = getClusterInfo(spikedatadir, recs);
    for recIdx = 1:numel(recs)
        %get ratio of data split based on inputs
        if strcmp(splittype, '8020') % this gives you a 4 s training and 1 s testing period
            tempstarttimes = 0:5:round(recInfo(recIdx).duration);
            tempendtimes = tempstarttimes-1;
            starttimes = tempstarttimes(1:end-1);
            endtimes = tempendtimes(2:end);
        elseif strcmp(splittype, '5050') % this gives you a 1 s training and 1 s testing period
            starttimes = 0:2:round(recInfo(recIdx).duration);
            endtimes = 1:2:round(recInfo(recIdx).duration);
        elseif strcmp(splittype,'5050large') % this gives you a 4 s training and 4 s testing period for bigger bins
            starttimes = 0:8:round(recInfo(recIdx).duration);
            endtimes = 4:8:round(recInfo(recIdx).duration);
        elseif strcmp(splittype,'5050verylarge') % this gives you a 10 s training and 10 s testing period for bigger bins
            starttimes = 0:20:round(recInfo(recIdx).duration);
            endtimes = 10:20:round(recInfo(recIdx).duration);
        end
        
        %split data
        if length(starttimes) > length(endtimes); starttimes = starttimes(1:end-1); end;
        trainingtimes{recIdx} = [starttimes', endtimes'];
        testingtimes{recIdx} = [endtimes(1:end-1)',starttimes(2:end)' ];
    end
else
    trainingtimes{1} = [];
    testingtimes{1} = [];
end

end