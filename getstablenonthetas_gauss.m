function thetaInd = getstablenonthetas_gauss(index, thetadir, FRInfo, clusterInfo, recInfo)
% this function finds stable theta periods and converts them into start and
% end indices that can be used for calculating mean/peak firing rate 
% NJ 05/05/18

for clst = 1:length(FRInfo.stabletimes) % for all clusters
    stabletemp = FRInfo.stabletimes{1,clst};
    
    for rec = 1:size(index,1) % for all recordings
        filename = [thetadir 'nonthetas' num2str(index(rec,3)) '.mat'];
        load(filename); % loads each theta file
                
        thetaInd{clst}{rec} = [];
                   
        % converting theta times, currently broken into multiple recording
        % files with restarting starttime at zero, to actual times with 
        % respect to entire recording day (stitching recorded files back to
        % back)

        if rec == 1            
            previousRecEnd = 0;% first file of the day does not need change
        else
            previousRecEnd = previousRecEnd + recInfo(rec-1).duration;

        end
        newThetaStart = nonthetas{index(rec,1)}{index(rec,2)}{end}.startind/20000 + previousRecEnd; %may need to subtract a second as well
        newThetaEnd =  nonthetas{index(rec,1)}{index(rec,2)}{end}.endind/20000 + previousRecEnd;
        newStableStart = stabletemp(rec,1) + previousRecEnd;
        newStableEnd = stabletemp(rec,2) + previousRecEnd;

        % start and end times of theta and stable periods in a n-by-2 matrix
        thetaTimes = [newThetaStart, newThetaEnd];
        stableTimes = [newStableStart, newStableEnd]; % picks up start and end stable times from each recording
        
        % 'both' uses the isExcluded function to assign values that are
        % marks periods that are both stable and during theta as 1
        both = [isExcluded(thetaTimes(:,1), stableTimes), isExcluded(thetaTimes(:,2),stableTimes)];
        
        if size(thetaTimes,1) ~= 0 % accounts for when there is no stable theta period within a cluster
            for jjj = 1:size(thetaTimes,1) % for all theta periods
                if both(jjj,1) == 1 && both(jjj,2) == 1 % ensures both start and end times of a theta period is within stable times
                    thetatemp = (thetaTimes(jjj,:)-1)*2000/20; % converts theta times to indices by multipling by theta sampling rate and dividing by conversion factor 
                    thetaInd{clst}{rec} = [thetaInd{clst}{rec}; thetatemp]; % stores stable theta indices per recording for future calculations
                else
                    continue
                end
            end
        else % if there is no stable theta period
        end
    end
end

