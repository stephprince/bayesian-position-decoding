function plotEstActPosition_100sblocks(positionEstimate, positionInfo, eventTimes, dayindex, savedfiguresdir, recIdx,PFbinsize)
%SP 9.27.18
% this function plots the blocks of decoded data above the animals actual
% position in 100s blocks in order to visualize how well the decoding is
% working on a trial by trial basis

%% set parameters
plotbinsize = 0.02
totaltime = 0:plotbinsize:100-plotbinsize;
posEstOverTime = nan(360/PFbinsize,length(totaltime));
counter = 0;
if isempty(positionEstimate)
    return
elseif isempty(positionEstimate{recIdx})
    return
end

%% loop through the events to compile 100s blocks of data
for eventIdx = 1:size(positionEstimate{recIdx},2)    
    if ~isempty(eventTimes{eventIdx})
        % if end of a recording
        endofrec = 0;
        if totaltime(end) > positionInfo(recIdx).timeSmooth(end)
            totaltime(end) = positionInfo(recIdx).timeSmooth(end);
            endofrec = 1;
        end
        
        % if still compiling data into 100s block
        if (totaltime(end) > eventTimes{eventIdx}(end)) && endofrec ~=1
            if ~isempty(positionEstimate{recIdx}{eventIdx})
                blockSize = eventTimes{eventIdx}(1):plotbinsize:eventTimes{eventIdx}(2)-plotbinsize;
                if ~isempty(blockSize)
                    posEstForBlock = repmat(positionEstimate{recIdx}{eventIdx},1,size(blockSize,2));
                    timeIdx = find(ismember(round(totaltime,2),round(eventTimes{eventIdx},2)));
                    if isempty(timeIdx)
                        timeIdx = lookup2(round(eventTimes{eventIdx},2),round(totaltime,2));
                        if size(posEstForBlock,2) < length(timeIdx(1):timeIdx(2)-1)
                            %this is a not great correction but idk if there's a better way
                            timeIdx(2) = timeIdx(2)-1;
                        elseif size(posEstForBlock,2) > length(timeIdx(1):timeIdx(2)-1)
                            timeIdx(2) = timeIdx(2)+1;
                        end
                    end
                    posEstOverTime(:,timeIdx(1):timeIdx(2)-1) = posEstForBlock;
                end
            end
            
            % have 100s of data so ready to plot
        elseif ((totaltime(end) < eventTimes{eventIdx}(end)) && (sum(sum(posEstOverTime)) ~= 0)) && endofrec ~=1
            counter = counter+1;
            
            figure; hold on;
            %plot the estimated position
            subplot(2,1,1)
            imagesc(totaltime,[0:PFbinsize:358],posEstOverTime);
            set(gca,'color',[1 1 1])
            set(gca,'YDir','normal')
            oldclims = caxis;
            caxis([(oldclims(1)+(oldclims(2)-oldclims(1))/10) oldclims(2)]);
            title('Decoded Position')
            ylabel('Degrees')
            
            %plot the actual position
            %%% NOTE, WEIRD STUFF IF YOUR SAMPLING RATE IS
            %%% OFF OR IF THE SMOOTH TIMES DON'T START AT ZERO - BE CAREFUL
            subplot(2,1,2)
            if totaltime(1) < positionInfo(recIdx).timeSmooth(1)
                startpos = 1;
                timebeforeposdata = find(round(totaltime,2) == round(positionInfo(recIdx).timeSmooth(1),2));
                if isempty(timebeforeposdata)
                    timebeforeposdata = lookup2(round(positionInfo(recIdx).timeSmooth(1),2),round(totaltime,2));
                end
            else
                startpos = find(round(positionInfo(recIdx).timeSmooth,2) == round(totaltime(1),2));
            end
            endpos = find(round(positionInfo(recIdx).timeSmooth,2) == round(totaltime(end),2));
            
            %catch if the sampling rates are slightly offset so there are no exact
            %matches
            if isempty(startpos)
                startpos = lookup2(round(totaltime(1),2),round(positionInfo(recIdx).timeSmooth,2));
            end
            if isempty(endpos)
                endpos = lookup2(round(totaltime(end),2), round(positionInfo(recIdx).timeSmooth,2));
            end
            
            posActualOverTime = positionInfo(recIdx).thetaSmooth(startpos:endpos);
            if totaltime(1) < positionInfo(recIdx).timeSmooth(1)
                posActualOverTime = [nan(timebeforeposdata-1,1); posActualOverTime];
            end
            if length(totaltime) ~= length(posActualOverTime)
                posActualOverTime_resampled = resample(posActualOverTime, length(totaltime), length(posActualOverTime));
            else
                posActualOverTime_resampled = posActualOverTime;
            end
            plot(totaltime,posActualOverTime_resampled,'k','LineWidth',2)
            title('Actual Position')
            ylabel('Degrees')
            xlabel('Time (s)')
            set(gcf,'units','normalize','outerpos',[0 0 1 1])
            suptitle(['Actual vs. Decoded Position for 100s periods - F' num2str(dayindex(1)) '_' num2str(dayindex(2)) ' - ' num2str(counter)])
            
            %save thd plot
            animaldir = [savedfiguresdir 'F' num2str(dayindex(1)) '_' num2str(dayindex(2)) '\'];
            filename = [animaldir 'decodedpos_100sblocks_rec' num2str(recIdx) 'num' num2str(counter)];
            saveas(gcf,filename,'png'); saveas(gcf,filename,'fig');
            
            %reset variables
            totaltime = totaltime + 100;
            posEstOverTime = nan(360/PFbinsize,length(totaltime));
            
            %if data runs into the end of the recording
        elseif endofrec == 1
            totaltime = totaltime+100;
            posEstOverTime = nan(360/PFbinsize,length(totaltime));
        end
    end
end

end