function [state, lindist] = getbehavestate(linpos, day, epoch, includeStates, varargin)
%[state, lindist] = GETBEHAVESTATE(linpos, day, epoch, includeStates, options)
%
%Finds the linear trajectory for each time point, and finds invalid behavioral times
%
%LINPOS is the output of linearizeposition for one day
%DAY is the day number for linpos
%EPOCH - which run epoch do you want to analyze?
%INCLUDESTATES is a number between 1 and 6.  A higher 
%number tells this program to include progressively more ambiguous behavior
%times into statevector.  (1 is least ambiguous). Here are the definitions
%for each number:
%
%1:  Include all times when the animal is in transit between the two wells
%    of a defined trajectory (either in the foreward or backward direction),  
%    and velocity is higher than minvelocity. This is defined as the times 
%    when the animal is in transit between the two different endpoints and is
%    on a track segment that is part of the trajectory. Also, head
%    direction and motiondir have to be in the same direction.
%2:  Include times not occurring during any of the defined well-to-well trajectories, 
%    and occur when the animal is not leaving and entering the same well,
%    and is in a track segment that only belongs to one of the defined linear
%    trajectories, and velocity > minvelocity. Also, head dir and motiondir
%    have to be in the same direction.
%3:  Fill all times that occur when the animal is on a segment 
%    that belongs to only one trajectory (even if it is leaving and entering the same well), 
%    and velocity > minvelocity. Also, head dir and motiondir
%    have to be in the same direction.
%4:  Include times that occur when the animal was on an amiguous segment (like the
%    home arm), and velocity > minvelocity. Also, head dir and motiondir
%    have to be in the same direction.
%5:  Include any remaining points occuring in the beginning and end of the
%    session with velocity > minvelocity. Also, head dir and motiondir
%    have to be in the same direction.
%6:  Include all points by choosing the closest defined trajectory
%
% Options - 
%           'minlinvelocity' - the minimum linear speed the animal must be
%           moving for the behavior to be counted as valid (default 6cm/sec)
%
%           'min2dvelocity' - the minimum 2d speed. Default [] (not in use)
%                             This requires an additional option entry of
%                             'pos' (the 2d position structure for the day)
%
%           'timerange'     - a 1 by 2 array with the minimum and maximum times
%                             to consider (in seconds). Everything outside this range will be given
%                             a -1 in the 'state' output. Default: [-inf inf]
%
%           'headdir'       - 1 (default) if the animal's head direction and motion
%                             direction must be int he same direction
%                             0 if head direction and motion direction do
%                             not have to agree
%
%State gives a value for which trajectory the animal was on
%for each time in statematrix.  Odd trajectories are when the animal
%is moving in a positive linear direction, and even trajectories are for
%negative directions. If the trajectory is unknown then it gives a -1
%value.  
% NOTE: STATE IS NOT ALWAYS CORRECT.  use plotstatevector to visualize and
% check.  State especially has problems with velocity = 0 or animal tilting
% head so head direction is unclear (e.g. at the well and perhaps at turns, 
% depending on your track).
%
%lindist gives the linear distance traveled in the current trajectory



minvelocity = 6;
min2dvelocity = [];
timerange = [-inf inf];
headdir = 1;

%set variable options
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'minlinvelocity'
            minvelocity = varargin{option+1}; 
        case 'min2dvelocity'
            min2dvelocity = varargin{option+1};
        case 'pos'
            pos = varargin{option+1};
        case 'timerange'
            timerange = varargin{option+1};
        case 'headdir'
            headdir = varargin{option+1};
        otherwise
            error(['Option', varargin{option},'not defined']);
        
    end
end


statematrix = linpos{day}{epoch}.statematrix;
segmenttable = linpos{day}{epoch}.segmenttable;
trajwells = linpos{day}{epoch}.trajwells;

traject = trajwells;
if ~isempty(min2dvelocity)
    vel2d = pos{day}{epoch}.data(:,5);
end
includeStates = 1:includeStates;

%if no reference well is assigned, it is not a valid point
%also, only points inside the designated time range are valid
validPoints = ((statematrix.referenceWell > 0) & (statematrix.time >= timerange(1)) & (statematrix.time <= timerange(2))); 

validPointsIndex = find(validPoints);
trajvector = ones(length(validPointsIndex),1)*-1;
statevector = ones(length(validPointsIndex),1)*-1;
referencewell = statematrix.referenceWell(validPointsIndex);
welltraj = statematrix.wellExitEnter(validPointsIndex,:); %exit and enter wells
segmentindex = statematrix.segmentIndex(validPointsIndex);
samewell = (welltraj(:,1) == welltraj(:,2)); %Which times are during trajectories between the same wells 
diffwell = (welltraj(:,1) ~= welltraj(:,2)); %Which times are during trajectories between different wells
trajcount = rowcount(statematrix.segmentIndex(validPointsIndex),segmenttable(:,3)); %how many of the linear trajectories does each segment belong to

matrixIndex = sub2ind(size(statematrix.linearDistanceToWells),validPointsIndex,referencewell);
lineardistance = statematrix.linearDistanceToWells(matrixIndex);

if headdir == 1  %if head and motion direction should be in same direction
    if ~isempty(min2dvelocity)
        forewarddir = ((statematrix.segmentHeadDirection(matrixIndex) > 0) & (statematrix.linearVelocity(matrixIndex) >= minvelocity) & (vel2d(matrixIndex) > min2dvelocity)); %facing positve dir and moving positive dir
        backwarddir = ((statematrix.segmentHeadDirection(matrixIndex) < 0) & (statematrix.linearVelocity(matrixIndex) <= -minvelocity) & (vel2d(matrixIndex) > min2dvelocity));%facing negative dir and moving negative dir
    else
        forewarddir = ((statematrix.segmentHeadDirection(matrixIndex) > 0) & (statematrix.linearVelocity(matrixIndex) >= minvelocity)); %facing positve dir and moving positive dir
        backwarddir = ((statematrix.segmentHeadDirection(matrixIndex) < 0) & (statematrix.linearVelocity(matrixIndex) <= -minvelocity));%facing negative dir and moving negative dir
    end
elseif headdir == 0 %if head and motion direction do not need to be in same direction 
    if ~isempty(min2dvelocity)
        forewarddir = ( (statematrix.linearVelocity(matrixIndex) >= minvelocity) & (vel2d(matrixIndex) > min2dvelocity)); %facing positve dir and moving positive dir
        backwarddir = ( (statematrix.linearVelocity(matrixIndex) <= -minvelocity) & (vel2d(matrixIndex) > min2dvelocity));%facing negative dir and moving negative dir
    else
        forewarddir = ( (statematrix.linearVelocity(matrixIndex) >= minvelocity)); %facing positve dir and moving positive dir
        backwarddir = ((statematrix.linearVelocity(matrixIndex) <= -minvelocity));%facing negative dir and moving negative dir
    end
end

%LEVEL 1   
for trajnum = 1:size(traject,1)
    currtraj = traject(trajnum,:);
    
    %which times occur when the animal is intransit between the two wells
    %for the current trajectory
    inDefinedTraj(:,trajnum) = (((welltraj(:,1)==currtraj(1))&(welltraj(:,2)==currtraj(2))) | ((welltraj(:,1)==currtraj(2))&(welltraj(:,2)==currtraj(1))));        
    %which segments are in the current trajectory
    segmentsInTraj = segmenttable(find(segmenttable(:,1) == trajnum),3);
    
    %find all times when the animal was on the current trajectory (either
    %in the foreward or backward direction).  This is defined as the times 
    %when the animal is in transit between the two correct endpoints and is
    %on a track segment that is part of the trajectory.
    
    %first find the foreward direction times (both head and motion
    %direction in the positive direction, and velocity greater than
    %minvelocity)       
    findindex = find( inDefinedTraj(:,trajnum)& ismember(segmentindex,segmentsInTraj) & forewarddir);   
    %assign the current trajectory to those indeces
    trajvector(findindex) = 2*(trajnum)-1;
    if ismember(1,includeStates)
        statevector(findindex) = trajvector(findindex);
    end
    %then find the negative direction times
    findindex = find( inDefinedTraj(:,trajnum) & ismember(segmentindex,segmentsInTraj) & backwarddir);        
    trajvector(findindex) = 2*(trajnum);
    if ismember(1,includeStates)
        statevector(findindex) = trajvector(findindex);
    end
end

undefinedindex = (trajvector == -1); %still undefined times
inAnyDefinedTraj = (sum(inDefinedTraj,2) > 0); %These times occur when the animal is in transit between two wells of a defined trajectory

%LEVEL 2
%for the times that are still undefined and not occurring during any of the above well-to-well trajectories, 
%and occur when the animal is not leaving and entering the same well,
%and is in a track segment that only belongs to one of the defined linear
%trajectories, assign the proper trajectory 
findindex = find((undefinedindex) & (diffwell) & (trajcount == 1) & (forewarddir) & (~inAnyDefinedTraj)); 
trajvector(findindex) = 2*(segmenttable(rowfind(segmentindex(findindex),segmenttable(:,3)),1))-1;
if ismember(2,includeStates)
    statevector(findindex) = trajvector(findindex);
end
findindex = find((undefinedindex) & (diffwell) & (trajcount == 1) & (backwarddir) & (~inAnyDefinedTraj)); 
trajvector(findindex) = 2*(segmenttable(rowfind(segmentindex(findindex),segmenttable(:,3)),1));
if ismember(2,includeStates)
    statevector(findindex) = trajvector(findindex);
end

if ((nargout == 1) & (max(includeStates) <= 2))
    return
end

%LEVEL 3
%fill all times when the animal is on a segment that belongs to only one
%trajectory (even if it is leaving and entering the same well)
undefinedindex = (trajvector == -1); %still undefined times
findindex = find((undefinedindex) & (trajcount == 1) & (forewarddir)); 
trajvector(findindex) = 2*(segmenttable(rowfind(segmentindex(findindex),segmenttable(:,3)),1))-1;
if ismember(3,includeStates)
    statevector(findindex) = trajvector(findindex);
end
findindex = find((undefinedindex) & (trajcount == 1) & (backwarddir)); 
trajvector(findindex) = 2*(segmenttable(rowfind(segmentindex(findindex),segmenttable(:,3)),1));
if ismember(3,includeStates)
    statevector(findindex) = trajvector(findindex);
end

%LEVEL 4
%Fill in the undefined times when the animal was on an amiguous segment (like the
%home arm in a w-maze) 
undefinedindex = (trajvector == -1); %the undefined times 
%if facing positive direction, then assign the traj of the next traj in the
%future
findindex = find((undefinedindex) & (trajcount > 1) & (forewarddir)); 
trajvector = vectorfill(trajvector, -1, 1, findindex);
findex2 = findindex(find(~mod(trajvector(findindex),2)));
trajvector(findex2) = trajvector(findex2)-1;
if ismember(4,includeStates)
    statevector(findindex) = trajvector(findindex);
end
%otherwise assign the closest traj that happened in the past 
findindex = find((undefinedindex) & (trajcount > 1) & (backwarddir)); 
trajvector = vectorfill(trajvector, -1, -1, findindex);
findex2 = findindex(find(mod(trajvector(findindex),2)));
trajvector(findex2) = trajvector(findex2)+1;
if ismember(4,includeStates)
    statevector(findindex) = trajvector(findindex);
end

%LEVEL 5
% now we may have some undefined points in the beginning and end of the
% session- we assign these to the first trajectory
undefinedindex = (trajvector < 0); %still undefined times
trajvector(find(undefinedindex)) = -1;
statevector(find(undefinedindex)) = -1;
findindex = find((undefinedindex) & (forewarddir)); 
trajvector(findindex) = 1;
if ismember(5,includeStates)
    statevector(findindex) = trajvector(findindex);
end
findindex = find((undefinedindex) & (backwarddir)); 
trajvector(findindex) = 2;
if ismember(5,includeStates)
    statevector(findindex) = trajvector(findindex);
end

%LEVEL 6
% Fill in all remaining points, regardless of running velocity
undefinedindex = find(trajvector == -1); %still undefined times
trajvector = vectorfill(trajvector, -1, 0, undefinedindex);
if ismember(6,includeStates)
    statevector = trajvector;
end


state = ones(length(statematrix.time),1)*-1;
lindist = ones(length(statematrix.time),1)*-1;
state(validPointsIndex) = statevector;
lindist(validPointsIndex) = lineardistance;
    
    
    
    