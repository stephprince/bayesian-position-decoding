function [out] = calccorrecttraj(index, linpos, task)
% [out] = calccorrecttraj(index, linpos, task)
%
% assumes animal is performing an alternation task, with a middle well and
% 2, alternating outside wells, though it does not have to be a W-shaped
% track.  assumes task{d}{e}.wellseq exists and gives correct well order
% such that the middle number corresponds to the middle well, the first and 
% third ones correspond to outside wells,
% for instance [2 1 3] would mean the middle well is well 1, the outside
% wells are 2 and 3.  these well numbers are determined in createtaskstruct
%
%   OUT: nx3 matrix, each row is a trajectory.  column 1 index into state
%   or linpos.statematrix for start of each trajectory column 2 is sequence of wells
%   at end of trajectory (well animal is heading to), column 3 is 1 if trajecotry 
%   was correct, 0 if incorrect 

%define correct seq with wellseq from task struct
midwell= task{index(1)}{index(2)}.wellseq(2); %middle well
owell1 = task{index(1)}{index(2)}.wellseq(1); %1st outside well
owell2 = task{index(1)}{index(2)}.wellseq(3); %2nd outside well
%get seq of wells visited
wellExitEnter = linpos{index(1)}{index(2)}.statematrix.wellExitEnter;
wells = wellExitEnter(:,2);
wellchange = find(diff(wells))+1;
wellchange = [1;wellchange]; %vector of indexes where exit well and traj changes

wellSequence = wells(wellchange); %seq of wells visited
correct = 0; %tracks number correct traj
correctSequence = []; %tracks if each traj correct (1) or incorrect(0)
total = 0;%tracks total number traj counted
lastoutwell = []; %tracks last outward bound well

%makes a vector: correctSeq of same length as wellSeq that is 0 if traj is
%incorrect, 1 if traj is incorrect
t = 1; %first traj
if ( ((wellSequence(t) == owell1) || (wellSequence(t) == owell2)) || (wellSequence(t) == midwell) ) %is an outbound traj
    correct = correct+1;
    tmpEval = 1;
    total = total+1;
    correctSequence = [correctSequence;tmpEval];
else
    total = total+1;
    tmpEval = 0;
    correctSequence = [correctSequence;tmpEval];
end

for t = 2:length(wellSequence)
    tmpEval = 0;
    if ( ((wellSequence(t) == owell1) || (wellSequence(t) == owell2)) && (wellSequence(t-1) == midwell) ) %is an outbound traj
        if t==2
            correct = correct+1;
            tmpEval = 1;
            total = total+1;
            correctSequence = [correctSequence;tmpEval];
        elseif ( (wellSequence(t-2) ~= wellSequence(t)) ) %if last outward bound well is different from current outward bound well
            correct = correct+1;
            tmpEval = 1;
            total = total+1;
            correctSequence = [correctSequence;tmpEval];
        else
            total = total+1;
            correctSequence = [correctSequence;tmpEval];
        end
    elseif ( (wellSequence(t) == midwell) & (wellSequence(t-1) ~= midwell) ) %inward traj
        correct = correct+1;
        tmpEval = 1;
        total = total+1;
        correctSequence = [correctSequence;tmpEval];
    else %incorrect
        total = total+1;
        correctSequence =[correctSequence; tmpEval];
    end

    if ( ((wellSequence(t) ==owell1) || (wellSequence(t) == owell2)) )
        lastoutwell = wellSequence(t); %set last outward well to current outwell for next trial
    end
end

performance = correct/(length(wellSequence)-1);

seqMatrix= [wellchange wellSequence correctSequence]; %index from linpos, well entering, in/correct
seqMatrix = [seqMatrix; length(wells) seqMatrix(end, 2:3)]; %to include last trajectory
out = seqMatrix;
end