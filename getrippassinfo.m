function [traj correct passnum] = getrippassinfo(index, linpos, task, riptimes)
%[traj correct] = getrippassinfo(index, linpos, task, riptimes)
%
% riptimes is Nx2 [starttimes endtimes] of each ripple
%
% traj is traj for each ripple of riptimes
% correct is correct/incorrect for each ripple of riptimes

%get state to pass info including traj and lindistance
[state, lindist] = getbehavestate(linpos, index(1,1), index(1,2), 6, 'minlinvelocity', 2);  %filters linpos.statematrix.traj by includestates

% get correct traj
if ~isfield(task{index(1,1)}{index(1,2)}, 'wellseq')
    [seqMatrix] = calccorrecttraj(index(1,[1 2]), linpos, task, 'correctorder', [2 1 3]);
else
    [seqMatrix] = calccorrecttraj(index(1,[1 2]), linpos, task);
end
statevector = getcorrectstate(seqMatrix, state, 1,0);

%get pass info
wellchange = seqMatrix(:,1);
startind = [wellchange(1:end-1)]; %vector of indexes where exit well and traj changes
endind = [wellchange(2:end-1)-1; wellchange(end)];
passtimes = [linpos{index(1,1)}{index(1,2)}.statematrix.time(startind) linpos{index(1,1)}{index(1,2)}.statematrix.time(endind)];
passcorrect = seqMatrix(1:end-1,3);

%calc trajectory for each pass
wellvis = linpos{index(1,1)}{index(1,2)}.statematrix.wellExitEnter(startind,:);  %start and end wells visited, each row is a pass
passtraj = zeros(size(wellvis,1), 1);
trajmat = linpos{index(1,1)}{index(1,2)}.trajwells;
trajnum = reshape([1:size(trajmat,1)*2], 2,size(trajmat,1))';
for t = 1:size(trajmat,1) %for each traj and its reverse
    passtraj( ismember(wellvis,trajmat(t,:),'rows') ) = trajnum(t,1);
    passtraj( ismember(wellvis,fliplr(trajmat(t,:)),'rows') ) = trajnum(t,2);  %trajectory for each pass
end
passinfo = [passtimes passcorrect passtraj]; %passstarttime passendtime in/correct traj


%for each ripple get traj
currentpasstrj = []; nextpasstrj = []; trajc = zeros(size(riptimes,1),3);
for n = 1:size(riptimes,1);% get wellexitenter that ripple occured on
    % get info on next and last pass
    currentp = lookup(riptimes(n,1), passtimes(:,1), -1);%findcurrent pass row, where start of ripple occurs
    currentpasstrj = [currentpasstrj; passinfo(currentp,3:4) currentp];
    if currentp ~= size(passinfo,1)
        nextpasstrj = [nextpasstrj ; passinfo(currentp+1,3:4) currentp+1]; %[passcorrect passtraj passnum];
    else
        nextpasstrj = [nextpasstrj ; NaN NaN NaN];
    end
end
%sort rippsikes by current/future traj
eventraj = rem(currentpasstrj(:,2),2)==0; %if it is inbound
oddtraj = rem(currentpasstrj(:,2),2)~=0; %outbound
trajc(eventraj,:) = nextpasstrj(eventraj,:); %if inbound traj, select next pass
trajc(oddtraj,:) = currentpasstrj(oddtraj,:); %if outbound traj, select current pass

traj = trajc(:,2);
correct = trajc(:,1);
passnum = trajc(:,3);
end