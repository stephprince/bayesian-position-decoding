function [out] = getripcodefutpast( index, linpos, spikes, riptimes, traj)

%get state to compute trajdata
[state, lindist] = getbehavestate(linpos, index(1,1), index(1,2), 6, 'minlinvelocity', 2);  %filters linpos.statematrix.traj by includestates

% trajdata & Peak PF info
pktrjs = zeros(1, size(index,1)); pkwelldist = pktrjs;
for c = 1:size(index,1) %for each cells
    if ~isempty(spikes{index(c,1)}{index(c,2)}{index(c,3)}{index(c,4)}) && ~isempty(spikes{index(c,1)}{index(c,2)}{index(c,3)}{index(c,4)}.data)
        
        %get occ normd firing rate
        trajdata{index(c,1)}{index(c,2)}{index(c,3)}{index(c,4)} = calclinfields(spikes,state,lindist,linpos,index(c,:));
        
        %find peak and distance from center well to peak for each cell in
        %pair
        pkmatrix1 = getmaxpeaktrjloc(trajdata{index(1,1)}{index(1,2)}{index(c,3)}{index(c,4)}); %[trj dist maxpk]
        [pk1 ind] = max(pkmatrix1(:,3));
        peakrate(c) = pk1;
        
        %get traj and location of PF peak for each cell
        pktrjs(c) = pkmatrix1(ind,1); %traj peak is on
        pkwelldist(c) = pkmatrix1(ind,2); %from center well
    end
end

%for each cell compute if active in each ripple and get trajdata
ripspikes = zeros(size(index,1), size(riptimes));
for c = 1:size(index,1) %for each cells
    if ~isempty(spikes{index(c,1)}{index(c,2)}{index(c,3)}{index(c,4)}) && ~isempty(spikes{index(c,1)}{index(c,2)}{index(c,3)}{index(c,4)}.data)
        %get activation in ripples
        ripactive = unique(getbinindex(spikes{index(c,1)}{index(c,2)}{index(c,3)}{index(c,4)}.data(:,1),[riptimes] )); %find all bins where there was a spikes
        ripactive = ripactive(ripactive>0); %exclude 0 bin index
        ripspikes(c, ripactive) = 1; %if spike in a bin for this cell, update ripspikes, each row is a cell, each column is a ripple
    end
end

%for each ripple compute if agrees with future or past traj
pkrip = repmat(pktrjs',1, size(ripspikes,2)) .* logical(ripspikes); %traj for each cell that fired in the ripple, each column is a ripple
pkrip = repmat(pkwelldist'>80,1, size(ripspikes,2)) .* pkrip; %only include if cells past CP
rip1 = sum( pkrip == 1 | pkrip == 2 ); % number of cells that have peaks on trj 1 or 2
rip3 = sum( pkrip == 3 | pkrip == 4 ); % number of cells that have peaks on trj 3 or 4
ripspktrj = zeros(1, size(ripspikes,2)); ripspktrjagree = ripspktrj;
ripspktrj(rip1>rip3 ) = 1 ; %if more cells code for trj 1, then assign ripple to coding for trj 1
ripspktrj(rip1<rip3 ) = 3;

ripspktrjagree(traj' == ripspktrj) = 1; %if traj animal is about to take equals traj most coded by cells in ripple
ripspktrjagree(traj' ~= ripspktrj & ripspktrj~=0 ) = -1; %agrees with past traj

out = ripspktrjagree';

end