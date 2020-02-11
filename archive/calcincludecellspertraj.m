function includecells = calcincludecellspertraj(minV, minPeak, trajgroup)
%gets all cells with peak > minpeak basedon occnormed firing at > minV
%assumes group is a structure and each element starts with the following columns: [animal day epoch tet
%cell ...], group{1} = [1 3 4 5 6 values of interest]
%
% minV -- minimun velocity, times when animal moves>min velocity will be included in analysis
% minPeak -- minimum peak, cells with peaks > min peak on specified
%   trajectory will be listed in include cells
% trajgroup -- strucutre of different trajectory requirements.
%   trajgroup{1}(1).include lists trajectories for animal 1 that can include a placefield --
%   at least one of these traj must include a placefeild,
%   trajgroup{1}(1).exclude lists trajectories for animal 1 that must not have any
%   placefields.  trajectories not listed in include or exclude will be
%   ignored.  If trajgroup is left empty, [], then finds all cells with
%   peak>minPeak on any traj (identical to calcincludecells)
%
% includecells -- list of cell indexes included based on criteria above,
% includecells{1} referst to trajgroup{1} above
%
% concatenates epochs, so that it finds cells with reuired placefields per
% day not per epoch in case animal sampled environment different on each
% day

if ~isempty(trajgroup)
    PKspertrajg = getpeakratespertraj(minV); % 7 columns: [an day epoch tetrode cell traj peakrate]
    allcells = unique(PKspertrajg{1}(:,1:5), 'rows'); %all cells: [an day epoch tet cell]

    for i = 1:length(trajgroup) %for each traj grouping
        includecells{i} = [];
        for j = 1:length(trajgroup{i}) %for each animal
            anPKspertrajg{1} = PKspertrajg{1}(PKspertrajg{1}(:,1)==j,:); %per animal
            includetraj = trajgroup{i}.include; %trajectories that should have a placefield
            includetraj = trajgroup{i}(j).include;
            excludetraj = trajgroup{i}(j).exclude;
            includes = anPKspertrajg{1}(ismember(anPKspertrajg{1}(:,6), includetraj) & anPKspertrajg{1}(:,7)>minPeak, :);
            excludes = anPKspertrajg{1}(ismember(anPKspertrajg{1}(:,6), excludetraj) & anPKspertrajg{1}(:,7)>minPeak, :);
            inclcells = unique(includes(:,[1 2 4 5]), 'rows'); %cells to include: [an day tet cell]
            exclcells = unique(excludes(:,[1 2 4 5]), 'rows'); %cells to exclude: [an day tet cell]

            cells = inclcells(~ismember(inclcells,exclcells,'rows'),:);
            %replace epoch info
            includecells{i} = [includecells{i}; allcells(ismember(allcells(:,[1 2 4 5]),cells,'rows'),:)];
        end
    end
else
    allrunPKgroup = getpeakratesallrun(minV);%Bar, Dwi, Cal
    includecells = allrunPKgroup{1}(allrunPKgroup{1}(:,6)>minPeak, 1:5);
end

end
