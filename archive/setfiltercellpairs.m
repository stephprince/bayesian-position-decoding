function f = setfiltercellspairs(f, filterStringArray)
% f = setfiltercells(f, filterStringArray)
% For each epoch in the filter F, this function finds the indices to the
% cells that satisfy each of the filter condtions in the 2nd and 3rd
% element of the filterStringArray.  
%
% This creates a two lists of cells, and the first element of 
% filterStringArray specifies the combination operator to use to create the final
% list of cells.  
%
% Valid operators are:
%	'allcomb'  - specifies that the resulting list contains all
%		possible combinations where each combination is one cell from
%		each list of cells.  Only combinates where each cell is listed
%		at most once are included.
%	'difftet'  - same as 'allcomb' but only cells from different tetrodes
%		are included.
%
% The final list is a double array with one row for each pair and six columns:
% dataset epoch tet1 cell1 tet2 cell2
%
% The syntax for each filterString is defined in EVALUATEFILTER.m. 
% The animal and desired epochs for the filter need to be predefined. 
% Assumes that each animal's data  folder contains a file 'cellinfo.mat' that 
% contains a cell structure with information about each cell.

if (length(filterStringArray) ~= 3) 
    error('filterStringArray must have 3 elements');
end

for an = 1:length(f)
    if isempty(f(an).animal)
        error(['You must define an animal for the filter before filtering the cells'])
    end
    if isempty(f(an).epochs)
        error(['You must define the desired epochs for the filter before filtering the cells'])
    end
    datadir = f(an).animal{2};
    animalprefix = f(an).animal{3};
    cellinfo = loaddatastruct(datadir,animalprefix,'cellinfo');

    for i = 1:length(f(an).epochs)
        if isempty(f(an).epochs{i})
            f(an).data{i} = [];
        end
        for j = 1:size(f(an).epochs{i},1)
	    for k = 2:length(filterStringArray)
		tmplist{k-1} = evaluatefilter(cellinfo{f(an).epochs{i}(j,1)}{f(an).epochs{i}(j,2)}, filterStringArray{k});
		tmpind{k-1} = 1:size(tmplist{k-1}, 1);
	    end
	    % combine the lists using the option in the first element of the
	    % filterStringArray
	    [i1 i2] = meshgrid(tmpind{1}, tmpind{2});
	    celllist = [tmplist{1}(i1,:) tmplist{2}(i2,:)];
	    switch filterStringArray{1}
	    case 'allcomb'
		% get rid of identical elements
		valid = (sum(celllist(:,1:2) == celllist(:,3:4), 2) ~= 2);
		celllist = celllist(valid,:);
	    case 'difftet'
		% get rid of elements with the same tetrode number
		valid = (celllist(:,1) ~= celllist(:,3));
		celllist = celllist(valid,:);
	    end
	    if (length(celllist))
		% check for pairs listed twice because the order was reversed
		tmpcelllist = [celllist(:,3:4) celllist(:,1:2)];
		[c i1 i2] = intersect(celllist, tmpcelllist, 'rows');
		valid = ones(size(celllist,1), 1);
		% for each member of tmpcelllist that is in celllist, we
		% need to delete the correponding member of celllist if
		% it's index is larger than the index in tmpcelllist
		valid(i2) = (i2 < i1);
		f(an).data{i}{j} = celllist(logical(valid),:);
	    else
		f(an).data{i}{j} = [];
	    end
        end
    end
end
