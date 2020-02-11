function f = setfiltercells(f, filterString)
% f = setfiltercells(f, filterString)
% For each epoch in the filter F, this function finds the indices to the
% cells that satisfy the given filter condtions in filterString.  The syntax
% for filterString is defined in EVALUATEFILTER.m. The animal and desired epochs
% for the filter need to be predefined. Assumes that each animal's data
% folder contains a file 'cellinfo.mat' that contains a cell structure with
% information about each cell.

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
            f(an).data{i}{j} = evaluatefilter(cellinfo{f(an).epochs{i}(j,1)}{f(an).epochs{i}(j,2)}, filterString);
        end
    end
end