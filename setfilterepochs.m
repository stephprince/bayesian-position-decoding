function f = setfilterepochs(f, filterInput)
% f = setfilterepochs(f, filterInput)
% Sets the wanted recording epochs for the data filter f.  filterInput is
% either a  cell array of filter strings (for multiple data groups), or 
% just one filter string.  Each filter string following the rules in EVALUATEFILTER.m
% Assumes that each animal's data folder has files named 'task##.mat' that
% contain cell structures with task information.


for an = 1:length(f)
    if isempty(f(an).animal)
        error(['You must define an animal for the filter before filtering the epochs'])
    end
    datadir = f(an).animal{2};
    animalprefix = f(an).animal{3};
    task = loaddatastruct(datadir,animalprefix,'task');
    f(an).epochs = [];
    f(an).data = [];
    f(an).excludetime = [];
    if iscell(filterInput) %if there are multiple filters in a cell array, create multiple epoch groups
        for j = 1:length(filterInput)
            if isstr(filterInput{j})
                f(an).epochs{j} = evaluatefilter(task,filterInput{j});
                f(an).data{j} = [];
                f(an).excludetime{j} = [];
                for k = 1:size(f(an).epochs{j},1)
                    f(an).data{j}{k} = [];
                    f(an).excludetime{j}{k} = [];
                end
            else
                error('Each cell in filterInput must contain a string');
            end
        end
    elseif isstr(filterInput) %if there is only one filter string, create just one epoch group
        f(an).epochs{1} = evaluatefilter(task,filterInput);
        f(an).data{1} = [];
        f(an).excludetime{1} = [];
        for k = 1:size(f(an).epochs{1},1)    
            f(an).data{1}{k} = [];
            f(an).excludetime{1}{k} = [];
        end
    else
        error('FILTERINPUT must either be a cell array or a string');
    end
end

