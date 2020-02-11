function f = setfilteranimal(f, animaldata)
% f = setfilteranimal(f, animaldata)
% Sets the animal data for a data filter.  
% If ANIMALDATA is a string, then the function ANIMALDEF is called to
% translate the string to the 3 required data elements for the animal.  If
% ANIMALDATA is a cell array with multiple animalnames, then multiple
% filters are created, one for each animal.

if isstr(animaldata)
    animalinfo = animaldef(animaldata);
    f.animal = animalinfo;
    f.epochs = [];
    f.data = [];
    f.excludetime = [];
    f.output = [];
elseif iscell(animaldata)
    for i = 1:length(animaldata)
        animalinfo = animaldef(animaldata{i});
        f(i).animal = animalinfo;
        f(i).epochs = [];
        f(i).data = [];
        f(i). excludetime = [];
        f(i).output = [];
    end
end



