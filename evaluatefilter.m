function resultindex = evaluatefilter(cellvar, filterString)
%
% resultindex = evaluatefilter(cellvar, filterString)
%
% Return the indices of the cell structure that have fields that satisfy
% the given constraints. Each row of the output is the index into one cell.
%
% CELLVAR - any cell structure, for example spikes{}{}{}{}.####
% FILTERSTRING - A string containing an expression with the field variables
%               of CELLVAR.  Each field variable must be preceded with a
%               '$'.
%
% Example: index = evaluatefilter(cellinfo,'(($meanrate < 10) && isequal($area,''CA3''))')
%          Returns all indices of cellinfo where the 'area' field is equal to 'CA3' and the 'meanrate' field is less than 10.



newstruct = [];
[variables, expression] = parsefilterstring(filterString);
for i = 1:length(variables)
    varfetch = cellfetch(cellvar, variables{i});
    for j = 1:length(varfetch.values)
        
        newstruct = setfield(newstruct,{j},variables{i},varfetch.values{j});
    end
end
results = [];
for i = 1:length(newstruct)
    structVar = newstruct(i);
    eval(['tmpresult = ',expression,';']);
    if isempty(tmpresult)
        tmpresult = 0;
    end
    results(i) = tmpresult;
end
    
cellindex = cellfetch(cellvar, '');
resultindex = cellindex.index(find(results),:);
