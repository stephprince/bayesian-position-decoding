
function output = cellfetch(cellinput, field, branchindex)
% searchresults = cellfetch(input, field)
%
% This function parces a nested cell variable 'input', and assumes nothing about
% the structure of the internal cell tree. It returns a structure with the value of the
% field variable in the field 'field' at the end of each cell branch.  The
% indices of each entry is put in the index field of the output.

if ~isstr(field)
    error('FIELD must be a string');
end
output.index = [];
output.values = {};
indexcount = 0;
if (nargin < 3)
    branchindex = [];
end

if iscell(cellinput)
    for i = 1:length(cellinput)
        if isstruct(cellinput{i})
            indexcount = indexcount + 1;
            output.index(indexcount,1) = i;
            if isfield(cellinput{i},field)
                output.values{indexcount,1} = getfield(cellinput{i},field);
            else
                output.values{indexcount,1} = [];
            end
        elseif iscell(cellinput{i})
            [branch] = cellfetch(cellinput{i}, field, i);
            numentries = size(branch.index,1);
            branch.index = [repmat(i,[numentries,1]) branch.index];
            output.index(indexcount+1:indexcount+numentries,1:size(branch.index,2)) = branch.index;
            indexcount = indexcount+numentries;
            output.values = [output.values; branch.values];
        end
    end
end
            
