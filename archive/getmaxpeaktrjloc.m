function maxpkmatrix = getmaxpeaktrjloc(trajdata)
%[trj dist pk ]= getmaxpeaktrjloc(trajdata)
% INPUT
%   trajdata{trj} is 7 columns, 5th column is occ norm's firing rate
% OUTPUT
%   maxpkmatrix is [trj dist pk]
%   trj is trajectory max peak found on
%   dist is distance along that peak
%   pk is peak rate 

maxpk = [];
for t = 1:length(trajdata)
    if ~isempty(trajdata{t})
    [y ind] = max(trajdata{t}(:,5));
    maxpk = [maxpk; t trajdata{t}(ind,1) y];
    else
       maxpk = [maxpk; t NaN NaN];
    end
end

maxpkmatrix = maxpk;