function [out] = getripples(index, ripples, cellinfo, varargin)
%
% out is [starttime endtime] of ripples
% only includes ripples if at least 50msec
%
% index is [day epoch]
%
% options are
%	'cellfilter', 'cellfilterstring'
%		     specifies a cell filter to select the tetrodes to use 
%	'minenergy', E
%		     specifies the minimum energy of a valid ripple event
%   'excludeperiods', excludeperiods
%   'minstd', S
%           specifies the minimum ripple threshold (in standand deviations)
%           to count as a ripple
%   'maxcell', 1
%           uses only the tetrode with the most cells on it
%   'minrip', n
%           requires that at least n tetrodes have a ripple
%
%

% assign the options
cellfilter = '';
minenergy = 0;
excludeperiods = [];
maxcell = 0;
minstd = 0;
minrip = 1;
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'cellfilter'
            cellfilter = varargin{option+1};
        case 'minenergy'
            minenergy = varargin{option+1};
        case 'excludeperiods'
            excludeperiods = varargin{option+1};
        case 'minstd'
            minstd = varargin{option+1};
        case 'maxcell'
            maxcell = varargin{option+1};
        case 'minrip'
            minrip = varargin{option+1};
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end
d = index(1);
e = index(2);
tetlist =  evaluatefilter(cellinfo{d}{e}, cellfilter);
tetlist = unique(tetlist(:,1))';
ripplestd = [];

if (maxcell)
    %find tetrode with max number of cells per day
    numcells = [0 0];
    for s = 1:length(tetlist)
        tet = tetlist(s); %tetrode
        numcells = [numcells; tet length(cellinfo{d}{e}{tet}) ];
    end
    [y i] = max(numcells(:,2));
    tetlist = numcells(i,1);


    r = ripples{d}{e}{tetlist(1)};
    validripples = find((r.energy >= minenergy) & (r.maxthresh >= minstd) & (~isExcluded(r.midtime,excludeperiods)));
    ripplestdout = r.maxthresh(validripples);
    out = [r.starttime(validripples) r.endtime(validripples)];
else

    % go through the tetlist and construct an an array where each element
    % represents the number of active tetrodes for each 1 ms timestep.

    r = ripples{d}{e}{tetlist(1)};
    times = r.timerange(1):0.001:r.timerange(end);
    nrip = zeros(size(times));
    ripplestd = zeros(size(times));
    for t = 1:length(tetlist)
        tmprip = ripples{d}{e}{tetlist(t)};
        %disp([min(tmprip.energy) median(tmprip.energy) max(tmprip.energy)]);

        % get the indeces for the ripples with energy above minenergy
        % and maxthresh above minstd
        rvalid = find((tmprip.energy >= minenergy) & (tmprip.maxthresh >= minstd) & (~isExcluded(tmprip.midtime,excludeperiods)));
        rtimes = [tmprip.starttime(rvalid) tmprip.endtime(rvalid)];
        tmpripplestd = [tmprip.maxthresh(rvalid) tmprip.maxthresh(rvalid)];
        % create another parallel vector with bordering times for zeros
        nrtimes = [(rtimes(:,1) - 0.00001) (rtimes(:,2) + 0.00001)];
        rtimes = reshape(rtimes', length(rtimes(:)), 1);
        rtimes(:,2) = 1;
        tmpriplestd = [rtimes(:,1) tmpripplestd(:)];
        nrtimes = [r.timerange(1) ; reshape(nrtimes', ...
            length(nrtimes(:)), 1) ; r.timerange(2)];
        nrtimes(:,2) = 0;
        % create a new list with all of the times in it
        tlist = sortrows([rtimes ; nrtimes]);
        [junk, ind] = unique(tlist(:,1));
        tlist = tlist(ind,:);
        
        %stdlist = sortrows([tmpripplestd ; nrtimes]);
        %[junk, ind2] = unique(stdlist(:,1));
        %stdlist = stdlist(ind2,:);
        
        % use interp to create a set of ones and zeros for each time
        % and add to nrip to get a cumulative count of the number of
        % ripples per timestep
        try
            nrip = nrip + interp1(tlist(:,1), tlist(:,2), times, 'nearest');
            %ripplestd = max([ripplestd interp1(stdlist(:,1), stdlist(:,2), times, 'nearest')]);
        catch
            keyboard
        end
    end

    %find the start and end borders of each ripple
    inripple = (nrip >= minrip);
    startrippleind = find(diff(inripple) == 1)+1;
    endrippleind = find(diff(inripple) == -1)+1;
    ripplestdout = [];
    if ((length(startrippleind) > 1) & (length(endrippleind) > 1))
        if (endrippleind(1) < startrippleind(1))
            endrippleind = endrippleind(2:end);
        end
        if (endrippleind(end) < startrippleind(end))
            startrippleind = startrippleind(1:end-1);
        end
        startripple = times(startrippleind);
        endripple = times(endrippleind);
%         for i = 1:length(startripple)
%             ripplestdout(i,1) = max(ripplestd(startrippleind(i):endrippleind(i)));
%         end
        
        out = [startripple(:) endripple(:)];
        
        out = out( ((out(:,2)-out(:,1))>.050),:); %only include if > 50msec long
    else
        out = [];
        ripplestdout = [];
    end
       
end


