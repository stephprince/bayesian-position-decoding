function out = calclinvelocity(linpos,index, varargin)
% out = calclinvelocity(linpos,index, varargin)
% Produces a cell structure with the fields:
% time, velocity (linearized)
%   INDEX - N by 2 vector [day epoch]
%
%   OPTION: 'smooth', default no smoothing
%                   to compute linear speed, can smooth linear position data.
%                   It is smoothed with a gaussian of length VSW and std VSW/4.
%                   default for lineardayprocess is 2 seconds

smooth = [];
if ~isempty(varargin)
    smooth = varargin{2};
end

out.time = linpos{index(1)}{index(2)}.statematrix.time;
timestep = mean(diff(linpos{index(1)}{index(2)}.statematrix.time));
if isempty(smooth)
    linv=diff([0; linpos{index(1)}{index(2)}.statematrix.lindist])/timestep;
    linv(linpos{index(1)}{index(2)}.statematrix.lindist==-1)=NaN;
elseif ~isempty(smooth)
    %to smooth
    npoints = smoothwidth/timestep;
    filtstd = smoothwidth/(4*timestep);
    % the default filter for smoothing motion direction is a n second long gaussian
    % with a n/4 second stdev
    filt = gaussian(filtstd, npoints);
    % smooth the linear distances with the filter and then go through the linear
    % distance positions and take all of the differences to get velocities
    smoothdist = smoothvect(linpos{index(1)}{index(2)}.statematrix.lindist, filt);
    linv = diff(smoothdist) / timestep;
end
    out.velocity=linv;

end