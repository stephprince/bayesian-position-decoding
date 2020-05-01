function out = applyplacefieldcriteria(spikeratenormocc, dayindex, behavioridx, sessindex, processeddatadir, PFbinsize, datadir, varargin)
% SP 9.27.18
% this function takes all single units from a recording and determines
% whether to include them in place field calculations depending on various
% criteria

% INPUTS:
% spikeratenormocc - firing rate x space normalized for time occupancy
% dayindex - [animal day rec]

% OPTIONS:
% pyr - includes only units that are classified as PYR cells
%     - input is the celltypedirectory
% peak - includes only units that have a peak firing rate above threshold
%     - input is the threshold
% std - includes only units that have a peak FR some number of standard deviations above the mean
%     - input is the threshold
% reliability - includes only units with a reliability index above threshold (see Saleem et al. 2013 or 2018)
%     - input is the threshold
% maxmean - includes only units that have a mean firing rate below some threshold (to get rid of potential interneurons)
%     - input is the threshold
% spatialinfo - includes only units that have spatial info above threshold
%     - input is the threshol

%% import cases
applycelltype = 0;
applypeak = 0;
applystd = 0;
applyreliability_stabletimes = 0;
applyreliability_alltimes = 0;
applymaxmean = 0;
applyspatialinfo = 0;
for option = 1:2:length(varargin)-1
    if isstr(varargin{option})
        switch(varargin{option})
            case 'pyr'
                applycelltype = 1;
                celltypedir = varargin{option+1};
            case 'peak'
                applypeak = 1;
                peakthreshold = varargin{option+1};
            case 'std'
                applystd = 1;
                stdthreshold = varargin{option+1};
            case 'reliability_alltimes'
                applyreliability_alltimes = 1;
                reliabilitythreshold_alltimes = varargin{option+1};
            case 'reliability_stabletimes'
                applyreliability_stabletimes = 1;
                reliabilitythreshold_stabletimes = varargin{option+1};
            case 'maxmean'
                applymaxmean = 1;
                maxmeanthreshold = varargin{option+1};
            case 'spatialinfo'
                applyspatialinfo = 1;
                spatialinfothreshold = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    else
        error('Options must be strings, followed by the variable');
    end
end

if isempty(spikeratenormocc)
    out = [];
    return
end

%% get exclusion or inclusion indices
%celltype
if applycelltype
    filename = [celltypedir 'cellTypeProps_F' num2str(dayindex(1)) '_' num2str(dayindex(2)) '.mat'];
    load(filename)
    cellTypeInfo = allrecProps{dayindex(1)}{dayindex(2)}.cellTypeInfo;
    tempexclcelltype = ones(1,size(spikeratenormocc.rate,1));
    tempexclcelltype(cellTypeInfo.PYRidx) = 0;
    exclcelltype = find(tempexclcelltype == 1);
else
    exclcelltype = [];
end
out.exclcelltype = exclcelltype;

%peak firing 
if applypeak
    peakrates = max(spikeratenormocc.ratenormocc,[],2); %SP changed 191014 not sure what happened here 
    exclpeak = find(peakrates < peakthreshold);     
else
    exclpeak = [];
    peakthreshold = [];
end
out.exclpeakidx = exclpeak;
out.peakthreshold = peakthreshold;

%maxmean
if applymaxmean
    meanrates = nanmean(spikeratenormocc.ratenormocc,2);
    exclmaxmean = find(meanrates > maxmeanthreshold);     
else
    exclmaxmean = [];
    maxmeanthreshold = [];
end
out.exclmaxmeanidx = exclmaxmean;
out.maxmeanthreshold = maxmeanthreshold;

%std
if applystd
    sd1 = nanstd(spikeratenormocc.ratenormocc,[],2);
    mean = nanmean(spikeratenormocc.ratenormocc,2);
    exclstd = find(peakrates < mean+(stdthreshold*sd1)); 
else
    exclstd = [];
    stdthreshold = [];
end
out.exclstdidx = exclstd;
out.stdthreshold = stdthreshold;

%reliability
if applyreliability_alltimes
    filename2 = [datadir 'reliabilityindex_alltimes.mat'];
    if exist(filename2)
        load(filename2)
    else
        reliabilityindex = getreliabilityindex_alltimes(spikeratenormocc, dayindex, behavioridx, reliabilitythreshold_alltimes, processeddatadir);
        save(filename2,'reliabilityindex');
    end
    excl1 = find(reliabilityindex < reliabilitythreshold_alltimes);
    excl2 = find(isnan(reliabilityindex));
    exclreliability_alltimes = unique([excl1 excl2]);
else
    exclreliability_alltimes = [];
    reliabilitythreshold_alltimes = [];
end
out.exclreliabilityidx_alltimes = exclreliability_alltimes;
out.reliabiltythreshold_alltimes = reliabilitythreshold_alltimes;

if applyreliability_stabletimes
    filename2 = [datadir 'reliabilityindex_stabletimes.mat'];
    if exist(filename2)
        load(filename2)
    else
        reliabilityindex = getreliabilityindex_stabletimes(spikeratenormocc, PFbinsize, dayindex, behavioridx, sessindex, reliabilitythreshold_stabletimes, processeddatadir);
        save(filename2,'reliabilityindex');
    end
    excl1 = find(reliabilityindex < reliabilitythreshold_stabletimes);
    excl2 = find(isnan(reliabilityindex));
    exclreliability_stabletimes = unique([excl1 excl2]);
else
    exclreliability_stabletimes = [];
    reliabilitythreshold_stabletimes = [];
end
out.exclreliabilityidx_stabletimes = exclreliability_stabletimes;
out.reliabiltythreshold_stabletimes = reliabilitythreshold_stabletimes;

%spatial info percentile
if applyspatialinfo
    filename2 = [datadir 'spatialinfo_' num2str(spatialinfothreshold) 'threshold.mat'];
    if exist(filename2)
        load(filename2)
    else
        [spatialinfo spatialinfoindex] = getspatialinfo(spikeratenormocc, PFbinsize, dayindex, behavioridx, sessindex, spatialinfothreshold, processeddatadir);
        save(filename2,'spatialinfo','spatialinfoindex');
    end
    excl1 = find(spatialinfoindex < 1); %spatial info index is a logical of units that crossed the threshold already, so find ones that didint 
    excl2 = find(isnan(spatialinfoindex));
    exclspatialinfo = unique([excl1 excl2]);
else
    exclspatialinfo = [];
    spatialinfothreshold = [];
end
out.exclspatialinfo = exclspatialinfo;
out.spatialinfothreshold = spatialinfothreshold;

%% combine multiple criteria
excl = unique([exclpeak' exclstd' exclreliability_alltimes exclreliability_stabletimes exclcelltype exclmaxmean' exclspatialinfo]);

%% apply criteria
fnames = fieldnames(spikeratenormocc);
for fieldIdx = 1:3
    temp = spikeratenormocc.(fnames{fieldIdx});
    temp(excl,:) = nan; %replace excl indices with nans
    out.(fnames{fieldIdx}) = temp;
    exclAll = find(isnan(temp(:,1)));
end
out.stabletimesall = spikeratenormocc.stabletimesall;
out.exclallcriteriaidx = excl;

end