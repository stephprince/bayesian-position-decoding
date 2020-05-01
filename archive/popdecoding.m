%Animal selection
%-----------------------------------------------------
%animals = {'Dudley','Miles','Conley','Ten','Bond','Frank','Nine'};
%animals = {'Dudley','Miles','Conley','Alex','Bond','Frank','Nine','Ten'};
%animals = {'Conley','Bond','Frank','Miles','Nine', 'Ten'};
%animals = {'Miles','Nine', 'Ten'};
%animals = {'Conley','Bond','Frank'};
%animals = {'Ten'};
%animals = {'Dudley','Miles','Conley', 'Ten'};
%animals = {'Bond','Frank','Nine'};
%animals = {'Miles'};
animals = {'Miles', 'Ten', 'Conley', 'Dudley', 'Eight', 'Five', 'Coriander'} ; %wtrack = 5

%-----------------------------------------------------


%Filter creation for training data, eg place field activity
%--------------------------------------------------------
%cellfilter = '(isequal($area, ''CA1'') && ($meanrate < 7))';
%cellfilter = '(($meanrate < 7))';
cellfilter =  ' ( (isequal($area, ''CA3'') | isequal($area, ''CA1'') ) && ($meanrate < 4) )';

maxd = 12;
for d = 1:maxd
    epochfilter{d} = ['(isequal($type, ''run''))  & ($exposure ==', num2str(d), ') '];
end
%epochfilter{1} = ['($dailyexposure == 1) & isequal($environment, ''TrackB'')'];
%epochfilter{2} = ['($dailyexposure == 1) & isequal($environment, ''TrackA'') & ($exposureday > 3)'];

%set up a time filter to calculate place fields
%timefilter = { {'getlinstate', '(($traj ~= -1) & (abs($velocity) >= 3))', 6},{'getriptimes','($nripples == 0)', [], 'cellfilter', '((isequal($area, ''CA3''))|(isequal($area, ''CA1'')))'} };
%timefilter = { {'getriptimes','($nripples == 0)', [], 'cellfilter', '(isequal($area, ''CA1''))'} };
timefilter = { {'getriptimes','($nripples == 0)', [], 'cellfilter', ...
    '((isequal($area, ''CA3''))|(isequal($area, ''CA1'')))'} };

% create the training filter for calculating the place fields
trainingfilter = createfilter('animal',animals,'epochs',epochfilter,'cells',cellfilter,'excludetime', timefilter);
%-----------------------------------------------------------


%create training data by calulating the linearized rates of all cells
%--------------------------------------------
iterator = 'multicellanal';
trainingfilter = setfilteriterator(trainingfilter,iterator);
trainingfilter = setfilterfunction(trainingfilter, 'calcpopulationlinfields', {'spikes','linpos'},2,3);
trainingfilter = runfilter(trainingfilter);
%-------------------------------------------------



%Filter creation for position decoding
%--------------------------------------------------------
%cellfilter = '(($meanrate < 7) & ($tetnum > 4))';
%cellfilter = '(($meanrate < 7))';
cellfilter =  ' ( (isequal($area, ''CA3'') | isequal($area, ''CA1'') ) && ($meanrate < 4) )';

maxd = 12;
for d = 1:maxd
    epochfilter{d} = ['(isequal($type, ''run'')) & ($exposure ==', num2str(d), ') '];
end
%epochfilter{1} = ['isequal($environment, ''TrackB'')'];
%epochfilter{1} = ['($dailyexposure == 1) & isequal($environment, ''TrackB'')'];
%epochfilter{1} = ['($dailyexposure == 2) & isequal($environment, ''TrackB'')'];
%epochfilter{1} = ['($dailyexposure == 1) & isequal($environment,
%''TrackA'') & ($exposureday > 3)'];

timefilter = {};
clear decodefilter
global decodefilter ;

decodefilter = createfilter('animal',animals,'epochs',epochfilter,'cells',cellfilter,'excludetime', timefilter);
%-----------------------------------------------------------



%Get the binned spike counts for all cells
%-----------------------------------------------------------
% the cell count thresh is the minimum number of cells per ripple
cellcountthresh = 2;
ripVeqn = '<1';
ripwelldist = 20;
iterator = 'multicellanal';
decodefilter = setfilteriterator(decodefilter,iterator);
decodefilter = setfilterfunction(decodefilter, 'getpopulationevents2', {'spikes','linpos','pos','ripples','cellinfo', 'task'},cellcountthresh, ripVeqn, ripwelldist);
%decodefilter = setfilterfunction(decodefilter, 'getpopulationevents_nolinpos', {'spikes','pos', 'ripples','cellinfo'},cellcountthresh);
%timebin = .2;
%decodefilter = setfilterfunction(decodefilter, 'getpopulationrates', {'spikes','linpos'},timebin);
decodefilter = runfilter(decodefilter);
%----------------------------------------------------------



%establish the list of animals and epochs to look at
%list = [1 1;1 2;3 2;3 3;3 4;3 5;3 6;3 7;3 8;3 9;3 10;3 11];
%list = [1 1;1 2;2 1;2 2;2 3;2 4;2 5;2 6;2 7;2 8;3 2;3 3;3 4;3 5;3 6;3 7;3 8;3 9;3 10;3 11];
%list = [1 1;1 2;3 3;3 4;3 5;3 6;3 7;3 8;3 9;3 10;3 11];
%list = [1 1;1 2;3 2;3 3;3 4];
%list = [4 1;5 1;6 1;6 2;6 3;6 4];
%list = [1 2];
%out = [];
%for i = 1:size(list,1)
% get the decoding statistics for each replay event
%   out = [out; calcepochreplaystats(trainingfilter, decodefilter)];
%out = [out; calcepochreplaystatsbothtrainers(list(i,:), trainingfilter, decodefilter)];
%end
out = calcepochreplaystats(trainingfilter, decodefilter);
%[1:an 2:day 3:epoch 4:group(exposure) 5:slope 6:R^2 7:p 8:immobiletime 9:numcellsactive? 10:in/correct 11:fut/past 12:activpastCP 13:passnum]
 %added later: 14-16:maxdist 17:pdfpkdist]
%---------------------------------------------------------


%% plot

pthresh = 0.05

for col = [17] %[5 11 14:17]
    if col == 5
        lbl = 'slope of regression';
        bins = -50:1:50
        ind1 = out(:,10)==0 & out(:,7)<pthresh; %incorrect, below pthres
        ind2 = out(:,10)==1  & out(:,7)<pthresh;
    elseif col == 11
        lbl = 'fut/past bias'
        bins = -1:0.01:1
        ind1 = out(:,10)==0 & out(:,12)==1 & out(:,7)<pthresh; %incorrect, code past CP
        ind2 = out(:,10)==1 & out(:,12)==1 & out(:,7)<pthresh;
    elseif ismember(col, [14:16])
        lbl = 'maxdist'
        bins = 0:1:200
        ind1 = out(:,10)==0 & out(:,col)>80 & out(:,col)<160 & out(:,7)<pthresh;
        ind2 = out(:,10)==1 & out(:,col)>80 & out(:,col)<160 & out(:,7)<pthresh;
    elseif ismember(col, [17])
        lbl = 'pdf peakdist'
        bins = 0:5:200
        ind1 = out(:,10)==0 & out(:,col)>80 & out(:,col)<160 & out(:,7)<pthresh;
        ind2 = out(:,10)==1 & out(:,col)>80 & out(:,col)<160 & out(:,7)<pthresh;
    end
    
    %     figure
    %     subplot(2,1,1)
    %     hist(out(ind1,col),bins)
    %     [h p] = ttest(out(ind1,col));
    %     title(['incorrect. diff from zero p', num2str(p)])
    %     subplot(2,1,2)
    %     hist(out(ind2,col),bins)
    %     xlabel(lbl)
    %     ylabel('count')
    %     [h p] = ttest(out(ind2,col));
    %     title(['correct. diff from zero p', num2str(p)])
    %
    data1 = out(ind1,col);
    data2 = out(ind2,col);
    rp = plotmeanhistcumsum(data1, data2, 1,lbl, 'inc v c', [], bins)
    
    plotcumsum(data1, data2, bins)
    xlabel(lbl)
    title(['p ', num2str(rp)])
    
    h1=histc(data1, bins);
    h2=histc(data2, bins);
    yl2 = max([max(h1/sum(h1)) max(h2/sum(h2))]);
    figure
    plot(bins, h1/sum(h1), 'k', 'linewidth', 3)
    hold on
    plot(bins, h2/sum(h2), 'r', 'linewidth', 3)
    set(gca, 'fontsize', 14)
    ylabel('fraction')
    xlabel(lbl)
    title(['p ', num2str(rp)])
    ylim([0 0.5])
    
    Ns = [size(data1,1) size(data2,1)]
    
end