function plotplacefields_sortedpeakfr(placefields,dayindex,savedfiguresdir)
%SP 9.27.18
%this function plots all place fields in order of peak FR in order to see
%the span of the place fields and properties across the track

if isempty(placefields)
    return
end

%% find and sort by peak firing rate
%find peak firing rate
for cellIdx = 1:size(placefields.ratenormocc,1)
    if ~isnan(placefields.ratenormocc(cellIdx,1))
        maxind(cellIdx) = find(placefields.ratenormocc(cellIdx,:) == max(placefields.ratenormocc(cellIdx,:)));
    end
end

%sort units and normalize by peak firing rate so max = 1
if ~exist('maxind'); maxind = []; end
[sorted sortedind] = sort(maxind);
sortedind(sorted == 0) = [];
for cellIdx = 1:length(sortedind)
    placefields.sortedratenormocc(cellIdx,:) = placefields.ratenormocc(sortedind(cellIdx),:);
    placefields.sortedratenormocc_normfr(cellIdx,:) = placefields.ratenormocc(sortedind(cellIdx),:)./max(placefields.ratenormocc(sortedind(cellIdx),:));
end

%% plot the units
if isfield(placefields,'sortedratenormocc_normfr')
    %normalized by peak fr
    figure; hold on;
    imagesc(placefields.sortedratenormocc_normfr)
    title(['Place cells - F' num2str(dayindex(1)) num2str(dayindex(2))])
    animaldir = [savedfiguresdir 'F' num2str(dayindex(1)) '_' num2str(dayindex(2)) '\'];
    if ~exist(animaldir); mkdir(animaldir); end;
    filename = [animaldir 'placefields_sortedpeakfr'];
    saveas(gcf,filename,'png'); saveas(gcf,filename,'fig');
    
   %non-normalized FR to get idea of distribution of units
   figure; hold on;
    imagesc(placefields.sortedratenormocc)
    title(['Place cells raw FR - F' num2str(dayindex(1)) num2str(dayindex(2))])
    colorbar
    animaldir = [savedfiguresdir 'F' num2str(dayindex(1)) '_' num2str(dayindex(2)) '\'];
    if ~exist(animaldir); mkdir(animaldir); end;
    filename = [animaldir 'placefields_sortedpeakfr_rawFR'];
    saveas(gcf,filename,'png'); saveas(gcf,filename,'fig');
end

end