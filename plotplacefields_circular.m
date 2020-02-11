function out = plotplacefields_circular(placefields,dayindex,savedfiguresdir,binsize)
% SP 9.28.18
%this function plots place field firing rates around an annular track

if isempty(placefields)
    return
end

%% plot the units
for cellIdx = 1:size(placefields.ratenormocc)
    if ~isnan(placefields.ratenormocc(cellIdx,1))
        theta = 0:binsize:360;
        radius = 70:2:100;
        PFtoplot = [placefields.ratenormocc(cellIdx,:) placefields.ratenormocc(cellIdx,end)]; 
        FRtoplot = repmat(PFtoplot,16,1);
        
        figure; hold on;
        polarugcolor(radius,theta,FRtoplot);
        
        colormap(jet(256));
        brighten(0.3);
        title(['Place cell firing rate - F' num2str(dayindex(1)) num2str(dayindex(2)) ' - unit ' num2str(cellIdx)])
        animaldir = [savedfiguresdir 'F' num2str(dayindex(1)) '_' num2str(dayindex(2)) '\'];
        if ~exist(animaldir); mkdir(animaldir); end;
        filename = [animaldir 'placefields_annulartrack_clus' num2str(cellIdx)];
        saveas(gcf,filename,'png'); saveas(gcf,filename,'fig');
    end
end

end