
tmp = zscore(datafiles{1,1}')';
load('roi.mat');
ROI = cell(4,1);
ROI{1}=roi.hmtL';
ROI{2}=roi.hmtR';
ROI{3}=roi.mstL';
ROI{4}=roi.mstR';

%%
block = 12;
close all
figure(1); clf
for iplot = 1:4
    subplot(4,1,iplot)
    hold on
    for TR = 1:size(tmp,2)
        betas = tmp(ROI{iplot},TR);
        plot(TR*ones(size(betas)),betas,'.','Color',[0.5 0.5 0.5])
        plot(TR,mean(betas),'r.','markersize',20)
        ylim([-3 3])
        xlim([1 size(tmp,2)])
        %drawnow
    end
    plot(1:size(tmp,2),mean(tmp(ROI{iplot},:)),'r-','linewidth',2)
    plot([block*2:block*2:size(tmp,2);block*2:block*2:size(tmp,2)],[-3 3],'linewidth',2,'Color',[0 0 0])
end

%% average
figure(2); clf
for iplot = 1:4
    subplot(4,1,iplot)
    hold on
    meanB = zeros(block*6,1);
    for TR = 1:block*6
        betas = mean(tmp(ROI{iplot},TR:block*6:290),2);
        meanB(TR) = mean(betas);
        plot(TR*ones(size(betas)),betas,'.','Color',[0.5 0.5 0.5])
        plot(TR,meanB(TR),'r.','markersize',20)
        ylim([-3 3])
        xlim([1 72])
    end
    plot(1:block*6,meanB,'r-','linewidth',2)
    plot([block*2:block*2:block*6;block*2:block*2:block*6],[-3 3],'linewidth',2,'Color',[0 0 0])
end


%%
figure(2); clf

for TR = 72:120
betas = tmp(:,TR);

betas(betas==0) = -50;
betas(isnan(betas)) = -50;



bins = -3:0.01:3
cmap0 = cmaplookup(bins,min(bins),max(bins),[],(cmapsign4));
[rawimg,Lookup,rgbimg] = cvnlookup(subj,1,betas,[min(bins) max(bins)],cmap0,min(bins),[],0,{'roiname',{'MT_exvivo'},'roicolor',{'w'},'drawrpoinames',0,'roiwidth',{5},'fontsize',20});

color = [0.5];
[r,c,t] = size(rgbimg);
[i j] = find(all(rgbimg == repmat(color,r,c,3),3));

for ii = 1: length(i)
    rgbimg(i(ii),j(ii),:) = ones(1,3);
end

a = imagesc(rgbimg); axis image tight;

    
set(gcf,'Position',[ 277         119        1141         898])
axis off
hold on
% subplot(2,1,2)
plot(0,0);
colormap(cmap0);
hcb=colorbar('SouthOutside');
hcb.Ticks = [0 1];
hcb.TickLabels = {'-0.5';'0.5'}
hcb.FontSize = 25
hcb.Label.String = ['Moving dots vs. stationary (' num2str(TR) ' TR)']
hcb.TickLength = 0.001;
drawnow

end