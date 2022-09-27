clc
clear
close all
addpath(genpath('/Users/jk7127/matlab/toolboxes/GLMDenoise'));
addpath(genpath('/Users/jk7127/matlab/toolboxes/cvncode'));
addpath(genpath('/Users/jk7127/matlab/toolboxes/knkutils'));
addpath(genpath('/Applications/freesurfer/matlab/'))

datadir = '/Volumes/server/Projects/akinetopsia/derivatives/fmriprep'
subjects = dir(sprintf('%s/sub*',datadir));
subjects = subjects([subjects.isdir]);


hemi = {'L';'R'}
load dms
%%
s = 2
surf = 'fsnative'
% to convert from gii to mgz use this -  > mri_convert yourfile.gii yourfile.mgz
subj = subjects(s).name;
subj_dir = sprintf('%s/%s/ses-nyu3t01/func',datadir,subj);
d_L = dir(sprintf('%s/*mt*%s_hemi-%s*.mgz',subj_dir,surf,hemi{1}))
d_R = dir(sprintf('%s/*mt*%s_hemi-%s*.mgz',subj_dir,surf,hemi{2}))

files2run = [1:length(d_L)];

datafiles = cell(1,length(files2run));
matrices = cell(1,length(files2run));


for runs = files2run
    
    
    tmp = MRIread(sprintf('%s/%s',subj_dir,d_L(runs).name));
    data_L = squeeze(tmp.vol);

%     if contains('fsaverage',surf)
%         [surfvals_LH] = cvnsurfsmooth('fsaverage',data_L,fwhm,'lh','inflated','','iterative');
%         l = read_label('fsaverage',sprintf('lh.%s',label2mask));
%         
%     else
%         
%         [surfvals_LH] = cvnsurfsmooth(subj,data_L,fwhm,'lh','inflated','','iterative');
%         l = read_label(subj,sprintf('lh.%s',label2mask));
%         
%     end
%     
%     roi =  l(:,1) + 1;
%     roi_mask = zeros(size(surfvals_LH,1),1);
%     roi_mask(roi) = 1;
% %     surfvals_LH = surfvals_LH .* double(roi_mask);
% %     surfvals_LH = surfvals_LH;

    tmp = MRIread(sprintf('%s/%s',subj_dir,d_R(runs).name));
    data_R = squeeze(tmp.vol);
    
%     if contains('fsaverage',surf)
%         [surfvals_RH] = cvnsurfsmooth('fsaverage',data_R,fwhm,'rh','inflated','','iterative');
%         l = read_label('fsaverage',sprintf('rh.%s',label2mask));
%         
%     else
%         [surfvals_RH] = cvnsurfsmooth(subj,data_R,fwhm,'rh','inflated','','iterative');
%         l = read_label(subj,sprintf('rh.%s',label2mask));
%         
%     end
%     roi =  l(:,1) + 1;
%     roi_mask = zeros(size(surfvals_RH,1),1);
%     roi_mask(roi) = 1;
% %     surfvals_RH = surfvals_RH .* double(roi_mask);
% %     surfvals_RH = surfvals_RH;
    datafiles{runs} = cat(1,data_L,data_R);
    
    
    
    
end
datafiles(cellfun(@isempty,datafiles)) = [];


%% remove first volume

for runs = files2run
    
    matrices{runs} = (dms{runs});
    
end

%%
% opt = struct('numboots',100,'numpcstotry',5);  % to use denoising, could set numpcstotry to 20
% myresults = GLMdenoisedata(matrices,datafiles,2,3.4,[],[]);
% myresults = GLMestimatemodel(matrices,datafiles,1,1,'assume',[],0)
myresults = GLMdenoisedata(matrices,datafiles,1,1);

%%
setenv('SUBJECTS_DIR','/Volumes/server/Projects/akinetopsia/derivatives/freesurfer')

conditions = {'central_moving';'central_stationary';'left_moving';'left_stationary';'right_moving';'right_stationary'}
subject = subj
bidsfolder = '/Volumes/server/Projects/akinetopsia/'
resultsdir = sprintf('%s/derivatives/GLMdenoise/%s/ses-nyu3t01/',bidsfolder,subject)
fspth = fullfile(bidsfolder, 'derivatives', 'freesurfer', [subject]);
lcurv = read_curv(fullfile(fspth, 'surf', 'lh.curv'));
rcurv = read_curv(fullfile(fspth, 'surf', 'rh.curv'));
mgz = MRIread(fullfile(fspth, 'mri', 'orig.mgz'));
leftidx  = 1:numel(lcurv);
rightidx = (1:numel(rcurv))+numel(lcurv);
%%

% findnan = isnan(myresults.modelmd{2}(:));
% myresults.modelmd{2}(findnan) = 0;

assert(isequal(numel(lcurv) + numel(rcurv), numel(myresults.R2)), ...
        'The number of vertices in the aprf results and the l&r curv files do not match;');



for b = 1 : size(myresults.modelmd{2},2)
    
mgz.vol = myresults.modelmd{2}(leftidx,b);
MRIwrite(mgz, fullfile(resultsdir, sprintf('lh.%s.mgz',conditions{b})));
mgz.vol = myresults.modelmd{2}(rightidx,b);
MRIwrite(mgz, fullfile(resultsdir, sprintf('rh.%s.mgz',conditions{b})));
end

% findnan = isnan(myresults.R2(:));
% myresults.R2(findnan)= 0;
mgz.vol = double(myresults.R2(leftidx,:));
MRIwrite(mgz, fullfile(resultsdir, sprintf('lh.%s.mgz','vexpl_glm')));
mgz.vol = double(myresults.R2(rightidx,:));
MRIwrite(mgz, fullfile(resultsdir, sprintf('rh.%s.mgz','vexpl_glm')));


mgz.vol = double(myresults.R2(leftidx,:))>10;
MRIwrite(mgz, fullfile(resultsdir, sprintf('lh.%s.mgz','vexpl_mask')));
mgz.vol = double(myresults.R2(rightidx,:))>10;
MRIwrite(mgz, fullfile(resultsdir, sprintf('rh.%s.mgz','vexpl_mask')));


pairs =[[1 2];[3 4];[5 6]]


C = [1 -1]'

betas = myresults.modelmd{2};
for p = 1 : size(pairs,1)
    
motion = nanmean(betas(:,[pairs(p,1)]),2);
stationary = nanmean(betas(:,[pairs(p,2)]),2);    
contrast = (C' * [motion stationary]')';

mgz.vol = double(contrast(leftidx,:));
MRIwrite(mgz, fullfile(resultsdir, sprintf('lh.%s_vs_%s.mgz',conditions{pairs(p,1)},conditions{pairs(p,2)})))
mgz.vol = double(contrast(rightidx,:));
MRIwrite(mgz, fullfile(resultsdir, sprintf('rh.%s_vs_%s.mgz',conditions{pairs(p,1)},conditions{pairs(p,2)})))
end

%%
%%
close all
figure(1);clf
bins = 0:1:100;
datatoplot = myresults.R2 ;
cmap0 = cmaplookup(bins,min(bins),max(bins),[],hot);

datatoplot(datatoplot==0) = 0;
datatoplot(isnan(datatoplot)) = 0;

[rawimg,Lookup,rgbimg] = cvnlookup(subj,1,datatoplot,[min(bins) max(bins)],cmap0,0,[],0,{'roiname',{'MT_exvivo'},'roicolor',{'w'},'drawrpoinames',0,'roiwidth',{5},'fontsize',20});

color = [0.5];
[r,c,t] = size(rgbimg);
[i j] = find(all(rgbimg == repmat(color,r,c,3),3));

for ii = 1: length(i)
    rgbimg(i(ii),j(ii),:) = ones(1,3);
end

a = imagesc(rgbimg); axis image tight;
axis off
hold on
% subplot(2,1,2)
plot(0,0);
colormap(cmap0);
hcb=colorbar('SouthOutside');
hcb.Ticks = [0:0.25:1];
% hcb.TickLabels = {'}
hcb.FontSize = 25
hcb.Label.String = 'R2%'
hcb.TickLength = 0.001;

title(subj)

%%
figure(2); clf
betas = myresults.modelmd{2};

motion = nanmean(betas(:,[3]),2);
stationary = nanmean(betas(:,[4]),2);

C = [1 -1]'
contrast = C' * [motion stationary]';

alltcs = nanmean(cat(3,datafiles{:}),3);
predttcs = dms{1} * myresults.modelmd{2}';


% pvals = 1-tcdf(tmp.vol(:),length(subjects)-1);
mymask = double(myresults.R2>5);



datatoplot = contrast' .* double(mymask);

datatoplot(datatoplot==0) = -50;
datatoplot(isnan(datatoplot)) = -50;



bins = -0.5:0.01:0.5
cmap0 = cmaplookup(bins,min(bins),max(bins),[],(cmapsign4));
%[rawimg,Lookup,rgbimg] = cvnlookup(subj,1,datatoplot,[min(bins) max(bins)],cmap0,min(bins),[],0,{'roiname',{'MT_exvivo'},'roicolor',{'w'},'drawrpoinames',0,'roiwidth',{5},'fontsize',20});
[rawimg,Lookup,rgbimg] = cvnlookup(subj,1,datatoplot,[min(bins) max(bins)],...
    cmap0,min(bins),[],0,{'roiname',{'Kastner*';' '},'roicolor',{'w';'k'},...
    'drawroinames',0,'roiwidth',{2,2},'fontsize',20,'surfshading',false});



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
hcb.Label.String = 'Moving dots vs. stationary'
hcb.TickLength = 0.001;





%%
figure(3);clf
mymask = double(myresults.R2>70);
sum(mymask)

tcs = nanmean(cat(3,datafiles{:}),3);
ObsResp = nanmean(tcs(logical(mymask),:),1);
dc = nanmean(ObsResp)
ObsResp = 100 * (ObsResp - dc) / dc;

plot(ObsResp)
% tcs
hold on

stem(dms{1}(:,1))

legend({'Average timecourse from 100 most responsible voxels in MT';'Centrally moving dots predictior'})
legend box off
ylabel('%BOLD')
xlabel('TRs')
set(gca,'Fontsize',25)