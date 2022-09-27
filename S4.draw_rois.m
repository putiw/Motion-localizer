%setenv('SUBJECTS_DIR','/Volumes/server/Projects/akinetopsia/derivatives/freesurfer/');
addpath(genpath('./../nsdcode'))
roilabels = {'hMT';'MST';}
rng    = [0 2]
subjid = 'sub-0201';
ses = 'ses-01';
mgznames = {'left_moving_vs_left_stationary';'right_moving_vs_right_stationary';}; %{'angle_adj';'eccen';}
for zz = 1:3

cmap = cmapsign4;
bins = [-0.5:0.1:0.5];
cmaps{zz} = cmaplookup(bins,min(bins),max(bins),[],cmap);
crngs{zz} = [min(bins) max(bins)];
threshs{zz} = min(bins);
end
roivals = [];
% 
bidsDir = '/Users/pw1246/Desktop/MRI/CueIntegration';
conditions = {'central_moving';'central_stationary';'left_moving';'left_stationary';'right_moving';'right_stationary'}
resultsdir = sprintf('%s/derivatives/GLMdenoise/%s/%s/',bidsDir,subjid,ses);
githubDir = '/Users/pw1246/Documents/GitHub';
codeDir = pwd;
user = 'puti';
projectName = 'Localizer';
addpath(genpath(fullfile(githubDir, 'wpToolbox')));
setup_user(user,projectName,bidsDir,githubDir);

p = 2;
% thr_map_lh = MRIread(fullfile(resultsdir, sprintf('lh.%s_vs_%s.mgz',conditions{pairs(p,1)},conditions{pairs(p,2)})));
% thr_map_rh = MRIread(fullfile(resultsdir, sprintf('rh.%s_vs_%s.mgz',conditions{pairs(p,1)},conditions{pairs(p,2)})));
% thr_map_lh = MRIread('/Volumes/server/Projects/akinetopsia/derivatives/prfvista/sub-wlsubj140/ses-nyu3t01/lh.vexpl_prf.mgz')
% thr_map_rh = MRIread('/Volumes/server/Projects/akinetopsia/derivatives/prfvista/sub-wlsubj140/ses-nyu3t01/rh.vexpl_prf.mgz')

%mythr = double([thr_map_lh.vol thr_map_rh.vol] > 0.01);
%cvndefinerois
drawMyRois
%%
fileIn = '/Users/pw1246/Documents/GitHub/Cue_decoding_pipeline/lh.MT.mgz';
fileOut = '/Users/pw1246/Documents/GitHub/Cue_decoding_pipeline/lh.MT.nii';

%!mri_surf2vol(--surfval --hemi lh fileIn --fillribbon --o fileOut)
% mri_surf2vol --surfval /Users/pw1246/Documents/GitHub/Cue_decoding_pipeline/lh.MT.mgz --hemi lh --fillribbon --o /Users/pw1246/Documents/GitHub/Cue_decoding_pipeline/lh.MT.nii --identity sub-0201 --template /Users/pw1246/Desktop/MRI/CueIntegration/derivatives/freesurfer/sub-0201/mri/brain.mgz
%freeview /Users/pw1246/Desktop/MRI/CueIntegration/derivatives/freesurfer/sub-0201/mri/T1.mgz  /Users/pw1246/Documents/GitHub/Cue_decoding_pipeline/importantMGZs/lh.MT.nii
%%
close all
figure(1);clf
datatoplot = roivals;

cmap = hsv
bins = 1 : 6;
cmap0 = cmaplookup(bins,min(bins),max(bins),[],cmap);
[rawimg,Lookup,rgbimg] = cvnlookup(subjid,1,datatoplot,[min(bins) max(bins)],cmap0,0,[],0,{'roiname',{'Kastner2015';' '},'roicolor',{'w';'k'},'drawroinames',0,'roiwidth',{2,2},'fontsize',20,'surfshading',false});
imagesc(rgbimg)
axis image off


%%
figure(2);clf
hemi = {'lh';'rh'}

filename = 'left_moving_vs_left_stationary'

datatoplot = []
for h = 1 : length(hemi)
    
    tmp = MRIread(sprintf('%s/%s.%s.mgz',resultsdir,hemi{h},filename))
    datatoplot = [datatoplot;tmp.vol];
    
end

cmap = cmapsign4
bins = [-0.5:0.1:0.5];
cmap0 = cmaplookup(bins,min(bins),max(bins),[],cmap);
% [rawimg,Lookup,rgbimg] = cvnlookup(subjid,1,datatoplot,[min(bins) max(bins)],cmap0,min(bins),[],0,{'roiname',{'Kastner2015';' '},'roicolor',{'w';'k'},'drawroinames',0,'roiwidth',{2,2},'fontsize',20,'surfshading',false});
[rawimg,Lookup,rgbimg] = cvnlookup(subjid,{ {{[-90 30 0] [90 30 0]} {'lh' 'rh'}} 'inflated' 0 1500 0 []},datatoplot,[min(bins) max(bins)],cmap0,min(bins),[],0,{'roiname',{'Kastner2015';' '},'roicolor',{'w';'k'},'drawroinames',0,'roiwidth',{2,2},'fontsize',20,'surfshading',false});



imagesc(rgbimg)
axis image off
% colormap(cmap0(1:length(eccbins_lower),:));
colormap(cmap)
c=colorbar;
c.Ticks=[0 1];
c.TickLabels = {num2str(min(bins))},{num2str(min(bins))};
c.Label.String = 'Contrast'
set(gca,'FontSize',15)



%%  Central
figure(3);clf
hemi = {'lh';'rh'}
filename = 'central_moving_vs_central_stationary'
%filename = 'MT'

datatoplot = []
for h = 1 : length(hemi)
        tmp = MRIread(sprintf('%s/%s.%s.mgz',resultsdir,hemi{h},filename))
    datatoplot = [datatoplot;reshape(tmp.vol,numel(tmp.vol),1)];
    
end

thr = []
varexpl = 'vexpl_glm'
for h = 1 : length(hemi)
        tmp = MRIread(sprintf('%s/%s.%s.mgz',resultsdir,hemi{h},varexpl))
    thr = [thr;tmp.vol];
    
end


datatoplot = datatoplot .* double(thr > 0);
datatoplot(datatoplot == 0) = -50;

cmap = cmapsign4
bins = [-0.5:0.1:0.5];

%bins = [0:0.1:2];
cmap0 = cmaplookup(bins,min(bins),max(bins),[],cmap);
% [rawimg,Lookup,rgbimg] = cvnlookup(subjid,1,datatoplot,[min(bins) max(bins)],cmap0,min(bins),[],0,{'roiname',{'Kastner2015';' '},'roicolor',{'w';'k'},'drawroinames',0,'roiwidth',{2,2},'fontsize',20,'surfshading',false});
[rawimg,Lookup,rgbimg] = cvnlookup(subjid,1,datatoplot,[min(bins) max(bins)],cmap0,min(bins),[],0,{'roiname',{'Kastner2015';' '},'roicolor',{'w';'k'},'drawroinames',0,'roiwidth',{3,10},'fontsize',20,'surfshading',false});



imagesc(rgbimg)
axis image off
% colormap(cmap0(1:length(eccbins_lower),:));
colormap(cmap)
c=colorbar;
c.Ticks=[0 1];
c.TickLabels = {num2str(min(bins))},{num2str(min(bins))};
c.Label.String = 'Contrast'
set(gca,'FontSize',15)
%% Left
figure(4);clf
hemi = {'lh';'rh'}
filename = 'left_moving_vs_left_stationary'

datatoplot = []
for h = 1 : length(hemi)
        tmp = MRIread(sprintf('%s/%s.%s.mgz',resultsdir,hemi{h},filename))
    datatoplot = [datatoplot;tmp.vol];
    
end

thr = []
varexpl = 'vexpl_glm'
for h = 1 : length(hemi)
        tmp = MRIread(sprintf('%s/%s.%s.mgz',resultsdir,hemi{h},varexpl))
    thr = [thr;tmp.vol];
    
end


datatoplot = datatoplot .* double(thr > 20);
datatoplot(datatoplot == 0) = -50;

cmap = cmapsign4
bins = [-0.5:0.1:0.5];
cmap0 = cmaplookup(bins,min(bins),max(bins),[],cmap);
% [rawimg,Lookup,rgbimg] = cvnlookup(subjid,1,datatoplot,[min(bins) max(bins)],cmap0,min(bins),[],0,{'roiname',{'Kastner2015';' '},'roicolor',{'w';'k'},'drawroinames',0,'roiwidth',{2,2},'fontsize',20,'surfshading',false});
[rawimg,Lookup,rgbimg] = cvnlookup(subjid,1,datatoplot,[min(bins) max(bins)],cmap0,min(bins),[],0,{'roiname',{'Kastner2015';' '},'roicolor',{'w';'k'},'drawroinames',0,'roiwidth',{2,2},'fontsize',20,'surfshading',false});



imagesc(rgbimg)
axis image off
% colormap(cmap0(1:length(eccbins_lower),:));
colormap(cmap)
c=colorbar;
c.Ticks=[0 1];
c.TickLabels = {num2str(min(bins))},{num2str(min(bins))};
c.Label.String = 'Contrast'
set(gca,'FontSize',15)

%% Right
figure(5);clf
hemi = {'lh';'rh'}
filename = 'right_moving_vs_right_stationary'

datatoplot = []
for h = 1 : length(hemi)
        tmp = MRIread(sprintf('%s/%s.%s.mgz',resultsdir,hemi{h},filename))
    datatoplot = [datatoplot;tmp.vol];
    
end

thr = []
varexpl = 'vexpl_glm'
for h = 1 : length(hemi)
        tmp = MRIread(sprintf('%s/%s.%s.mgz',resultsdir,hemi{h},varexpl))
    thr = [thr;tmp.vol];
    
end


datatoplot = datatoplot .* double(thr > 0);
datatoplot(datatoplot == 0) = -50;

cmap = cmapsign4
bins = [-0.5:0.1:0.5];
cmap0 = cmaplookup(bins,min(bins),max(bins),[],cmap);
% [rawimg,Lookup,rgbimg] = cvnlookup(subjid,1,datatoplot,[min(bins) max(bins)],cmap0,min(bins),[],0,{'roiname',{'Kastner2015';' '},'roicolor',{'w';'k'},'drawroinames',0,'roiwidth',{2,2},'fontsize',20,'surfshading',false});
[rawimg,Lookup,rgbimg] = cvnlookup(subjid,1,datatoplot,[min(bins) max(bins)],cmap0,min(bins),[],0,{'roiname',{'Kastner2015';' '},'roicolor',{'w';'k'},'drawroinames',0,'roiwidth',{2,2},'fontsize',20,'surfshading',false});



imagesc(rgbimg)
axis image off
% colormap(cmap0(1:length(eccbins_lower),:));
colormap(cmap)
c=colorbar;
c.Ticks=[0 1];
c.TickLabels = {num2str(min(bins))},{num2str(min(bins))};
c.Label.String = 'Contrast'
set(gca,'FontSize',15)
%% central
figure(6);clf
hemi = {'lh';'rh'}
filename = 'central_moving_vs_central_stationary'

datatoplot = []
for h = 1 : length(hemi)
        tmp = MRIread(sprintf('%s/%s.%s.mgz',resultsdir,hemi{h},filename))
    datatoplot = [datatoplot;tmp.vol];
    
end

thr = []
varexpl = 'vexpl_glm'
for h = 1 : length(hemi)
        tmp = MRIread(sprintf('%s/%s.%s.mgz',resultsdir,hemi{h},varexpl))
    thr = [thr;tmp.vol];
    
end


datatoplot = datatoplot .* double(thr > 0);
datatoplot(datatoplot == 0) = -50;

cmap = cmapsign4
bins = [-0.5:0.1:0.5];
cmap0 = cmaplookup(bins,min(bins),max(bins),[],cmap);
% [rawimg,Lookup,rgbimg] = cvnlookup(subjid,1,datatoplot,[min(bins) max(bins)],cmap0,min(bins),[],0,{'roiname',{'Kastner2015';' '},'roicolor',{'w';'k'},'drawroinames',0,'roiwidth',{2,2},'fontsize',20,'surfshading',false});
[rawimg,Lookup,rgbimg] = cvnlookup(subjid,1,datatoplot,[min(bins) max(bins)],cmap0,min(bins),[],0,{'roiname',{'Kastner2015';' '},'roicolor',{'w';'k'},'drawroinames',0,'roiwidth',{2,2},'fontsize',20,'surfshading',false});



imagesc(rgbimg)
axis image off
% colormap(cmap0(1:length(eccbins_lower),:));
colormap(cmap)
c=colorbar;
c.Ticks=[0 1];
c.TickLabels = {num2str(min(bins))},{num2str(min(bins))};
c.Label.String = 'Contrast'
set(gca,'FontSize',15)
