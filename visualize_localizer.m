clc; clear all; close all;
bidsDir = '/Users/pw1246/Desktop/MRI/CueIntegration';
githubDir = '/Users/pw1246/Documents/GitHub';
user = 'puti';
projectName = 'Localizer';
addpath(genpath(fullfile(githubDir, 'wpToolbox')));
setup_user(user,projectName,bidsDir,githubDir);
roilabels = {'hMT';'MST';}
rng    = [0 2]
subjid = 'sub-0201';
ses = 'ses-01';
resultsdir = sprintf('%s/derivatives/GLMdenoise/%s/%s/',bidsDir,subjid,ses);
mgznames = {'central_moving_vs_central_stationary';'left_moving_vs_left_stationary';'right_moving_vs_right_stationary';'MT'}; %{'angle_adj';'eccen';}
for zz = 1:3
cmap = cmapsign4;
bins = [-0.5:0.1:0.5];
cmaps{zz} = cmaplookup(bins,min(bins),max(bins),[],cmap);
crngs{zz} = [min(bins) max(bins)];
threshs{zz} = min(bins);
end
zz = 4;
cmap = jet;
bins = [0:0.1:2];
cmaps{zz} = cmaplookup(bins,min(bins),max(bins),[],cmap);
crngs{zz} = [min(bins) max(bins)];
threshs{zz} = 0.5;
roivals = [];

drawMyRois

