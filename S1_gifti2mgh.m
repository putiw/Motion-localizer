bidsDir = '/Users/pw1246/Desktop/MRI/CueIntegration';
bidsDir = '/Users/pw1246/mnt/CBIUserData/rokerslab/CueIntegration';
addpath(genpath(pwd));
datadir = [bidsDir '/derivatives/fmriprep'];
subjects = dir(sprintf('%s/sub*',datadir));
subjects = subjects([subjects.isdir]);
   
setenv('PATH', ['/Applications/freesurfer/7.2.0/bin' getenv('PATH')]);
setenv('SUBJECTS_DIR','/Applications/freesurfer/7.2.0/subjects/')
setenv('FREESURFER_HOME','/Applications/freesurfer/7.2.0')
setenv('FSLOUTPUTTYPE','NIFTI_GZ')
addpath '/Users/pw1246/Documents/GitHub/wpToolbox'
tempPath = '/Applications/freesurfer/7.2.0/bin/';
surf = {'fsnative'}
ses = 'ses-01';

surf = surf(1);
for s = 1% : length(subjects)
    
    subj_dir = sprintf('%s/%s',datadir,subjects(s).name)
    
    for srf = 1 : length(surf)

        d = dir(sprintf('%s/ses-01/func/*%s*.gii',subj_dir,surf{srf}))
        %dd = dir(sprintf('%s/func/*%s*.gii',subj_dir,surf{srf}))
        freesurfer_init
        for run = 1 : length(d)
            
            mgzfile = sprintf('%s.mgh',d(run).name(1:end-4))
            
            %system(sprintf('%smri_convert %s/ses-01/func/%s  %s/ses-nyu3t01/func/%s',freesurfer_string,subj_dir,d(run).name,subj_dir,mgzfile))
            system(sprintf(['%s' tempPath 'mri_convert %s/%s/func/%s  %s/%s/func/%s'],freesurfer_string,subj_dir,ses,d(run).name,subj_dir,ses,mgzfile))
            
        end
        
    end
end

% file1 = '/Users/pw1246/mnt/CBIUserData/rokerslab/CueIntegration_nosession/derivatives/fmriprep/sub-0201/func/sub-0201_task-mt_run-2_space-fsnative_hemi-L_bold.func.mgz';
% file2 = '/Users/pw1246/mnt/CBIUserData/rokerslab/CueIntegration_nosession/derivatives/fmriprep/sub-0201/func/sub-0201_task-mt_run-2_space-fsnative_hemi-L_bold.func.mgh';
% system(sprintf(['%s' tempPath 'mri_convert %s  %s'],freesurfer_string,file1,file2))

% 
% for s = 1% : length(subjects)
%     
%     subj_dir = sprintf('%s/%s',datadir,subjects(s).name)
%     
%     for srf = 1 : length(surf)
% 
%         d = dir(sprintf('%s/func/*%s*.gii',subj_dir,surf{srf}))
%         freesurfer_init
%         for run = 1 : length(d)
%             
%             
%             mgzfile = sprintf('%s.mgz',d(run).name(1:end-4))
%             
%             %system(sprintf('%smri_convert %s/ses-01/func/%s  %s/ses-nyu3t01/func/%s',freesurfer_string,subj_dir,d(run).name,subj_dir,mgzfile))
%             system(sprintf(['%s' tempPath 'mri_convert %s/func/%s  %s/func/%s'],freesurfer_string,subj_dir,d(run).name,subj_dir,mgzfile))
%             
%         end
%         
%     end
% end
% % 
% % surf = surf(1);
% % for s = 1% : length(subjects)
% %     
% %     subj_dir = sprintf('%s/%s',datadir,subjects(s).name)
% %     
% %     for srf = 1 : length(surf)
% % 
% %         d = dir(sprintf('%s/ses-01/func/*%s*.gii',subj_dir,surf{srf}))
% %         freesurfer_init
% %         for run = 1 : length(d)
% %             
% %             
% %             mgzfile = sprintf('%s.mgz',d(run).name(1:end-4))
% %             
% %             %system(sprintf('%smri_convert %s/ses-01/func/%s  %s/ses-nyu3t01/func/%s',freesurfer_string,subj_dir,d(run).name,subj_dir,mgzfile))
% %             system(sprintf(['%s' tempPath 'mri_convert %s/ses-01/func/%s  %s/ses-01/func/%s'],freesurfer_string,subj_dir,d(run).name,subj_dir,mgzfile))
% %             
% %         end
% %         
% %     end
% % end