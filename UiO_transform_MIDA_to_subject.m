function [] = UiO_transform_MIDA_to_subject( folderMIDA, subjectMRI, folderSubject, folderWrite, subjID )
% Parameters should be strings with filename of the corresonding
% *.nii files. 

cd(folderWrite);
MIDA = 'MIDA_withoutair.nii';

spm_jobman('initcfg');

job_number = 0;

%% coregister and reslice

job_number = job_number+1;

fprintf('Creating coregistration job...\n');
matlabbatch{job_number}.spm.spatial.coreg.estwrite.ref = {[folderMIDA MIDA ',1']};
matlabbatch{job_number}.spm.spatial.coreg.estwrite.source = {[folderSubject subjectMRI ',1']};
matlabbatch{job_number}.spm.spatial.coreg.estwrite.other = {''};
matlabbatch{job_number}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{job_number}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{job_number}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{job_number}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{job_number}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{job_number}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{job_number}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{job_number}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

fprintf('Running coregistration...\n');
spm_jobman('run',matlabbatch(job_number));

%% 

job_number = job_number+1;

fprintf('Creating segmentation job...\n');
matlabbatch{job_number}.spm.spatial.preproc.channel.vols = {[folderSubject 'r' subjectMRI ',1']};
matlabbatch{job_number}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{job_number}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{job_number}.spm.spatial.preproc.channel.write = [0 0];
matlabbatch{job_number}.spm.spatial.preproc.tissue(1).tpm = {[folderWrite 'TPM_gray.nii,1']};
matlabbatch{job_number}.spm.spatial.preproc.tissue(1).ngaus = 2;
matlabbatch{job_number}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{job_number}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{job_number}.spm.spatial.preproc.tissue(2).tpm = {[folderWrite 'TPM_white.nii,1']};
matlabbatch{job_number}.spm.spatial.preproc.tissue(2).ngaus = 2;
matlabbatch{job_number}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch{job_number}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{job_number}.spm.spatial.preproc.tissue(3).tpm = {[folderWrite 'TPM_CSF.nii,1']};
matlabbatch{job_number}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{job_number}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{job_number}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{job_number}.spm.spatial.preproc.tissue(4).tpm = {[folderWrite 'TPM_bone.nii,1']};
matlabbatch{job_number}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{job_number}.spm.spatial.preproc.tissue(4).native = [1 0];
matlabbatch{job_number}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{job_number}.spm.spatial.preproc.tissue(5).tpm = {[folderWrite 'TPM_soft.nii,1']};
matlabbatch{job_number}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{job_number}.spm.spatial.preproc.tissue(5).native = [1 0];
matlabbatch{job_number}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{job_number}.spm.spatial.preproc.tissue(6).tpm = {[folderWrite 'TPM_background.nii,1']};
matlabbatch{job_number}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{job_number}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{job_number}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{job_number}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{job_number}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{job_number}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{job_number}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{job_number}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{job_number}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{job_number}.spm.spatial.preproc.warp.write = [1 1];

fprintf('Running segmentation...\n');
spm_jobman('run',matlabbatch(job_number));

%%

job_number = job_number+1;

fprintf('Creating normalization job...\n');
matlabbatch{job_number}.spm.spatial.normalise.write.subj.def = {[folderSubject 'iy_r' subjectMRI]};
matlabbatch{job_number}.spm.spatial.normalise.write.subj.resample = {[folderMIDA MIDA  ',1']};
matlabbatch{job_number}.spm.spatial.normalise.write.woptions.bb = spm_get_bbox([folderSubject 'r' subjectMRI],'nn');%nan(2,3); %[-120 -100 -125
                                                         % 120 100 125];
%matlabbatch{job_number}.spm.spatial.normalise.write.woptions.bb = [-90 -125 -150
%                                                          90 120 80];
matlabbatch{job_number}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
matlabbatch{job_number}.spm.spatial.normalise.write.woptions.interp = 0;
matlabbatch{job_number}.spm.spatial.normalise.write.woptions.prefix = subjID;

fprintf('Running normalization...\n');
spm_jobman('run',matlabbatch(job_number));

disp("Done. MIDA warped to subject space");

end

