% This script loads the already processed MIDA and transforms it into
% subject space.

% load MIDA and MRI and define tissues
path_ft = uigetdir([],'Give me the Field Trip folder!');
[mri_data,mri_path] = uigetfile('*.*','Feed the subjects MRI data');
MIDA_folder = uigetdir([],'Provide folder with processed MIDA and TPMs?');
save_folder = uigetdir([],'where do you want to save the subject space MIDA?');

restoredefaultpath
addpath(path_ft)
ft_defaults

% load MIDA and subject mri
MIDA = ft_read_mri([MIDA_folder '\MIDA_processed.nii']);


%% read a subject MRI and MIDA
 
mriSubjectBase = ft_read_mri([mri_path mri_data]); 
cfg = [];
cfg.spmversion = 'spm12';
cfg.method = 'interactive';
cfg.coordsys = 'acpc';
cfg.unit = 'mm';
[mriSubjectBase] = ft_volumerealign(cfg,mriSubjectBase);
 
% plot
cfg=[]; ft_sourceplot(cfg, mriSubjectBase);
ft_determine_coordsys(mriSubjectBase, 'interactive', 'no');
 
%% Coregister and reslice the subject to the MIDA

mriSubject = mriSubjectBase;
 
cfg = [];
cfg.method = 'flip';
[mriSubjectResliced] = ft_volumereslice(cfg, mriSubject);

%plot
cfg=[]; ft_sourceplot(cfg, mriSubjectResliced);
 
cfg = [];
cfg.method = 'nearest';
cfg.resolution = 1;
cfg.dim=round(mriSubjectResliced.dim.*mriSubjectResliced.hdr.volres./(MIDA.hdr.volres*3));
cfg.coordsys = 'acpc';
cfg.unit = 'mm';
[mriSubjectResliced] = ft_volumereslice(cfg, mriSubjectResliced);
 
% plot
cfg=[]; ft_sourceplot(cfg, mriSubjectResliced);
 

cfg = [];
cfg.spmversion = 'spm12';
cfg.method = 'spm';
cfg.coordsys = 'acpc';
cfg.unit = 'mm';
cfg.parameter = 'anatomy';
cfg.viewresult = 'yes';
[mriSubjectRealigned] = ft_volumerealign(cfg,mriSubjectResliced,MIDA);
 
% plot
cfg=[]; ft_sourceplot(cfg, mriSubjectRealigned);
 
% write subjects MRI. This is needed for interpolation later
save([save_folder '\MRI_processed.mat'],'mriSubjectRealigned','-v7.3');
 
%% Segment the subject using the MIDA TPMs
cfg = [];
cfg.spmversion = 'spm12';
cfg.spmmethod = 'new';
cfg.output = {'tpm'};
cfg.tpm = {[MIDA_folder '\TPM_gray.nii'], [MIDA_folder '\TPM_white.nii'], [MIDA_folder '\TPM_CSF.nii'], [MIDA_folder '\TPM_bone.nii'], [MIDA_folder '\TPM_soft.nii'], [MIDA_folder '\TPM_background.nii']};
opts.lkp = [1,1,2,2,3,3,4,4,4,4,5,5,6,6];
cfg.opts=opts;
cfg.write = 'yes';
cfg.name = 'segmentedSubject';
 
segmentedSubject = ft_volumesegment(cfg, mriSubjectRealigned);
 
%% Apply transformation to the MIDA to get it into subject space

[mri] = ft_volumenormalise(cfg, MIDA);

%% Start making the mesh and so forth