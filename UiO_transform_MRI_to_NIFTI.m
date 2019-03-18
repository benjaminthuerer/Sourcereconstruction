%Transform DICOM to nifti

function UiO_transform_MRI_to_NIFTI()

path_ft = uigetdir([],'Give me the Field Trip folder!');
[subjectMRI_DICOM, subjectMRIpath_DICOM] = uigetfile('*.*', 'Pick subject MRI in DICOM');
save_folder_nifti = uigetdir([],'where do you want to save NIFTI file?');

restoredefaultpath;
addpath(path_ft);
ft_defaults;

mri = ft_read_mri([subjectMRIpath_DICOM subjectMRI_DICOM]);

cfg = [];
cfg.filename  = [save_folder_nifti '\' subjectMRI_DICOM];
cfg.filetype  = 'nifti';
cfg.parameter = 'anatomy';
ft_sourcewrite(cfg, mri);

end