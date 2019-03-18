subjID = 's12';

path_ft = uigetdir([],'Give me the Field Trip folder!');
path_spm = uigetdir([],'Give me the SPM folder!');
save_folder = uigetdir([],'where do you want to save the MIDA template and TPM?');

if ispc
    save_folder = [save_folder '\'];
else
    save_folder = [save_folder '/'];
end

[MIDA, MIDApath] = uigetfile('*.nii', 'Pick unmodified MIDA template');
[subjectMRI, subjectMRIpath] = uigetfile('*.nii', 'Pick subject MRI');

[loc_file,loc_path] = uigetfile('*.*','Feed me the EEG electrode locations');

%% if subject MIR only a DICOM run this section
restoredefaultpath;
addpath(path_ft);
ft_defaults;

[subjectMRI_DICOM, subjectMRIpath_DICOM] = uigetfile('*.*', 'Pick subject MRI in DICOM');
save_folder_nifti = uigetdir([],'where do you want to save NIFTI file?');

mri = ft_read_mri([subjectMRIpath_DICOM subjectMRI_DICOM]);

cfg = [];
cfg.filename  = [save_folder_nifti '\' subjectMRI_DICOM];
cfg.filetype  = 'nifti';
cfg.parameter = 'anatomy';
ft_sourcewrite(cfg, mri);

% [subjectMRI, subjectMRIpath] = uigetfile('*.nii', 'Pick subject MRI');
subjectMRI = [subjectMRI_DICOM '.nii'];
subjectMRIpath = [save_folder_nifti '\'];

%% Process MIDA with fieldtrip
restoredefaultpath;
addpath(path_ft);
ft_defaults;

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

UiO_process_MIDA(MIDA, MIDApath, save_folder);

%% Warp to subject space with SPM
restoredefaultpath;
addpath(genpath(path_spm));

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

UiO_transform_MIDA_to_subject(save_folder, subjectMRI, subjectMRIpath, save_folder, subjID);

%% Create 12 tissue FEM model
restoredefaultpath;
addpath(path_ft);
ft_defaults;

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

UiO_forward_model_12T_FEM( [subjID 'MIDA_withoutair.nii'], save_folder, save_folder, loc_file, loc_path, path_ft);

%% Calculate inverse solution
restoredefaultpath;
addpath(path_ft);
ft_defaults;

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

[LF_Head_path] = uigetdir([],'Feed the leadfield and headmodel folder');

if ispc
    LF_Head_path = [LF_Head_path '\'];
else
    LF_Head_path = [LF_Head_path '/'];
end

%[MRI_processed] = uigetdir([],'Feed the MIDA and MRI processed folder');
% [mri_data,mri_path] = uigetfile('*.*','Feed the subjects MRI data');
[EEG_file,EEG_path] = uigetfile('*.*','Feed me the EEG data');

LF_data = [LF_Head_path,'leadfield_12T_FEM_gray-only.mat'];
Head_data = [LF_Head_path,'headmodel_12T_FEM_prepared_sens_vol.mat'];

UiO_inverse_model_eLORETA(['r' subjectMRI], subjectMRIpath, EEG_file, EEG_path, LF_data, Head_data, LF_Head_path, path_ft, save_folder);


