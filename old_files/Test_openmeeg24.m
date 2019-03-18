path_ft = uigetdir([],'Give me the Field Trip folder!');
[mri_data,mri_path] = uigetfile('*.*','Feed the MRI data');
[loc_file,loc_path] = uigetfile('*.*','Feed me the EEG electrode locations');
save_folder = uigetdir([],'where do you want to save MRI, headmodel, and leadfield?');

path_mri = [mri_path,mri_data];
path_loc_spec = [loc_path,loc_file];

%MRI alignment
met_align = 'interactive'; % method

%Head model preparation
met_head = 'openmeeg'; % You can also specify 'openmeeg', 'simbio', 'bemcp', or another method.
conductivity = [0.1429 0.3333 1.5385 0.0217 0.4348]; % order: white, gray, csf, skull, scalp. (skull: Oostendorp 2000 in Hoekema 2003; rest: Haueisen, 1997 in Liu et al. 2018)

%Sensor locations fieldtrip
path_loc_std = [path_ft '\template\electrode\standard_1005.elc'];

%Tissue probability mask of spm12
% tpm_file = {'Z:\Matlab_Scripts\spm12\spm12\tpm\TPM.nii';'Z:\Matlab_Scripts\spm12\spm12\tpm\TPM.nii';'Z:\Matlab_Scripts\spm12\spm12\tpm\TPM.nii';'Z:\Matlab_Scripts\spm12\spm12\tpm\TPM.nii';'Z:\Matlab_Scripts\spm12\spm12\tpm\TPM.nii'};

% Add path to field trip and set general settings
restoredefaultpath
addpath(path_ft)
ft_defaults
addpath([path_ft '\external\iso2mesh']);

%Channel customization (hack for Attentional load data)
path_loc_attload = path_loc_spec; %Path to file
elec_attload = ft_read_sens(path_loc_attload);
elec_attload.label{strcmpi(elec_attload.label,'AF7')} = 'AFz';
elec_attload.label{strcmpi(elec_attload.label,'AF8')} = 'FCz';
idx_eog = cellfun(@(x)contains(x,'EOG','IgnoreCase',true),elec_attload.label);

%Channels to keep
channels = ['LPA';'RPA';'Nz'; elec_attload.label(~idx_eog)];

%Options for leadfield matrix
resolution = 6; %mm

% save all parameters to logFile
logFile = [];

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

%%

mri = ft_read_mri(path_mri);

cfg            = [];
cfg.downsample = 1;
mri = ft_volumenormalise(cfg, mri);

%%

cfg = [];
cfg.output = {'brain','skull','scalp'};
segmentedmri = ft_volumesegment(cfg, mri);
% 
cfg = [];
cfg.output = {'scalp'};
segmentedmriS = ft_volumesegment(cfg, mri);

segmentedmri.scalp = segmentedmriS.scalp;
cfg = [];
cfg.tissue = {'scalp','skull','brain'};
cfg.numvertices = 2000;
bnd=ft_prepare_mesh(cfg,segmentedmri);


%%
subjId = 'test_1001';

vol = [];
vol.bnd = bnd;
vol.cond = [0.33 0.0041 0.33]; 
vol.type = 'openmeeg';
% vol.basefile = subjId;
% vol.path = ['./' subjId '/hm/openmeeg_out'];

% load Z:\Matlab_Scripts\fieldtrip\fieldtrip-20180610\template\sourcemodel\standard_sourcemodel3d5mm.mat


cfg = [];
cfg.reducerank = 2;
cfg.vol = ft_convert_units(vol,'mm');
cfg.elec = elec_aligned;
% cfg.grid.pos = ft_warp_apply(mri, sourcemodel.pos, 'sn2individual');
grid = ft_prepare_leadfield(cfg);