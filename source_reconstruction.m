%% Source reconstruction

%% Constants and options
path_ft = uigetdir([],'Give me the Field Trip folder!');
[mri_data,mri_path] = uigetfile('','Feed the MRI data [enable all data]');
[loc_file,loc_path] = uigetfile('','Feed me the EEG electrode locations [enable all data]');

path_mri = [mri_path,mri_data];
path_loc_spec = [loc_path,loc_file];

%MRI alignment
met_align = 'interactive'; % method

%MRI segmentation
tissue_seg = {'scalp', 'skull', 'csf', 'gray','white'};
vertices_mesh = [50000 50000 50000 50000 50000];
compress_mesh = [1800 1800 2200 2600 2200];

%Head model preparation
met_head = 'openmeeg'; % You can also specify 'openmeeg', 'simbio', 'bemcp', or another method.
conductivity = [0.4348 0.0217 1.5385 0.3333 0.1429]; % skull: Oostendorp 2000 in Hoekema 2003; rest: Haueisen, 1997 in Liu et al. 2018

%Sensor locations fieldtrip
path_loc_std = [path_ft '\template\electrode\standard_1005.elc'];

%Tissue probability mask of spm12
tpm_file = {'C:\Users\benjat\Documents\Matlab_Scripts\spm12\spm12\tpm\TPM.nii';'C:\Users\benjat\Documents\Matlab_Scripts\spm12\spm12\tpm\TPM.nii';'C:\Users\benjat\Documents\Matlab_Scripts\spm12\spm12\tpm\TPM.nii';'C:\Users\benjat\Documents\Matlab_Scripts\spm12\spm12\tpm\TPM.nii';'C:\Users\benjat\Documents\Matlab_Scripts\spm12\spm12\tpm\TPM.nii'};

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

%% Add path to field trip and set general settings
restoredefaultpath
addpath(path_ft)
ft_defaults
addpath([path_ft '\external\iso2mesh']);
%% Read MRI file 
mri = ft_read_mri(path_mri);

cfg = [];
cfg.resolution = 1;
cfg.dim = [256 256 256];
mri = ft_volumereslice(cfg,mri);

logFile.dim = cfg.dim;

cfg = [];
cfg.method = met_align;
cfg.coordsys = 'ctf';
[mri] = ft_volumerealign(cfg,mri);
logFile.coordsys = cfg.coordsys;


%% Segment MRI into 5 different tissue types

cfg = [];
cfg.spmversion     = 'spm12';
cfg.spmmethod      = 'mars'; %generate 6 tissues
cfg.write          = 'no';
cfg.output = tissue_seg;
cfg.tpm = tpm_file;
segmentedmri = ft_volumesegment(cfg, mri);

logFile.spmversion = cfg.spmversion;
logFile.segmentmethod = cfg.spmmethod;
logFile.segments = cfg.output;


% segment only skull because this runs better alone!
cfg = [];
cfg.spmversion     = 'spm12';
cfg.spmmethod      = 'mars'; %generate 6 tissues
cfg.write          = 'no';
cfg.output = {'skull'};
cfg.tpm = tpm_file;
segmentedmriS = ft_volumesegment(cfg, mri);

%% plot segementation
seg_i = ft_datatype_segmentation(segmentedmri,'segmentationstyle','indexed');

cfg              = [];
cfg.funparameter = 'seg';
cfg.location     = 'center';
cfg.funcolormap  = lines(5); % distinct color per tissue
cfg.atlas        = seg_i;    % the segmentation can also be used as atlas 
ft_sourceplot(cfg, seg_i);

% other plot of white and gray matter
segmentedmri.transform = mri.transform;
segmentedmri.anatomy   = mri.anatomy;

cfg              = [];
cfg.funparameter = 'gray';
ft_sourceplot(cfg,segmentedmri);


%% Prepare mesh

% add the seperate segmented skull to the other tissues
segmentedmri.skull = segmentedmriS.skull;

% for better mesh computation, fill out each tissue and make a smooth
% surface
segmentedmri.scalp = imfill(segmentedmri.scalp,'holes');
segmentedmri.scalp = imdilate(segmentedmri.scalp,strel_bol(3));
segmentedmri.skull = imfill(segmentedmri.skull,'holes');
segmentedmri.skull = imdilate(segmentedmri.skull,strel_bol(3));
segmentedmri.csf = imfill(segmentedmri.csf,'holes');
segmentedmri.csf = imdilate(segmentedmri.csf,strel_bol(3));
segmentedmri.gray = imfill(segmentedmri.gray,'holes');
segmentedmri.gray = imdilate(segmentedmri.gray,strel_bol(3));
segmentedmri.white = imfill(segmentedmri.white,'holes');
segmentedmri.white = imdilate(segmentedmri.white,strel_bol(3));

% prepare mesh for 5 tissues with high number of vertices (better results)
cfg = [];
cfg.tissue = tissue_seg;
cfg.numvertices = vertices_mesh;
cfg.spmversion   = 'spm12';
% cfg.method = 'iso2mesh';
bnd=ft_prepare_mesh(cfg,segmentedmri); %other methods 'iso2mesh' and 'isosurface' give better results but do not work with openmeeg
logFile.numvertices = cfg.numvertices;

% compress number of vertices (otherwise matlab runs out of memory)
[bnd(1).pos, bnd(1).tri] = meshresample(bnd(1).pos, bnd(1).tri,compress_mesh(1)/size(bnd(1).pos,1));
[bnd(2).pos, bnd(2).tri] = meshresample(bnd(2).pos, bnd(2).tri,compress_mesh(2)/size(bnd(2).pos,1));
[bnd(3).pos, bnd(3).tri] = meshresample(bnd(3).pos, bnd(3).tri,compress_mesh(3)/size(bnd(3).pos,1));
[bnd(4).pos, bnd(4).tri] = meshresample(bnd(4).pos, bnd(4).tri,compress_mesh(4)/size(bnd(4).pos,1));
[bnd(5).pos, bnd(5).tri] = meshresample(bnd(5).pos, bnd(5).tri,compress_mesh(5)/size(bnd(5).pos,1));

logFile.compressvertices = compress_mesh;

%% decouple intersected surfaces an repair mesh if necessery (if data quality is bad, this step has to be done several times)
bnd = decouplesurf(bnd);
 
[bnd(1).pos, bnd(1).tri] = meshcheckrepair(bnd(1).pos, bnd(1).tri,'dup');
[bnd(1).pos, bnd(1).tri] = meshcheckrepair(bnd(1).pos, bnd(1).tri,'isolated');
[bnd(1).pos, bnd(1).tri] = meshcheckrepair(bnd(1).pos, bnd(1).tri,'deep');
[bnd(1).pos, bnd(1).tri] = meshcheckrepair(bnd(1).pos, bnd(1).tri,'meshfix');

[bnd(2).pos, bnd(2).tri] = meshcheckrepair(bnd(2).pos, bnd(2).tri,'dup');
[bnd(2).pos, bnd(2).tri] = meshcheckrepair(bnd(2).pos, bnd(2).tri,'isolated');
[bnd(2).pos, bnd(2).tri] = meshcheckrepair(bnd(2).pos, bnd(2).tri,'deep');
[bnd(2).pos, bnd(2).tri] = meshcheckrepair(bnd(2).pos, bnd(2).tri,'meshfix');

[bnd(3).pos, bnd(3).tri] = meshcheckrepair(bnd(3).pos, bnd(3).tri,'dup');
[bnd(3).pos, bnd(3).tri] = meshcheckrepair(bnd(3).pos, bnd(3).tri,'isolated');
[bnd(3).pos, bnd(3).tri] = meshcheckrepair(bnd(3).pos, bnd(3).tri,'deep');
[bnd(3).pos, bnd(3).tri] = meshcheckrepair(bnd(3).pos, bnd(3).tri,'meshfix');

[bnd(4).pos, bnd(4).tri] = meshcheckrepair(bnd(4).pos, bnd(4).tri,'dup');
[bnd(4).pos, bnd(4).tri] = meshcheckrepair(bnd(4).pos, bnd(4).tri,'isolated');
[bnd(4).pos, bnd(4).tri] = meshcheckrepair(bnd(4).pos, bnd(4).tri,'deep');
[bnd(4).pos, bnd(4).tri] = meshcheckrepair(bnd(4).pos, bnd(4).tri,'meshfix');

[bnd(5).pos, bnd(5).tri] = meshcheckrepair(bnd(5).pos, bnd(5).tri,'dup');
[bnd(5).pos, bnd(5).tri] = meshcheckrepair(bnd(5).pos, bnd(5).tri,'isolated');
[bnd(5).pos, bnd(5).tri] = meshcheckrepair(bnd(5).pos, bnd(5).tri,'deep');
[bnd(5).pos, bnd(5).tri] = meshcheckrepair(bnd(5).pos, bnd(5).tri,'meshfix');

bnd = decouplesurf(bnd);

[bnd(1).pos, bnd(1).tri] = meshcheckrepair(bnd(1).pos, bnd(1).tri,'dup');
[bnd(1).pos, bnd(1).tri] = meshcheckrepair(bnd(1).pos, bnd(1).tri,'isolated');
[bnd(1).pos, bnd(1).tri] = meshcheckrepair(bnd(1).pos, bnd(1).tri,'deep');
[bnd(1).pos, bnd(1).tri] = meshcheckrepair(bnd(1).pos, bnd(1).tri,'meshfix');

[bnd(2).pos, bnd(2).tri] = meshcheckrepair(bnd(2).pos, bnd(2).tri,'dup');
[bnd(2).pos, bnd(2).tri] = meshcheckrepair(bnd(2).pos, bnd(2).tri,'isolated');
[bnd(2).pos, bnd(2).tri] = meshcheckrepair(bnd(2).pos, bnd(2).tri,'deep');
[bnd(2).pos, bnd(2).tri] = meshcheckrepair(bnd(2).pos, bnd(2).tri,'meshfix');

[bnd(3).pos, bnd(3).tri] = meshcheckrepair(bnd(3).pos, bnd(3).tri,'dup');
[bnd(3).pos, bnd(3).tri] = meshcheckrepair(bnd(3).pos, bnd(3).tri,'isolated');
[bnd(3).pos, bnd(3).tri] = meshcheckrepair(bnd(3).pos, bnd(3).tri,'deep');
[bnd(3).pos, bnd(3).tri] = meshcheckrepair(bnd(3).pos, bnd(3).tri,'meshfix');

[bnd(4).pos, bnd(4).tri] = meshcheckrepair(bnd(4).pos, bnd(4).tri,'dup');
[bnd(4).pos, bnd(4).tri] = meshcheckrepair(bnd(4).pos, bnd(4).tri,'isolated');
[bnd(4).pos, bnd(4).tri] = meshcheckrepair(bnd(4).pos, bnd(4).tri,'deep');
[bnd(4).pos, bnd(4).tri] = meshcheckrepair(bnd(4).pos, bnd(4).tri,'meshfix');

[bnd(5).pos, bnd(5).tri] = meshcheckrepair(bnd(5).pos, bnd(5).tri,'dup');
[bnd(5).pos, bnd(5).tri] = meshcheckrepair(bnd(5).pos, bnd(5).tri,'isolated');
[bnd(5).pos, bnd(5).tri] = meshcheckrepair(bnd(5).pos, bnd(5).tri,'deep');
[bnd(5).pos, bnd(5).tri] = meshcheckrepair(bnd(5).pos, bnd(5).tri,'meshfix');

bnd = decouplesurf(bnd);

% check mesh quality
figure;
ft_plot_mesh(bnd(3), 'facecolor',[0.4 0.2 0.4], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
hold on;
ft_plot_mesh(bnd(1),'edgecolor','none','facecolor',[0.4 0.3 0.4], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
hold on;
ft_plot_mesh(bnd(2),'edgecolor','none','facecolor',[0.6 0.4 0.2], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
hold on;
ft_plot_mesh(bnd(4),'edgecolor','none','facecolor',[0.2 0.8 0.8], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
hold on;
ft_plot_mesh(bnd(5),'edgecolor','none','facecolor',[0.8 0.2 0.8], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);

%% create headmodel with conductivity using openmeeg
% if this runs into an error, this is most likely the case because surfaces
% are still intersecting. Retry previous section.
save_folder = uigetdir([],'where do you want to save headmodel and leadfield?');

cfg = [];
cfg.method = met_head;
cfg.tissue = tissue_seg;
cfg.conductivity = conductivity;
vol = ft_prepare_headmodel(cfg, bnd);

logFile.conductivity = cfg.conductivity;

disp('!save headmodel. This will take a while!');

save([save_folder '\headmodel.mat'],'vol','-v7.3');





%% convert gray matter into MNE
bnd_MNE = ft_transform_geometry(segmentedmri.transform, bnd(4)); %this is used later for Leadfield


%% Read sensor locations (Currently somewhat of a hack)

elec = ft_read_sens(path_loc_std); %Read layout file

idx_keep = ismember(upper(elec.label),upper(channels)); %Find indices of channels to keep
%Updating fields

elec.chanpos = elec.chanpos(idx_keep,:);
elec.chantype = elec.chantype(idx_keep);
elec.chanunit = elec.chanunit(idx_keep);
elec.elecpos = elec.elecpos(idx_keep,:);
elec.label = elec.label(idx_keep);

logFile.EEGchanlabels = elec.label;

%% Plot electrode positions before realignment 

figure;
% head surface (scalp)
ft_plot_mesh(vol.bnd(3), 'edgecolor','none','facealpha',0.8,'facecolor',[0.6 0.6 0.8]); 
hold on;
% electrodes
ft_plot_sens(elec);
title('Electrode positions before realignment');

%% Align electrode positions to head model

%Find landmark-positions(nas, lpa, and rpa) in mri.hdr.fiducial.mri 
%(WARNING: Might have to use mri.cfg.fiducial for other files!!!)
nas = mri.cfg.fiducial.nas;
lpa = mri.cfg.fiducial.lpa;
rpa = mri.cfg.fiducial.rpa;

%Get fiducial positions in the ctf coordinate system using the
%transformation matrix of the mri and the ft_warp_apply function
transm = mri.transform;
nas = ft_warp_apply(transm,nas,'homogenous');
lpa = ft_warp_apply(transm,lpa,'homogenous');
rpa = ft_warp_apply(transm,rpa,'homogenous');

%Create a structure similar to a template set of electrodes
fid.elecpos = [nas; lpa; rpa]; %New coordinates of landmarks
fid.label = {'Nz','LPA','RPA'}; %Same labels as in elec 
fid.unit = 'mm'; %Same units as MRI
 
%Align electrodes based on landmarks
cfg = [];
cfg.method = 'fiducial';            
cfg.target = fid; %See above
cfg.elec = elec;
cfg.fiducial = {'Nz','LPA','RPA'}; %Labels of landmarks in fid and elec
elec_aligned = ft_electroderealign(cfg);

%% Plot electrode positions after realignment 

figure;
% head surface (scalp)
ft_plot_mesh(vol.bnd(1), 'edgecolor','none','facealpha',0.8,'facecolor',[0.6 0.6 0.8]); 
hold on;
% electrodes
ft_plot_sens(elec_aligned,'label','label');
title('Electrode positions after realignment');

%% Adjust final position of electrodes by eye

cfg           = [];
cfg.method    = 'interactive';
cfg.elec      = elec_aligned;
cfg.headshape = vol.bnd(1);
elec_aligned  = ft_electroderealign(cfg);

%% Prepare leadfield matrix

cfg = [];
cfg.elec = elec_aligned;
cfg.vol = vol;
cfg.grid.resolution = resolution;  
cfg.grid.unit = 'mm'; %same unit as above
cfg.channel = {'EEG','-LPA', '-RPA', '-Nz'};
% cfg.normalize = 'yes'; %only if data should not be compared later against a baseline or other conditions!
[grid,cfg] = ft_prepare_leadfield(cfg);

logFile.leadfieldResolution = resolution;
logFile.uni = cfg.grid.unit;

save([save_folder '\leadfield.mat'],'grid','-v7.3');
save([save_folder '\logFile.mat'],'logFile','-v7.3');
save([save_folder '\LFcfg.mat'],'cfg','-v7.3');



