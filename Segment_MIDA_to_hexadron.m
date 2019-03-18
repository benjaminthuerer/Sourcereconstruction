path_ft = uigetdir([],'Give me the Field Trip folder!');
[mri_data,mri_path] = uigetfile('*.*','Feed the subjects MRI data');
save_folder = uigetdir([],'where do you want to save MRI-template, headmodel, and leadfield?');

restoredefaultpath
addpath(path_ft)
ft_defaults

% load MIDA and subject mri
MIDA = ft_read_mri('Z:\07_fNetworks_rest-state\MIDA_head-model\MIDAv1.0\MIDA_v1.0\MIDA_v1_voxels\MIDA_v1.nii');
mri = ft_read_mri([mri_path mri_data]);

% Define tissue types
Tissues = {'skin','eyes','muscle','fat','spongybone','compactbone','gray','cerebellargray','white','cerebellarwhite','csf','brainstem'};

skin = [1,33:35,37,39,51,85,86];
eyes = [55:59];
muscle = [38,42,60,61,63:84,88:96,98];
fat = [43,62];
spongybone = [52];
compactbone = [36,40,41,44:49,53,54,87];
gray = [3:5,7,8,10,16,17,19,20,21,99,116];
cerebellargray = [2];
white = [12,18,22,23,100:115];
cerebellarwhite = [9];
csf = [6,24,25,32];
brainstem = [11,13:15];

% special tissue for GrayMatter
SegGrayMatter = sort([gray cerebellargray]);

%not segmented
air = [26:31,97];
bg = 50;
% hypophysis = [19];
% blood = [24,25];

bg = ismember(MIDA.anatomy,bg);
MIDA.anatomy(bg) = 0;

%% preprocess subjects mri and MIDA for subsequent warping

cfg = [];
cfg.spmversion     = 'spm12';
cfg.resolution = 1;
cfg.dim = [256 256 256];
mri = ft_volumereslice(cfg,mri);

cfg = [];
cfg.spmversion     = 'spm12';
cfg.method = 'interactive';
cfg.coordsys = 'acpc';
[mri] = ft_volumerealign(cfg,mri);




cfg = [];
cfg.spmversion     = 'spm12';
cfg.method = 'nearest';
cfg.downsample = 2;
cfg.dim = [256 256 256];
MIDA = ft_volumereslice(cfg,MIDA);

% ft_determine_coordsys(mri, 'interactive', 'no')

cfg = [];
cfg.method = 'interactive';
cfg.coordsys = 'acpc';
MIDA = ft_volumerealign(cfg, MIDA);

% ft_determine_coordsys(MIDA, 'interactive', 'no')

cfg = [];
cfg.method = 'spm';
cfg.coordsys = 'acpc';
cfg.parameter = 'anatomy';
cfg.viewresult = 'yes';
MIDA = ft_volumerealign(cfg,MIDA,mri);

save([save_folder '\MRI_MIDA_norm.mat'],'MIDA','-v7.3');
save([save_folder '\MRI_subject_norm.mat'],'mri','-v7.3');
disp('MRIs saved');

% ft_determine_coordsys(MIDA, 'interactive', 'no')

% cfg = [];
% cfg.method = 'interactive';
% cfg.coordsys = 'ctf';
% MIDA = ft_volumerealign(cfg, MIDA);

% cfg = [];
% cfg.write          = 'yes';
% cfg.spmversion     = 'spm12';
% cfg.spmmethod = 'new';
% % cfg.output = 'skullstrip';
% cfg.name = [save_folder '\TPM_template'];
% TPM = ft_volumesegment(cfg, TPM);
% 
% 
% 
% cfg = [];
% cfg.spmversion    = 'spm12';
% % cfg.coordsys = 'mni';
% cfg.resolution = 1;
% cfg.nonlinear = 'yes';
% test = ft_volumenormalise(cfg, MIDA);
% 
% cfg = [];
% cfg.filename  = [save_folder '\source_activity'];
% cfg.filetype  = 'nifti';
% cfg.parameter = 'anatomy';
% cfg.precision = 'single';
% ft_sourcewrite(cfg, TPM);
% 
% 
% copyMIDA = MIDA;
% 
% mri = copyMRI;
% 
% MIDA = copyMIDA;
% 
% 
% cfg = [];
% cfg.spmversion    = 'spm12';
% cfg.resolution = 1;
% cfg.nonlinear = 'no';
% cfg.coordsys = 'acpc';
% cfg.template = [save_folder '\TPM_template.nii'];
% MIDA = ft_volumenormalise(cfg,MIDA);
% 
% [MIDA] = ft_convert_units(MIDA, 'mm');
% 
% ft_determine_coordsys(MIDA, 'interactive', 'no')
% 
% cfg = [];
% cfg.method = 'spm';
% cfg.coordsys = 'acpc';
% cfg.parameter = 'anatomy';
% cfg.viewresult = 'yes';
% MIDA = ft_volumerealign(cfg,MIDA,mri);
% 
% 
% cfg = [];
% cfg.spmversion     = 'spm12';
% % cfg.resolution = 1;
% cfg.dim = MIDA.dim;
% mri = ft_volumereslice(cfg,mri);
% 
% 
% 
% %check if warping worked
% 
% ft_determine_coordsys(mri, 'interactive', 'no')
% ft_determine_coordsys(copyMIDA, 'interactive', 'no')

%% use warped MIDA to define segmented tissues

for i = 1:length(Tissues)
    Segment = ismember(MIDA.anatomy,eval(Tissues{i}));
    assignin('base', Tissues{i}, Segment);
end


Seg = struct;
Seg.skin = skin;
Seg.eyes = eyes;
Seg.muscle = muscle;
Seg.fat = fat;
Seg.spongybone = spongybone;
Seg.compactbone = compactbone;
Seg.gray = gray;
Seg.cerebellargray = cerebellargray;
Seg.white = white;
Seg.cerebellarwhite = cerebellarwhite;
Seg.csf = csf;
Seg.brainstem = brainstem;

Seg.dim = MIDA.dim;
Seg.coordsys = MIDA.coordsys;
Seg.transform = MIDA.transform;
Seg.unit = MIDA.unit;
Seg.transformorig = MIDA.transformorig;

%% Seg Gray matter

grayMatter = ismember(MIDA.anatomy,SegGrayMatter);

SegGray = struct;
SegGray.gray = grayMatter;
SegGray.dim = MIDA.dim;
SegGray.coordsys = MIDA.coordsys;
SegGray.transform = MIDA.transform;
SegGray.unit = MIDA.unit;
SegGray.transformorig = MIDA.transformorig;

%% plot segments
% seg_i = ft_datatype_segmentation(Seg,'segmentationstyle','indexed');
% cfg              = [];
% cfg.funparameter = 'seg';
% cfg.funcolormap  = lines(12); % distinct color per tissue
% cfg.location     = 'center';
% cfg.atlas        = seg_i;    % the segmentation can also be used as atlas 
% ft_sourceplot(cfg, seg_i);


%% compute hexahedral meshes
cfg = [];
cfg.spmversion = 'spm12';
cfg.tissues = Tissues;
cfg.shift  = 0;
cfg.method = 'hexahedral';
cfg.downsample = 1;
cfg.smooth = 'no';
mesh = ft_prepare_mesh(cfg,Seg);

cfg = [];
cfg.spmversion = 'spm12';
cfg.tissues = 'gray';
cfg.shift  = 0;
cfg.method = 'hexahedral';
cfg.downsample = 1;
cfg.smooth = 'no';
meshGray = ft_prepare_mesh(cfg,SegGray);

% ft_plot_mesh(mesh,'surfaceonly','yes')

%% headmodel
cfg        = [];
cfg.method ='simbio';
cfg.conductivity = [0.4348 0.5 0.1 0.04 0.04 0.0063 0.3333 0.2564 0.1429 0.1099 1.5385 0.1538];   % order follows mesh.tissuelabel
vol = ft_prepare_headmodel(cfg, mesh);   

cfg        = [];
cfg.method ='simbio';
cfg.conductivity = [0.3333];
volGray = ft_prepare_headmodel(cfg, meshGray);

save([save_folder '\headmodel_12T_FEM.mat'],'vol','-v7.3');
save([save_folder '\headmodel_12T_FEM_gray-only.mat'],'volGray','-v7.3');

%% Read sensor locations (Currently somewhat of a hack)
[loc_file,loc_path] = uigetfile('*.*','Feed me the EEG electrode locations');
path_loc_spec = [loc_path,loc_file];
path_loc_std = [path_ft '\template\electrode\standard_1005.elc'];
path_loc_attload = path_loc_spec; %Path to file
elec_attload = ft_read_sens(path_loc_attload);
elec_attload.label{strcmpi(elec_attload.label,'AF7')} = 'AFz';
elec_attload.label{strcmpi(elec_attload.label,'AF8')} = 'FCz';
idx_eog = cellfun(@(x)contains(x,'EOG','IgnoreCase',true),elec_attload.label);

%Channels to keep
channels = ['LPA';'RPA';'Nz'; elec_attload.label(~idx_eog)];

elec = ft_read_sens(path_loc_std); %Read layout file

elec = ft_determine_coordsys(elec);

idx_keep = ismember(upper(elec.label),upper(channels)); %Find indices of channels to keep
%Updating fields

elec.chanpos = elec.chanpos(idx_keep,:);
elec.chantype = elec.chantype(idx_keep);
elec.chanunit = elec.chanunit(idx_keep);
elec.elecpos = elec.elecpos(idx_keep,:);
elec.label = elec.label(idx_keep);

realignIdx = [];
for i = 1:numel(elec.label)
    realignIdx = [realignIdx, find(strcmp(upper(elec.label(:)),upper(channels(i))))];
end

elec.label = elec.label(realignIdx);
elec.chanpos = elec.chanpos(realignIdx,:);
elec.chantype = elec.chantype(realignIdx);
elec.chanunit = elec.chanunit(realignIdx);
elec.elecpos = elec.elecpos(realignIdx,:);

logFile.EEGchanlabels = elec.label;

%% Plot electrode positions before realignment 

figure
hold on
ft_plot_mesh(mesh,'surfaceonly','yes','vertexcolor','none','edgecolor','none','facecolor',[0.5 0.5 0.5],'face alpha',0.5)

ft_plot_sens(elec,'label','label'); 

%% Align electrode positions to head model

% %Find landmark-positions(nas, lpa, and rpa) in mri.hdr.fiducial.mri 
% %(WARNING: Might have to use mri.cfg.fiducial for other files!!!)
% nas = MIDA.cfg.fiducial.nas;
% lpa = MIDA.cfg.fiducial.lpa;
% rpa = MIDA.cfg.fiducial.rpa;
% 
% %Get fiducial positions in the ctf coordinate system using the
% %transformation matrix of the mri and the ft_warp_apply function
% transm = MIDA.transform;
% nas = ft_warp_apply(transm,nas,'homogenous');
% lpa = ft_warp_apply(transm,lpa,'homogenous');
% rpa = ft_warp_apply(transm,rpa,'homogenous');
% 
% %Create a structure similar to a template set of electrodes
% fid.elecpos = [nas; lpa; rpa]; %New coordinates of landmarks
% fid.label = {'Nz','LPA','RPA'}; %Same labels as in elec 
% fid.unit = 'mm'; %Same units as MRI
%  
% %Align electrodes based on landmarks
% cfg = [];
% cfg.method = 'fiducial';            
% cfg.target = fid; %See above
% cfg.elec = elec;
% cfg.fiducial = {'Nz','LPA','RPA'}; %Labels of landmarks in fid and elec
% elec_aligned = ft_electroderealign(cfg);

elec_aligned = elec;
%% Plot electrode positions after realignment 

figure;
hold on;
ft_plot_mesh(mesh,'surfaceonly','yes','vertexcolor','none','edgecolor','none','facecolor',[0.5 0.5 0.5],'face alpha',0.5)
% camlight
% electrodes
ft_plot_sens(elec_aligned,'label','label');
title('Electrode positions after realignment');


save([save_folder '\headmodel_12T_FEM_prepared_sens_vol.mat'],'vol','-v7.3');
save([save_folder '\headmodel_12T_FEM_prepared_sens_elec.mat'],'elec_aligned','-v7.3');

%% Adjust final position of electrodes by eye
% this is not necessary when using Simbio since it takes the closest vertex
% of the outer skin as electrode positions...
% 
% cfg           = [];
% cfg.method    = 'interactive';
% cfg.elec      = elec_aligned;
% cfg.headshape = vol;
% elec_aligned  = ft_electroderealign(cfg);

%% prepere sens (seems to be done in ft_prepare_leadfield anyway...)

channels = ft_channelselection({'EEG','-LPA', '-RPA', '-Nz'}, elec_aligned);
[vol, elec_aligned] = ft_prepare_vol_sens(vol, elec_aligned, 'channel', channels);
[volGray, elec_aligned] = ft_prepare_vol_sens(volGray, elec_aligned, 'channel', channels);

save([save_folder '\headmodel_12T_FEM_prepared_sens_vol.mat'],'vol','-v7.3');
save([save_folder '\headmodel_12T_FEM_prepared_sens_elec.mat'],'elec_aligned','-v7.3');

%% sourcemodel

cfg                = [];
cfg.mri            = SegGray;
cfg.grid.unit      ='mm';
cfg.grid.resolution = 7;
cfg.spmversion = 'spm12';
gridGray           = ft_prepare_sourcemodel(cfg);

ft_plot_mesh(gridGray.pos(gridGray.inside,:));

%% leadfield

cfg = [];
cfg.spmversion = 'spm12';
cfg.elec = elec_aligned;
cfg.headmodel = vol;
% cfg.resolution = 1;
cfg.grid.unit = 'mm'; %same unit as above
cfg.grid.pos = gridGray.pos;
cfg.grid.dim = gridGray.dim;
cfg.grid.inside = gridGray.inside;
cfg.normalize = 'yes'; %only if data should not be compared later against a baseline or other conditions!
[grid,cfg] = ft_prepare_leadfield(cfg);

save([save_folder '\leadfield_12T_FEM_gray-only.mat'],'grid','-v7.3');
save([save_folder '\LFcfg_12T_FEM.mat'],'cfg','-v7.3');
disp('saved leadfield and logFile');

