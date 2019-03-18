%% Source reconstruction

%% Constants and options
path_ft = uigetdir([],'Give me the Field Trip folder!');
[mri_data,mri_path] = uigetfile('*.*','Feed the MRI data');
[loc_file,loc_path] = uigetfile('*.*','Feed me the EEG electrode locations');
save_folder = uigetdir([],'where do you want to save MRI, headmodel, and leadfield?');

path_mri = [mri_path,mri_data];
path_loc_spec = [loc_path,loc_file];

%MRI alignment
met_align = 'interactive'; % method

%MRI segmentation
tissue_seg = {'white','gray','csf','skull','scalp'};
vertices_mesh = [50000 50000 50000 50000 50000];
compress_mesh = [2000 2000 2400 1200 1000];


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


%% Read MRI file 
mri = ft_read_mri(path_mri);
    
cfg = [];
cfg.spmversion     = 'spm12';
cfg.resolution = 1;
cfg.dim = [256 256 256];
mri = ft_volumereslice(cfg,mri);

logFile.dim = cfg.dim;

% ft_determine_coordsys(mri, 'interactive', 'no')

cfg = [];
cfg.spmversion     = 'spm12';
cfg.method = met_align;
cfg.coordsys = 'ctf';
[mri] = ft_volumerealign(cfg,mri);
logFile.coordsys = cfg.coordsys;

%warp to MNI not working!!! leads to bad segmentations!
% cfg = [];
% cfg.spmversion     = 'spm12';
% mri = ft_volumenormalise(cfg, mri);

disp('!save MRI. This will take a while!');

save([save_folder '\MRI_MNI.mat'],'mri','-v7.3');

disp('MRI saved!');

%% Segment MRI into 5 different tissue types

cfg = [];
cfg.spmversion     = 'spm8';
cfg.spmmethod      = 'mars'; %generate 6 tissues
cfg.write          = 'no';
cfg.output = tissue_seg;
segmentedmri = ft_volumesegment(cfg, mri);

logFile.spmversion = cfg.spmversion;
logFile.segmentmethod = cfg.spmmethod;
logFile.segments = cfg.output;


% segment only skull because this runs better alone!
cfg = [];
cfg.spmversion     = 'spm8';
cfg.spmmethod      = 'mars'; %generate 6 tissues
cfg.write          = 'no';
cfg.output = {'skull'};
segmentedmriSkull = ft_volumesegment(cfg, mri);


%% Prepare mesh

% add the seperate segmented skull to the other tissues
% copy = segmentedmri;
% segmentedmri = copy;
segmentedmri.skull = segmentedmriSkull.skull;
% 
% segmentedmri.white = smooth3(segmentedmri.white);


% segmentedmri.gray = segmentedmriGray.gray;

% for better mesh computation, fill out each tissue and make a smooth
% surface around it
segmentedmri.scalp = imfill(segmentedmri.scalp,'holes');
segmentedmri.scalp = imdilate(segmentedmri.scalp,strel_bol(3));
segmentedmri.scalp = imfill(segmentedmri.scalp,'holes');

segmentedmri.skull = imfill(segmentedmri.skull,'holes');
segmentedmri.skull = imdilate(segmentedmri.skull,strel_bol(2));
segmentedmri.skull = imfill(segmentedmri.skull,'holes');

segmentedmri.csf = imfill(segmentedmri.csf,'holes');
segmentedmri.csf = imdilate(segmentedmri.csf,strel_bol(3));
segmentedmri.csf = imfill(segmentedmri.csf,'holes');

segmentedmri.gray = imfill(segmentedmri.gray,'holes');
segmentedmri.gray = imdilate(segmentedmri.gray,strel_bol(1));
segmentedmri.gray = imfill(segmentedmri.gray,'holes');

segmentedmri.white = imfill(segmentedmri.white,'holes');
segmentedmri.white = imdilate(segmentedmri.white,strel_bol(0));
segmentedmri.white = imfill(segmentedmri.white,'holes');


% prepare mesh for 5 tissues with high number of vertices (better results)
cfg = [];
cfg.tissue = tissue_seg;
% cfg.tissue = {'cortex','csf','skull','scalp'};
cfg.numvertices = vertices_mesh;
% cfg.numvertices = 2800;
cfg.spmversion   = 'spm12';
cfg.method = 'iso2mesh';
% cfg.radbound = 1;
% cfg.maxsurf = 0;
bnd = ft_prepare_mesh(cfg,segmentedmri); %other methods 'iso2mesh' and 'isosurface' give better results but do not work with openmeeg
logFile.numvertices = cfg.numvertices;

%%

% compress number of vertices (otherwise matlab runs out of memory)
[bnd(1).pos, bnd(1).tri] = meshresample(bnd(1).pos, bnd(1).tri,compress_mesh(1)/size(bnd(1).pos,1));
[bnd(2).pos, bnd(2).tri] = meshresample(bnd(2).pos, bnd(2).tri,compress_mesh(2)/size(bnd(2).pos,1));
[bnd(3).pos, bnd(3).tri] = meshresample(bnd(3).pos, bnd(3).tri,compress_mesh(3)/size(bnd(3).pos,1));
[bnd(4).pos, bnd(4).tri] = meshresample(bnd(4).pos, bnd(4).tri,compress_mesh(4)/size(bnd(4).pos,1));
[bnd(5).pos, bnd(5).tri] = meshresample(bnd(5).pos, bnd(5).tri,compress_mesh(5)/size(bnd(5).pos,1));

% [bnd(1).pos, bnd(1).tri] = meshresample(bnd(1).pos, bnd(1).tri,1800/size(bnd(1).pos,1));
% [bnd(2).pos, bnd(2).tri] = meshresample(bnd(2).pos, bnd(2).tri,1800/size(bnd(2).pos,1));
% [bnd(3).pos, bnd(3).tri] = meshresample(bnd(3).pos, bnd(3).tri,1800/size(bnd(3).pos,1));
% [bnd(4).pos, bnd(4).tri] = meshresample(bnd(4).pos, bnd(4).tri,1800/size(bnd(4).pos,1));
% [bnd(5).pos, bnd(5).tri] = meshresample(bnd(5).pos, bnd(5).tri,1800/size(bnd(5).pos,1));
logFile.compressvertices = compress_mesh;

%%
% check mesh quality
figure;
ft_plot_mesh(bnd(1),'edgecolor','none', 'facecolor',[0.4 0.2 0.4], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
hold on
ft_plot_mesh(bnd(2),'edgecolor','none','facecolor',[0.4 0.3 0.4], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
ft_plot_mesh(bnd(3),'facecolor',[0.6 0.4 0.2], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
ft_plot_mesh(bnd(4),'edgecolor','none','facecolor',[0.2 0.8 0.8], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
ft_plot_mesh(bnd(5),'edgecolor','none','facecolor',[0.1 0.1 0.1], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);

%% decouple intersected surfaces an repair mesh if necessery (if data quality is bad, this step has to be done several times)

bnd_idx = 1:length(bnd);
bndCopy = [];
for i = 1:length(bnd)
    bndCopy = [bndCopy bnd((length(bnd)+1)-i)];
end
bnd = decouplesurf(bnd);
bnd = [];
for i = 1:length(bndCopy)
    bnd = [bnd bndCopy((length(bndCopy)+1)-i)];
end

for i = 1:4
    [bnd(i).pos,bnd(i).tri] = meshcheckrepair(bnd(i).pos,bnd(i).tri);
end

[newnode,newelem] = surfboolean(bnd(4).pos,bnd(4).tri,'decouple',bnd(3).pos,bnd(3).tri);
bnd(4).pos = newnode(:,1:3);
bnd(4).tri = newelem(:,1:3);
[newnode,newelem] = surfboolean(bnd(3).pos,bnd(3).tri,'decouple',bnd(2).pos,bnd(2).tri);
bnd(3).pos = newnode(:,1:3);
bnd(3).tri = newelem(:,1:3);
for i = 1:5
    [newnode,newelem] = surfboolean(bnd(2).pos,bnd(2).tri,'decouple',bnd(1).pos,bnd(1).tri);
    bnd(2).pos = newnode(:,1:3);
    bnd(2).tri = newelem(:,1:3);
end



%% create headmodel with conductivity using openmeeg
% if this runs into an error, this is most likely the case because surfaces
% are still intersecting. Retry previous section.
cfg = [];
cfg.method = met_head;
% cfg.tissue = tissue_seg;
cfg.tissue = {'cortex','csf','skull','scalp'};
% cfg.conductivity = conductivity;
cfg.conductivity = [0.3333 1.5385 0.0217 0.4348];
vol = ft_prepare_headmodel(cfg, bnd);

logFile.conductivity = cfg.conductivity;
logFile.fiducial = mri.cfg.fiducial;
logFile.transform = mri.transform;

disp('!save headmodel. This will take a while!');

save([save_folder '\headmodel_5T_openmeeg_MNI.mat'],'vol','-v7.3');

disp('Headmodel saved!');


%% Read sensor locations (Currently somewhat of a hack)

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

figure;
% head surface (scalp)
ft_plot_mesh(vol.bnd(1), 'edgecolor','none','facealpha',0.8,'facecolor',[0.6 0.6 0.8]); 
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
addpath('Z:\Matlab_Scripts\Fieldtrip\new_fieldtrip\external\openmeeg');

load(fullfile('Z:\Matlab_Scripts\Fieldtrip\new_fieldtrip\template\sourcemodel\standard_sourcemodel3d10mm'));
template_grid = sourcemodel;
clear sourcemodel;

cfg                = [];
cfg.grid.warpmni   = 'yes';
cfg.grid.template  = template_grid;
cfg.grid.nonlinear = 'yes';
cfg.mri            = mri;
cfg.grid.unit      ='mm';
grid               = ft_prepare_sourcemodel(cfg);

cfg = [];
cfg.spmversion = 'spm12';
cfg.elec = elec_aligned;
cfg.headmodel = vol;
% cfg.grid.resolution = resolution; %because resolution set for gray matter
cfg.grid.unit = 'mm'; %same unit as above
cfg.grid.pos = grid.pos;
cfg.grid.dim = grid.dim;
cfg.grid.inside = grid.inside;
cfg.channel = {'EEG','-LPA', '-RPA', '-Nz'};
cfg.normalize = 'yes'; %only if data should not be compared later against a baseline or other conditions!
[grid,cfg] = ft_prepare_leadfield(cfg);

%%


% define grid of gray matter only
cfg                = [];
cfg.mri            = segmentedmriGray;
cfg.grid.unit      ='mm';
cfg.grid.resolution = 8;
cfg.spmversion = 'spm12';
gridGray           = ft_prepare_sourcemodel(cfg);

cfg                = [];
cfg.mri            = segmentedmriWhite;
cfg.grid.unit      ='mm';
cfg.grid.resolution = 8;
cfg.spmversion = 'spm12';
gridWhite           = ft_prepare_sourcemodel(cfg);

figure;
hold on
ft_plot_mesh(gridGray.pos(gridGray.inside,:))
ft_plot_mesh(gridWhite.pos(gridWhite.inside,:),'vertexcolor','r')

%%
cfg                = [];
cfg.mri            = mri;
cfg.grid.unit      ='mm';
cfg.grid.resolution = 8;
cfg.spmversion = 'spm12';
cfg.grid.nonlinear = 'yes'; 
cfg.grid.warpmni = 'yes';
gridGray           = ft_prepare_sourcemodel(cfg);



cfg                = [];
cfg.mri            = segmentedmriGray;
cfg.grid.unit      ='mm';
cfg.grid.template = gridGray;
cfg.spmversion = 'spm12';
gridGray2           = ft_prepare_sourcemodel(cfg);
ft_plot_mesh(gridGray2.pos(gridGray2.inside,:))


cfg                = [];
cfg.mri            = segmentedmriWhite;
cfg.grid.unit      ='mm';
cfg.grid.template = gridGray;
cfg.spmversion = 'spm12';
gridWhite           = ft_prepare_sourcemodel(cfg);


%%
% How warp grid, vol, elec to mni space???


%compute leadfield
cfg = [];
cfg.spmversion = 'spm12';
cfg.elec = elec_aligned;
cfg.headmodel = vol;
% cfg.grid.resolution = resolution; %because resolution set for gray matter
cfg.grid.unit = 'mm'; %same unit as above
cfg.grid.pos = gridGray.pos;
cfg.grid.dim = gridGray.dim;
cfg.grid.inside = gridGray.inside;
cfg.channel = {'EEG','-LPA', '-RPA', '-Nz'};
cfg.normalize = 'yes'; %only if data should not be compared later against a baseline or other conditions!
[grid,cfg] = ft_prepare_leadfield(cfg);

logFile.leadfieldResolution = resolution;
logFile.unit = grid.unit;

save([save_folder '\leadfield_5T_openmeeg_MNI.mat'],'grid','-v7.3');
save([save_folder '\logFile_5tissue.mat'],'logFile','-v7.3');
save([save_folder '\LFcfg_5tissue.mat'],'cfg','-v7.3');

disp('saved leadfield and logFile');
