%% Source reconstruction

%% Constants and options
path_ft = uigetdir([],'Give me the Field Trip folder!');
path_mri = uigetfile('','Feed the MRI data [enable all data]');
path_loc_spec = uigetfile('','Feed me the EEG electrode locations [enable all data]');

%MRI alignment
met_align = 'interactive'; % method

%MRI segmentation
tissue_seg = {'brain','skull','scalp'}; % tissue types

%MRI mesh preparation
tissue_mesh = {'brain','skull','scalp'}; % tissue types
vertices_mesh = [3000 2000 1000];

%Head model preparation
met_head = 'singleshell'; % You can also specify 'openmeeg', 'bemcp', or another method.

%Sensor locations fieldtrip
path_loc_std = [path_ft '\template\electrode\standard_1005.elc'];


%% Add path to field trip and set general settings
addpath(path_ft)
ft_defaults

%% Read MRI file 
mri = ft_read_mri(path_mri);

%% Inspect the axes, align the anatomical MRI to desired coordinate system/template, and inspect again if correct
ft_determine_coordsys(mri, met_align, 'no')

cfg = [];
cfg.method = met_align;
[mri] = ft_volumerealign(cfg,mri);

ft_determine_coordsys(mri, met_align, 'no')

%% Segment MRI into different tissue types
cfg = [];
cfg.output = tissue_seg;
segmentedmri = ft_volumesegment(cfg, mri);

%% Prepare mesh
cfg = [];
cfg.tissue = tissue_mesh;
cfg.numvertices = vertices_mesh;
bnd=ft_prepare_mesh(cfg,segmentedmri);

%% convert brain into MNE
bnd_MNE = ft_transform_geometry(segmentedmri.transform, bnd(1)); %this is used later for Leadfield

%% Prepare head model
cfg = [];
cfg.method = met_head;
vol = ft_prepare_headmodel(cfg, bnd);

%% Plot segmented head model

% figure;
% ft_plot_mesh(vol.bnd(3),'facecolor','none'); %scalp
% figure;
% ft_plot_mesh(vol.bnd(2),'facecolor','none'); %skull
% figure;
% ft_plot_mesh(vol.bnd(1),'facecolor','none'); %brain

figure;
ft_plot_mesh(vol.bnd(3), 'facecolor',[0.2 0.2 0.2], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
hold on;
ft_plot_mesh(vol.bnd(2),'edgecolor','none','facealpha',0.4);
hold on;
ft_plot_mesh(vol.bnd(1),'edgecolor','none','facecolor',[0.4 0.6 0.4]);

%% Read sensor locations (Currently somewhat of a hack)

%Loading electrode location from experiment and from a FT location file
elec_specific = ft_read_sens(path_loc_spec);
elec_standard = ft_read_sens(path_loc_std);

%Hack for Andres amazing cap!!!!!! Please change this for a future experiment
elec_specific.label{strcmpi(elec_specific.label,'AF7')} = 'AFz';
elec_specific.label{strcmpi(elec_specific.label,'AF8')} = 'FCz';

%Finding indices of EOG channels
idx_eog = cellfun(@(x)contains(x,'EOG','IgnoreCase',true),elec_specific.label);

%Finding indices of landmarks and experiment electrodes in FT location file 
idx_elec = ismember(upper(elec_standard.label),upper(['LPA';'RPA';'Nz';elec_specific.label(~idx_eog)]));

%Updating all the electrode info based on the FT location file
elec = [];
elec.chanpos = elec_standard.chanpos(idx_elec,:);
elec.chantype = elec_standard.chantype(idx_elec);
elec.chanunit = elec_standard.chanunit(idx_elec);
elec.elecpos = elec_standard.elecpos(idx_elec,:);
elec.label = elec_standard.label(idx_elec);
elec.unit = elec_standard.unit;
elec.type = elec_standard.type;

%% Plot electrode positions before realignment 

figure;
% head surface (scalp)
ft_plot_mesh(vol.bnd(3), 'edgecolor','none','facealpha',0.8,'facecolor',[0.6 0.6 0.8]); 
hold on;
% electrodes
ft_plot_sens(elec);
title('Electrode positions before realignment');

%% Align electrode positions to head model

%Find landmarks (nas, lpa, and rpa) positions in mri.hdr!
nas = mri.cfg.fiducial.nas;
lpa = mri.cfg.fiducial.lpa;
rpa = mri.cfg.fiducial.rpa;

%Transformation matrix
transm = mri.transform;

%Transform landmark coordinates
nas = ft_warp_apply(transm,nas, 'homogenous');
lpa = ft_warp_apply(transm,lpa, 'homogenous');
rpa = ft_warp_apply(transm,rpa, 'homogenous');

fid.elecpos = [nas; lpa; rpa]; %New coordinates of landmarks
fid.label = {'Nz','LPA','RPA'}; %Same labels as in elec 
fid.unit = 'mm'; %Same units as MRI
 
%Align electrodes based on landmarks
cfg = [];
cfg.method = 'fiducial';            
cfg.target = fid; %See above
cfg.elec = elec;
cfg.fiducial = {'Nz', 'LPA','RPA' }; %Labels of landmarks in fid and elec
elec_aligned = ft_electroderealign(cfg);

%% Plot electrode positions after realignment 

figure;
% head surface (scalp)
ft_plot_mesh(vol.bnd(3), 'edgecolor','none','facealpha',0.8,'facecolor',[0.6 0.6 0.8]); 
hold on;
% electrodes
ft_plot_sens(elec_aligned);
title('Electrode positions after realignment');

%% Adjust final position of electrodes by eye

cfg           = [];
cfg.method    = 'interactive';
cfg.elec      = elec_aligned;
cfg.headshape = vol.bnd(3);
elec_aligned  = ft_electroderealign(cfg);

%% now compute leadfield matrix!!!!!!!!!!!!

cfg            = [];
cfg.elec       = elec;
cfg.vol        = vol;
cfg.grid.resolution = 6;  
cfg.grid.unit       = 'mm';% same unit as above, i.e. in cm
cfg.channel         = {'EEG','-LPA', '-RPA', '-Nz'};
% cfg.normalize = 'yes'; %only if data should not be compared later against a baseline or other conditions!
[grid,cfg] = ft_prepare_leadfield(cfg);

% transform grid.leadfield into leadfield matrix
LFM = [];
LFMx = [];
LFMy = [];
LFMz = [];
for i = 1:length(grid.leadfield)
    if ~isempty(grid.leadfield{i})
        LFMx(:,end+1) = grid.leadfield{i}(:,1);
    end
end
for i = 1:length(grid.leadfield)
    if ~isempty(grid.leadfield{i})
        LFMy(:,end+1) = grid.leadfield{i}(:,1);
    end
end
for i = 1:length(grid.leadfield)
    if ~isempty(grid.leadfield{i})
        LFMz(:,end+1) = grid.leadfield{i}(:,1);
    end
end

[LFM] = [LFMx,LFMy,LFMz];

% generates one matrix containing one number for each interaction between
% the sources and the electrodes
Ls = sqrt(LFMx.^2 + LFMy.^2 + LFMz.^2);



