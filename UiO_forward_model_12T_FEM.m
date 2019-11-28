function [] = UiO_forward_model_12T_FEM(mri_data, mri_path, save_folder, loc_file, loc_path, path_ft)


%% load MIDA and segment
MIDA = ft_read_mri([mri_path mri_data]);

% cfg = [];
% cfg.spmversion     = 'spm12';
% cfg.resolution = 1;
% cfg.dim = [256 256 256];
% MIDA = ft_volumereslice(cfg,MIDA);

cfg = [];
cfg.spmversion = 'spm12';
cfg.method = 'interactive';
cfg.coordsys = 'acpc';
cfg.unit = 'mm';
[MIDA] = ft_volumerealign(cfg,MIDA);

% Define tissue types
Tissues = {'skin','eyes','muscle','fat','spongybone','compactbone','gray','cerebellargray','white','cerebellarwhite','csf','brainstem'};

%% define MIDA tissues OLD:
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

%% define MIDA tissues NEW:
% skin = [1,51];
% eyes = [55:59];
% muscle = [38,60,63:84,91:96];
% fat = [43,62];
% spongybone = [52:54];
% compactbone = [40];
% gray = [10];
% cerebellargray = [2];
% white = [12];
% cerebellarwhite = [9];
% csf = [6,32];
% brainstem = [11,14,15];

% special tissue for GrayMatter
SegGrayMatter = sort([gray cerebellargray]);

bg = 50;
bg = ismember(MIDA.anatomy,bg);
MIDA.anatomy(bg) = 0;


%% use warped MIDA to define segmented tissues

for i = 1:length(Tissues)
    Segment = ismember(MIDA.anatomy,eval(Tissues{i}));
    eval([Tissues{i} ' = Segment;']);
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


%% compute hexahedral meshes (make sure Simbio is in your path: fieldtrip/external/simbio !!!!)
cfg = [];
cfg.spmversion = 'spm12';
cfg.tissues = Tissues;
cfg.shift  = 0;
cfg.method = 'hexahedral';
cfg.downsample = 1;
%cfg.smooth = []; % If you define something in smooth, it will not give
%you hexahedral or tetrahedral!
disp('start preparing meshes');
mesh = ft_prepare_mesh(cfg,Seg);
disp('done');


%% headmodel
cfg        = [];
cfg.method ='simbio';
cfg.conductivity = [0.4348 0.5 0.1 0.04 0.04 0.0063 0.3333 0.2564 0.1429 0.1099 1.5385 0.1538];   % order follows mesh.tissuelabel
disp('start preparing headmodel');
vol = ft_prepare_headmodel(cfg, mesh);   

%% Read sensor locations (Currently somewhat of a hack)
path_loc_spec = [loc_path,loc_file];
path_loc_std = [path_ft '\template\electrode\standard_1005.elc'];
path_loc_attload = path_loc_spec; %Path to file
elec_attload = ft_read_sens(path_loc_attload);
disp('done');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% hard coded passage for UiO (flip EEG channels AF7...)!
%%%%%%%%%%%%%%%%%%%%%%%%%% 

% old for UiO
% hard coded channel flip for UiO:
elec_attload.label{strcmpi(elec_attload.label,'AF7')} = 'AFz';
elec_attload.label{strcmpi(elec_attload.label,'AF8')} = 'FCz';

idx_eog = cellfun(@(x)contains(x,'EO','IgnoreCase',true),elec_attload.label) + cellfun(@(x)contains(x,'EM','IgnoreCase',true),elec_attload.label);

channels = [elec_attload.label(~idx_eog)];

elec = ft_read_sens(path_loc_std); %Read layout file

idx_keep = ismember(upper(elec.label),upper(channels)); %Find indices of channels to keep
%Updating fields

elec.chanpos = elec.chanpos(idx_keep,:);
elec.chantype = elec.chantype(idx_keep);
elec.chanunit = elec.chanunit(idx_keep);
elec.elecpos = elec.elecpos(idx_keep,:);
elec.label = elec.label(idx_keep);

% realign channel order to subject channel order
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

elec_attload = elec;

%% new for Sebastian data

% idx_eog = cellfun(@(x)contains(x,{'EO', 'EM'},'IgnoreCase',true),elec_attload.label);
% 
% channels = [elec_attload.label(~idx_eog)];
% 
% elec = elec_attload;
% 
% 
% 
% idx_keep = ismember(upper(elec.label),upper(channels)); %Find indices of channels to keep
% %Updating fields
% 
% elec.chanpos = elec.chanpos(idx_keep,:);
% elec.chantype = elec.chantype(idx_keep);
% elec.chanunit = elec.chanunit(idx_keep);
% elec.elecpos = elec.elecpos(idx_keep,:);
% elec.label = elec.label(idx_keep);
% 
% elec_aligned = elec;
% elec_aligned = ft_convert_units(elec_aligned,'mm');
% 
% elec_aligned.elecpos = [elec_aligned.elecpos(:,2)*-1 elec_aligned.elecpos(:,1) elec_aligned.elecpos(:,3)];
% elec_aligned.chanpos = [elec_aligned.chanpos(:,2)*-1 elec_aligned.chanpos(:,1) elec_aligned.chanpos(:,3)];

% elec = ft_determine_coordsys(elec_aligned);
%%%%%%%


%% Adjust final position of electrodes by eye
% this is not necessary when using Simbio since it takes the closest vertex
% of the outer skin as electrode positions...
% 
cfg           = [];
cfg.method    = 'interactive';
cfg.elec      = elec_attload;
cfg.headshape = vol;
disp('start realign electrodes');
elec_attload  = ft_electroderealign(cfg);


%% prepere sens (seems to be done in ft_prepare_leadfield anyway...)

% channels = ft_channelselection({'EEG'}, elec_aligned);
disp('prepare volume and sensor matrix. This take a while...');
[vol, elec_attload] = ft_prepare_vol_sens(vol, elec_attload, 'channel', elec_attload.label);
%[volGray, elec_attload] = ft_prepare_vol_sens(volGray, elec_attload, 'channel', channels);

save([save_folder '\headmodel_12T_FEM_prepared_sens_vol.mat'],'vol','-v7.3');
save([save_folder '\headmodel_12T_FEM_prepared_sens_elec.mat'],'elec_attload','-v7.3');
disp('headmodel and electrodes saved');

% %% Plot electrode positions after realignment 
% figure;
% hold on;
% ft_plot_mesh(mesh,'surfaceonly','yes','vertexcolor','none','edgecolor','none','facecolor',[0.5 0.5 0.5],'face alpha',0.5)
% % camlight
% % electrodes
% ft_plot_sens(elec_attload,'label','label');
% title('Electrode positions after realignment');

%% sourcemodel
cfg                = [];
cfg.mri            = SegGray;
cfg.elec = elec_attload;
cfg.vol = vol;
cfg.grid.unit      ='mm';
cfg.grid.resolution = 6; %change resolution???
cfg.spmversion = 'spm12';
disp('start preparing sourcemodel');
gridGray           = ft_prepare_sourcemodel(cfg, vol, elec_attload);


%% leadfield
cfg = [];
cfg.spmversion = 'spm12';
cfg.elec = elec_attload;
cfg.headmodel = vol;
cfg.grid = gridGray;
% cfg.normalize = 'yes'; %only if data should not be compared later against a baseline or other conditions!
cfg.reducerank = 3; % Make sure to check that because it is different for MEG data!!!!????
disp('start preparing leadfield');
[grid,cfg] = ft_prepare_leadfield(cfg);

disp('save leadfield');
save([save_folder '\leadfield_12T_FEM_gray-only.mat'],'grid','-v7.3');
save([save_folder '\LFcfg_12T_FEM.mat'],'cfg','-v7.3');
disp('saved leadfield and logFile');
end

