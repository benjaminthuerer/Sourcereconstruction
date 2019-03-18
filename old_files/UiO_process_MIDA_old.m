% This script loads and preprocesses MIDA and calculates tissue propability maps 
% (TPM) of the MIDA which will be used in a different script to normalize
% MIDA to subject space.

% load MIDA define TPMs
path_ft = uigetdir([],'Give me the Field Trip folder!');
save_folder = uigetdir([],'where do you want to save the MIDA template and TPM?');

restoredefaultpath
addpath(path_ft)
ft_defaults

% load MIDA and subject mri
MIDA = ft_read_mri('Z:\07_fNetworks_rest-state\MIDA_head-model\MIDAv1.0\MIDA_v1.0\MIDA_v1_voxels\MIDA_v1.nii');

bg = 50;
bg = ismember(MIDA.anatomy,bg);
MIDA.anatomy(bg) = 0;

%% realign, reslice, and write MIDA

cfg = [];
cfg.spmversion = 'spm12';
cfg.method = 'interactive';
cfg.coordsys = 'acpc';
cfg.unit = 'mm';        % before it was cm. This might be a problem!
[MIDA] = ft_volumerealign(cfg,MIDA);

cfg = [];
cfg.method = 'flip';
[MIDA] = ft_volumereslice(cfg, MIDA);

% plot
cfg=[]; ft_sourceplot(cfg, MIDA);

cfg=[];
cfg.downsample=2;
cfg.method='nearest';
cfg.smooth='no';
cfg.spmversion='spm12';
MIDA=ft_volumedownsample(cfg,MIDA);

%possibly also resample by same factor as TPMS here
ft_write_mri([save_folder '\MIDA_processed.nii'],MIDA,'dataformat','nifti');


%% Make TPM segments to segment the subject

Tissues = {'Soft','Background','Bone','Gray','White','CSF'};
 
Soft = [33:35,37,39,51,85,86,55:59,38,42,60,61,63:84,88:96,98,43,62];
Bone = [52,36,40,41,44:49,53,54,87];
Gray = [2,3:5,7,8,10,16,17,20,21,99,116];
White = [9, 12,18,22,23,100:115];
CSF = [6,32];
Background = [50];
 
for i = 1:length(Tissues)
    Segment = ismember(MIDA.anatomy,eval(Tissues{i}));
    assignin('base', Tissues{i}, Segment);
end
 
Seg = struct;
Seg.Soft = Soft;
Seg.Bone = Bone;
Seg.Gray = Gray;
Seg.White = White;
Seg.CSF = CSF;
Seg.Background = Background;

Seg.dim = size(MIDA.anatomy);
Seg.coordsys = 'acpc';
Seg.transform = MIDA.transform;
Seg.unit = 'mm';

%% Write and smoot TPM segments
myTPM=MIDA;
 
myTPM.anatomy=zeros(myTPM.dim);
myTPM.anatomy(Seg.Gray)=0.999;
myTPM.anatomy(~Seg.Gray)=0.001;
 
cfg=[];
cfg.smooth=500;
cfg.spmversion='spm12';
myTPM=ft_volumedownsample(cfg,myTPM);
ft_write_mri([save_folder '\TPM_gray.nii'],myTPM,'dataformat','nifti');
 
 
myTPM=MIDA;
 
myTPM.anatomy=zeros(myTPM.dim);
myTPM.anatomy(Seg.White)=0.999;
myTPM.anatomy(~Seg.White)=0.001;
 
cfg=[];
cfg.smooth=500;
cfg.spmversion='spm12';
myTPM=ft_volumedownsample(cfg,myTPM);
ft_write_mri([save_folder '\TPM_white.nii'],myTPM,'dataformat','nifti');
 
 
myTPM=MIDA;
 
myTPM.anatomy=zeros(myTPM.dim);
myTPM.anatomy(Seg.Soft)=0.999;
myTPM.anatomy(~Seg.Soft)=0.001;
 
cfg=[];
cfg.smooth=500;
cfg.spmversion='spm12';
myTPM=ft_volumedownsample(cfg,myTPM);
ft_write_mri([save_folder '\TPM_soft.nii'],myTPM,'dataformat','nifti');
 
 
myTPM=MIDA;
 
myTPM.anatomy=zeros(myTPM.dim);
myTPM.anatomy(Seg.Bone)=0.999;
myTPM.anatomy(~Seg.Bone)=0.001;
 
cfg=[];
cfg.smooth=500;
cfg.spmversion='spm12';
myTPM=ft_volumedownsample(cfg,myTPM);
ft_write_mri([save_folder '\TPM_bone.nii'],myTPM,'dataformat','nifti');
 
 
myTPM=MIDA;
 
myTPM.anatomy=zeros(myTPM.dim);
myTPM.anatomy(Seg.CSF)=0.999;
myTPM.anatomy(~Seg.CSF)=0.001;
 
cfg=[];
cfg.smooth=500;
cfg.spmversion='spm12';
myTPM=ft_volumedownsample(cfg,myTPM);
ft_write_mri([save_folder '\TPM_CSF.nii'],myTPM,'dataformat','nifti');
 
 
myTPM=MIDA;
 
myTPM.anatomy=ones(myTPM.dim)*0.999;
myTPM.anatomy(~Seg.Background)=0.001;
 
cfg=[];
cfg.smooth=500;
cfg.spmversion='spm12';
myTPM=ft_volumedownsample(cfg,myTPM);
ft_write_mri([save_folder '\TPM_background.nii'],myTPM,'dataformat','nifti');

disp('done. All TPMs and MIDA are written into nifti');