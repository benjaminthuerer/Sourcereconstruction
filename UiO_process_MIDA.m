function [] = UiO_process_MIDA( MIDA, MIDApath, MIDAwrite )
% Remove air, downsample MIDA (which should contain filename and/or path of MIDA model, and location to write to disk and write the TPMs for segmentation
% This script loads and preprocesses MIDA and calculates tissue propability maps 
% (TPM) of the MIDA which will be used in a different script to normalize
% MIDA to subject space.

% load MIDA define TPMs
MIDA = ft_read_mri([MIDApath MIDA]);

cfg = [];
cfg.spmversion = 'spm12';
cfg.method = 'interactive';
cfg.coordsys = 'acpc';
cfg.unit = 'mm';
[MIDA] = ft_volumerealign(cfg,MIDA);

%% Make TPM segments to segment the subject

Tissues = {'Soft','Background','Bone','Gray','White','CSF'};

Soft = [33:35,37,39,51,85,86,55:59,38,42,60,61,63:84,88:96,98,43,62];
Background = [50];
Bone = [52,36,40,41,44:49,53,54,87];
Gray = [2,3:5,7,8,10,16,17,20,21,99,116];
White = [9, 12,18,22,23,100:115];
CSF = [6,32];

for i = 1:length(Tissues)
    Segment = ismember(MIDA.anatomy,eval(Tissues{i}));
    eval([Tissues{i} ' = Segment;']);
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

%% write a downsampled version of MIDA without background

MIDA_withoutair = MIDA;
MIDA_withoutair.anatomy(Seg.Background)=0;
cfg=[];
cfg.downsample=2;
cfg.method='nearest';
cfg.smooth='no';
cfg.spmversion='spm12';
MIDA_withoutair=ft_volumedownsample(cfg,MIDA_withoutair);

ft_write_mri([MIDAwrite 'MIDA_withoutair.nii'],MIDA_withoutair,'dataformat','nifti');


%% Write and smooth TPM segments

myTPM=MIDA;

myTPM.anatomy=zeros(myTPM.dim);
myTPM.anatomy(Seg.Gray)=0.999;
myTPM.anatomy(~Seg.Gray)=0.001;

cfg=[];
cfg.downsample=2;
cfg.smooth=500;
cfg.spmversion='spm12';
myTPM=ft_volumedownsample(cfg,myTPM);
ft_write_mri([MIDAwrite 'TPM_gray.nii'],myTPM,'dataformat','nifti');


myTPM=MIDA;

myTPM.anatomy=zeros(myTPM.dim);
myTPM.anatomy(Seg.White)=0.999;
myTPM.anatomy(~Seg.White)=0.001;

cfg=[];
cfg.downsample=2;
cfg.smooth=500;
cfg.spmversion='spm12';
myTPM=ft_volumedownsample(cfg,myTPM);
ft_write_mri([MIDAwrite 'TPM_white.nii'],myTPM,'dataformat','nifti');


myTPM=MIDA;

myTPM.anatomy=zeros(myTPM.dim);
myTPM.anatomy(Seg.Soft)=0.999;
myTPM.anatomy(~Seg.Soft)=0.001;

cfg=[];
cfg.downsample=2;
cfg.smooth=500;
cfg.spmversion='spm12';
myTPM=ft_volumedownsample(cfg,myTPM);
ft_write_mri([MIDAwrite 'TPM_soft.nii'],myTPM,'dataformat','nifti');


myTPM=MIDA;

myTPM.anatomy=zeros(myTPM.dim);
myTPM.anatomy(Seg.Bone)=0.999;
myTPM.anatomy(~Seg.Bone)=0.001;

cfg=[];
cfg.downsample=2;
cfg.smooth=500;
cfg.spmversion='spm12';
myTPM=ft_volumedownsample(cfg,myTPM);
ft_write_mri([MIDAwrite 'TPM_bone.nii'],myTPM,'dataformat','nifti');


myTPM=MIDA;

myTPM.anatomy=zeros(myTPM.dim);
myTPM.anatomy(Seg.CSF)=0.999;
myTPM.anatomy(~Seg.CSF)=0.001;

cfg=[];
cfg.downsample=2;
cfg.smooth=500;
cfg.spmversion='spm12';
myTPM=ft_volumedownsample(cfg,myTPM);
ft_write_mri([MIDAwrite 'TPM_CSF.nii'],myTPM,'dataformat','nifti');


myTPM=MIDA;

myTPM.anatomy=ones(myTPM.dim)*0.999;
myTPM.anatomy(~Seg.Background)=0.001;

cfg=[];
cfg.downsample=2;
cfg.smooth=500;
cfg.spmversion='spm12';
myTPM=ft_volumedownsample(cfg,myTPM);
ft_write_mri([MIDAwrite 'TPM_background.nii'],myTPM,'dataformat','nifti');
disp('Done. All TPMs and MIDA are written into nifti.');
end



