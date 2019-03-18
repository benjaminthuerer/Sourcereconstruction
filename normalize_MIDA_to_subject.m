subjectMRIPath='C:\Users\sebashal\Dropbox\Work\Matlab\DATA\johnlennon_data\halder\matlab_data\AuditoryMRI\AMRI02\0007\b3337-0007-00001-000176-01.nii';
MIDAPath='C:\Users\sebashal\Dropbox\Work\Matlab\DATA\johnlennon_data\halder\matlab_data\MIDA\MIDA_v1_voxels\MIDA_v1.nii';
%% Read the MIDA
 
MIDA = ft_read_mri(MIDAPath);
 
cfg = [];
cfg.spmversion = 'spm12';
cfg.method = 'interactive';
cfg.coordsys = 'acpc';
cfg.unit = 'cm';
[MIDA] = ft_volumerealign(cfg,MIDA);
cfg = [];
cfg.method = 'flip';
[MIDA] = ft_volumereslice(cfg, MIDA);
cfg=[]; ft_sourceplot(cfg, MIDA);
 
 
 
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
    assignin('base', Tissues{i}, Segment);
end
 
 
 Seg = struct;
Seg.Soft = Soft;
% Seg.Eye = Eye;
% Seg.Muscle = Muscle;
% Seg.Fat = Fat;
% Seg.SpongyBone = SpongyBone;
 Seg.Bone = Bone;
 Seg.Gray = Gray;
% Seg.CerebellarGray = CerebellarGray;
 Seg.White = White;
% Seg.CerebellarWhite = CerebellarWhite;
 Seg.CSF = CSF;
  
 Seg.Background = Background;
% Seg.Brainstem = Brainstem;
 
Seg.dim = size(MIDA.anatomy);
Seg.coordsys = 'acpc';
Seg.transform = MIDA.transform;
Seg.unit = 'cm';
%% write a downsampled version of MIDA without air
 
MIDA_withoutair = MIDA;
MIDA_withoutair.anatomy(Seg.Background)=0;
cfg=[];
cfg.downsample=3;
cfg.method='nearest';
cfg.smooth='no';
cfg.spmversion='spm12';
MIDA_withoutair=ft_volumedownsample(cfg,MIDA_withoutair);
%possibly also resample by same factor as TPMS here
ft_write_mri('MIDA_withoutair.nii',MIDA_withoutair,'dataformat','nifti');
 
%% Write and smoot TPM segments
 
myTPM=MIDA;
 
myTPM.anatomy=zeros(myTPM.dim);
myTPM.anatomy(Seg.Gray)=0.999;
myTPM.anatomy(~Seg.Gray)=0.001;
 
cfg=[];
cfg.downsample=3;
cfg.smooth=500;
cfg.spmversion='spm12';
myTPM=ft_volumedownsample(cfg,myTPM);
ft_write_mri('TPM_gray.nii',myTPM,'dataformat','nifti');
 
 
myTPM=MIDA;
 
myTPM.anatomy=zeros(myTPM.dim);
myTPM.anatomy(Seg.White)=0.999;
myTPM.anatomy(~Seg.White)=0.001;
 
cfg=[];
cfg.downsample=3;
cfg.smooth=500;
cfg.spmversion='spm12';
myTPM=ft_volumedownsample(cfg,myTPM);
ft_write_mri('TPM_white.nii',myTPM,'dataformat','nifti');
 
 
myTPM=MIDA;
 
myTPM.anatomy=zeros(myTPM.dim);
myTPM.anatomy(Seg.Soft)=0.999;
myTPM.anatomy(~Seg.Soft)=0.001;
 
cfg=[];
cfg.downsample=3;
cfg.smooth=500;
cfg.spmversion='spm12';
myTPM=ft_volumedownsample(cfg,myTPM);
ft_write_mri('TPM_soft.nii',myTPM,'dataformat','nifti');
 
 
myTPM=MIDA;
 
myTPM.anatomy=zeros(myTPM.dim);
myTPM.anatomy(Seg.Bone)=0.999;
myTPM.anatomy(~Seg.Bone)=0.001;
 
cfg=[];
cfg.downsample=3;
cfg.smooth=500;
cfg.spmversion='spm12';
myTPM=ft_volumedownsample(cfg,myTPM);
ft_write_mri('TPM_bone.nii',myTPM,'dataformat','nifti');
 
 
myTPM=MIDA;
 
myTPM.anatomy=zeros(myTPM.dim);
myTPM.anatomy(Seg.CSF)=0.999;
myTPM.anatomy(~Seg.CSF)=0.001;
 
cfg=[];
cfg.downsample=3;
cfg.smooth=500;
cfg.spmversion='spm12';
myTPM=ft_volumedownsample(cfg,myTPM);
ft_write_mri('TPM_CSF.nii',myTPM,'dataformat','nifti');
 
 
myTPM=MIDA;
 
myTPM.anatomy=ones(myTPM.dim)*0.999;
myTPM.anatomy(~Seg.Background)=0.001;
%myTPM.anatomy(Seg.Background)=0.999;
%myTPM.anatomy(~Seg.CSF&~Seg.Bone&~Seg.Soft&~Seg.Gray&~Seg.White)=0.999;
 
cfg=[];
cfg.downsample=3;
cfg.smooth=500;
cfg.spmversion='spm12';
myTPM=ft_volumedownsample(cfg,myTPM);
ft_write_mri('TPM_background.nii',myTPM,'dataformat','nifti');
 
 
%% read a subject MRI
 
mriSubjectBase = ft_read_mri(mriSubjectPath); 
cfg = [];
cfg.spmversion = 'spm12';
cfg.method = 'interactive';
cfg.coordsys = 'acpc';
cfg.unit = 'mm';
[mriSubjectBase] = ft_volumerealign(cfg,mriSubjectBase);
 
cfg=[]; ft_sourceplot(cfg, mriSubjectBase);
ft_determine_coordsys(mriSubjectBase, 'interactive', 'no');
 
%% Coregister and reslice the subject to the MIDA
 
mriSubject = mriSubjectBase;
 
cfg = [];
cfg.method = 'flip';
[mriSubjectResliced] = ft_volumereslice(cfg, mriSubject);
cfg=[]; ft_sourceplot(cfg, mriSubjectResliced);
 
cfg = [];
cfg.method = 'nearest';
%cfg.dim = MIDA_withoutair.dim;
cfg.resolution = 1.5;
cfg.dim=round(mriSubjectResliced.dim.*mriSubjectResliced.hdr.volres./(MIDA_withoutair.hdr.volres*3));
% cfg.xrange = [-MIDA_withoutair.dim(3)/2*cfg.resolution MIDA_withoutair.dim(3)/2*cfg.resolution];
% cfg.yrange = [-MIDA_withoutair.dim(2)/2*cfg.resolution MIDA_withoutair.dim(2)/2*cfg.resolution];
% cfg.zrange = [-MIDA_withoutair.dim(1)/2*cfg.resolution MIDA_withoutair.dim(1)/2*cfg.resolution];
cfg.coordsys = 'acpc';
cfg.unit = 'mm';
[mriSubjectResliced] = ft_volumereslice(cfg, mriSubjectResliced);
%ft_determine_coordsys(mriSubjectResliced, 'interactive', 'no');
 
cfg=[]; ft_sourceplot(cfg, mriSubjectResliced);
 
 
 
 
cfg = [];
cfg.spmversion = 'spm12';
cfg.method = 'spm';
cfg.coordsys = 'acpc';
cfg.unit = 'mm';
cfg.parameter = 'anatomy';
cfg.viewresult = 'yes';
%cfg.spm.reslice = 'yes';
[mriSubjectRealigned] = ft_volumerealign(cfg,mriSubjectResliced,MIDA_withoutair);
 
cfg=[]; ft_sourceplot(cfg, mriSubjectRealigned);
 
 
%% Segment the subject using the MIDA TPMs
cfg = [];
cfg.spmversion = 'spm12';
cfg.spmmethod = 'new';
cfg.output = {'tpm'};
cfg.tpm = {'TPM_gray.nii', 'TPM_white.nii', 'TPM_CSF.nii', 'TPM_bone.nii', 'TPM_soft.nii', 'TPM_background.nii'};
opts.lkp = [1,1,2,2,3,3,4,4,4,4,5,5,6,6];
cfg.opts=opts;
cfg.write = 'yes';
cfg.name = 'segmentedSubject';
 
segmentedSubject = ft_volumesegment(cfg, mriSubjectRealigned);
 
%% Apply transformation to the MIDA to get it into subject space
[mri] = ft_volumenormalise(cfg, MIDA_withoutair);
%% Start making the mesh and so forth
 
