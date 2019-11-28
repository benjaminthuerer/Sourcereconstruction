function [] = UiO_inverse_model_eLORETA(MRI_data, MRI_path, EEG_data, EEG_path, LF_data, Head_data, LF_Head_path, path_ft, save_folder)


tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

path_fastICA = uigetdir([],'Give me the fastICA folder!');
path_ICASSO = uigetdir([],'Give me the ICASSO folder!');

%% load leadfield and EEG data
disp('Load MRI, leadfield, and headmodel. This may take a while');
mri = ft_read_mri([MRI_path MRI_data]);
disp('MRI loaded');

cfg = [];
cfg.spmversion = 'spm12';
cfg.method = 'interactive';
cfg.coordsys = 'acpc';
cfg.unit = 'mm';
[mri] = ft_volumerealign(cfg,mri);

load([LF_data],'grid');
disp('leadfield loaded')
load([Head_data],'vol');
disp('headmodel loaded')

% load EEG (EEGlab struct)
load([EEG_path EEG_data]);

srate = EEG.srate;

% filter to alpha band?
% EEG = pop_eegfiltnew(EEG, [], 8, [], true, [], 0);
% EEG = pop_eegfiltnew(EEG, 13, [], [], true, [], 0);

%% transform data for fieldtrip processing
addpath([path_ft '\external\eeglab']);

% just a hack for eeglab2fieldtrip
EEG.icachansind = randi(1,1,size(EEG.data,1));

fieldbox = 'timelockanalysis';
transform = 'none'; %or DIPTFIT transformation of channel locations
EEGdata = eeglab2fieldtrip(EEG, fieldbox,transform);


% this might have to be updated if frequency analysis are done
if length(size(EEG.data)) < 3
    EEGdata.dimord = 'chan_time';
else
    EEGdata.dimord = 'rpt_chan_time';
end

EEGdata = rmfield(EEGdata,'fsample');
EEGdata.label = EEGdata.label';

%% if continuous data: loop in 1s epochs through inverse model using eLORETA

% compute inverse model
cfg = [];
cfg.method = 'eloreta';
cfg.grid = grid;
% cfg.eloreta.normalize               = 'yes';
cfg.eloreta.lambda                  = 0.05;
cfg.headmodel = vol;
source = ft_sourceanalysis(cfg, EEGdata);

%% maybe you want to stop here and start a new function


%% some plotting  

% interpolate
cfg            = [];
cfg.spmversion = 'spm12';
%cfg.downsample = 20;
cfg.parameter = 'pow';
sourceInt  = ft_sourceinterpolate(cfg, source , mri);
    
% normalise to MNI otherwise can't be plotted
sourceIntNorm = ft_volumenormalise(cfg, sourceInt);
   
    % plot on brain
cfg = [];
cfg.spmversion = 'spm12';
cfg.method         = 'surface';
cfg.funparameter   = 'avg.pow';
cfg.maskparameter  = cfg.funparameter;
cfg.funcolormap    = 'jet';
cfg.opacitymap     = 'rampup';
cfg.projmethod     = 'nearest';
cfg.surffile       = 'surface_white_both.mat';
cfg.surfdownsample = 10;
ft_sourceplot(cfg, sourceIntNorm);
view ([90 0])
    
    % plot on MRI
cfg              = [];
cfg.spmversion = 'spm12';
cfg.method       = 'slice';
cfg.funparameter = 'pow';
cfg.maskparameter = cfg.funparameter;
cfg.opacitymap    = 'rampup';
ft_sourceplot(cfg,sourceIntNorm);
    
%%
%how do i get the power now for each sample?
tic;
lead = source.avg.mom(source.inside);
Rdata = zeros(length(lead),length(lead{1}(1,:)),'single');

percentages = floor(length(lead)/10:length(lead)/10:length(lead));
percentages2 = 10:10:100;
disp(['start computing the norm for each source and time bin. This may take a while...']);

i = 1;
while i < length(lead)+1
    Rdata(i,:) = abs(hilbert(sqrt(lead{i}(1,:).^2 + lead{i}(2,:).^2 + lead{i}(3,:).^2))); % transform to hilbert envelope
    if ~isempty(find(percentages == i))
        disp([num2str(percentages2(percentages == i)) ' % done']);
    end
    i = i+1;
end
toc;

clear lead
%% from here this is Liu 2017 specific to find the DMN!!!!!!!!!!!!
% norm each source at each time bin and compute the hilbert envelope
tic;
lead = source.avg.mom(source.inside);
Rdata = zeros(length(lead),length(lead{1}(1,:)),'single');
RdataPower = zeros(length(lead),length(lead{1}(1,:)),'single');

percentages = floor(length(lead)/10:length(lead)/10:length(lead));
percentages2 = 10:10:100;
disp(['start computing the norm for each source and time bin. This may take a while...']);

i = 1;
while i < length(lead)+1
    Rdata(i,:) = abs(hilbert(sqrt(lead{i}(1,:).^2 + lead{i}(2,:).^2 + lead{i}(3,:).^2))); % transform to hilbert envelope
    if ~isempty(find(percentages == i))
        disp([num2str(percentages2(percentages == i)) ' % done']);
    end
    i = i+1;
end
toc;

% resample / downsample / interpolate the data to 1 Hz
timeVector = 1:size(Rdata,2);
intervalVec = 1:srate:size(Rdata,2);
Rsampled = [];

for i = 1:size(Rdata,1)
    Rsampled(i,:) = interp1(timeVector,Rdata(i,:),intervalVec);
end

Rdata = Rsampled;
disp('Rdata is computed!');    
    

%% perform fastICA
addpath(path_fastICA)
addpath(path_ICASSO)

[ICAcalc] = icassoEst('both', Rdata, 10, 'g', 'tanh', 'approach', 'defl', 'lastEig', 30, 'maxNumIterations',1000); % 'randinit' instead of 'both'?
[rA] = icassoExp(ICAcalc);
[Iq, A, W, S] = icassoResult(rA);


%% correlate each IC with template DMN
%Unclear what is happening here, does the normalization work?
% % load and normalise template DMN
% mriTemplate = ft_read_mri('E:\Dropbox\Work\Matlab\DATA\johnlennon_data\halder\matlab_data\ForwardModelTests\PNAS_Smith09_rsn10_00004_DMN.nii');
% 
% ft_determine_coordsys(mriTemplate, 'interactive', 'no')
% 
% cfg = [];
% cfg.resolution = 1;
% cfg.dim = [256 256 256];
% mriT = ft_volumereslice(cfg,mriTemplate);
% 
% cfg            = [];
% cfg.spmversion = 'spm12';
% % cfg.downsample = 2;
% cfg.parameter = 'anatomy';
% MRIInt  = ft_sourceinterpolate(cfg, mriT, mriT);
% 
% path_SPM = uigetdir([],'Give me the SPM12 folder!');
% addpath(path_SPM)
% 
% MRIIntNorm = ft_volumenormalise(cfg, MRIInt);
% 
% 
% % correlate template with each IC
% 
% r = [];
% p = [];
% 
% for i = 1:size(A,2)
%     
%     Rdata = abs(A(:,i)); %abs? leave negative values?
%     source.time = 1;
%     source.avg.mom = [];
% 
% 
%     if numel(Rdata) > sum(grid.inside)
%         source.avg.pow(grid.inside) = mean(Rdata,2);
%     else
%         source.avg.pow(source.inside) = Rdata;
%     end
% 
%     cfg            = [];
%     cfg.spmversion = 'spm12';
% %     cfg.downsample = 2;
%     cfg.parameter = 'pow';
%     sourceInt  = ft_sourceinterpolate(cfg, source , mri);
% 
%     % normalize to template
%     sourceIntNorm = ft_volumenormalise(cfg, sourceInt);
%     
%     [r(i),p(i)] = corr(reshape(MRIIntNorm.anatomy,numel(MRIIntNorm.inside),1),reshape(sourceIntNorm.pow,numel(sourceIntNorm.pow),1));
% end

%% Interpolate and plot the relevant IC

%[~,icI] = max(r); % IC with highest correlation
%Just save all ICs for now
% put spatial IC data into source structure
for icI=1:size(A,2) %10 is enough for now
    Rdata = abs(A(:,icI)); % variable A is sICA
    source.time = 1;
    source.avg.mom = [];
    
    if numel(Rdata) > sum(grid.inside)
        source.avg.pow(grid.inside) = mean(Rdata,2);
    else
        source.avg.pow(source.inside) = Rdata;
    end
    
    % interpolate
    cfg            = [];
    cfg.spmversion = 'spm12';
    cfg.parameter = 'pow';
    sourceInt  = ft_sourceinterpolate(cfg, source , mri);
    
    % normalise to MNI otherwise can't be plotted
   sourceIntNorm = ft_volumenormalise(cfg, sourceInt);
   
    % plot on brain
%     cfg = [];
%     cfg.spmversion = 'spm12';
%     cfg.method         = 'surface';
%     cfg.funparameter   = 'avg.pow';
%     cfg.maskparameter  = cfg.funparameter;
%     cfg.funcolormap    = 'jet';
%     cfg.opacitymap     = 'rampup';
%     cfg.projmethod     = 'nearest';
%     cfg.surffile       = 'surface_white_both.mat';
%     cfg.surfdownsample = 10;
%     ft_sourceplot(cfg, sourceIntNorm);
%     view ([90 0])
    
    % plot on MRI
    cfg              = [];
    cfg.spmversion = 'spm12';
    cfg.method       = 'slice';
    cfg.funparameter = 'pow';
    cfg.maskparameter = cfg.funparameter;
    cfg.opacitymap    = 'rampup';
    ft_sourceplot(cfg,sourceIntNorm);
    
end
    % write to nifti for MRIcron plotting
    %save_folder = uigetdir([],'where do you want to save nifti of sources?');
    %%
    cfg = [];
    cfg.filename  = [save_folder '\ICs\source_activity' num2str(icI)];
    cfg.filetype  = 'nifti';
    cfg.parameter = 'pow';
    cfg.precision = 'single';
    ft_sourcewrite(cfg, sourceIntNorm);

end