%% Constants and options
path_ft = uigetdir([],'Give me the Field Trip folder!');
[LF_Head_path] = uigetdir([],'Feed the leadfield and headmodel folder');
[EEG_file,EEG_path] = uigetfile('*.*','Feed me the EEG data');
[MRI_data,MRI_path] = uigetfile('*.*','Feed me the MRI data');
path_fastICA = uigetdir([],'Give me the fastICA folder!');
path_ICASSO = uigetdir([],'Give me the ICASSO folder!');

LF_data = [LF_Head_path,'\leadfield_12T_FEM_gray-only.mat'];
Head_data = [LF_Head_path,'\headmodel_12T_FEM_prepared_sens_vol.mat'];

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

restoredefaultpath
addpath(path_ft)
ft_defaults

mriTemplate = ft_read_mri('Z:\07_fNetworks_rest-state\data_attentional-load\DMN_template_melodic_IC_sum.nii');

mri = ft_read_mri([MRI_path MRI_data]);
disp('MRI loaded');

cfg = [];
cfg.spmversion = 'spm12';
cfg.method = 'interactive';
cfg.coordsys = 'acpc';
cfg.unit = 'mm';
[mri] = ft_volumerealign(cfg,mri);

cfg = [];
cfg.resolution = 1;
cfg.dim = [256 256 256];
mriT = ft_volumereslice(cfg,mriTemplate);

cfg            = [];
cfg.spmversion = 'spm12';
cfg.parameter = 'anatomy';
MRIInt  = ft_sourceinterpolate(cfg, mriT, mri);

MRIIntNorm = ft_volumenormalise(cfg, MRIInt);
%% load leadfield and EEG data
disp('Load MRI, leadfield, and headmodel. This may take a while');
load(LF_data,'grid');
disp('leadfield loaded')
load(Head_data,'vol');
disp('headmodel loaded')


load([EEG_path EEG_file]);

srate = EEG.srate;

% filter to alpha band?
% EEG = pop_eegfiltnew(EEG, [], 8, [], true, [], 0);
% EEG = pop_eegfiltnew(EEG, 13, [], [], true, [], 0);

%% transform data for fieldtrip processing
addpath('Z:\Matlab_Scripts\Fieldtrip\new_fieldtrip\external\eeglab')

EEG.icachansind = 1:size(EEG.data,1);

fieldbox = 'timelockanalysis';
transform = 'none'; %or DIPTFIT transformation of channel locations
EEGdata = eeglab2fieldtrip(EEG, fieldbox,transform);

EEGdata.dimord = 'chan_time';

% EEGdata.label = EEGdata.label'; %check if this is ok!

%% eloreta for each 1s epoch
EpochLength = 1; %epoch length in seconds
Elength = EpochLength*srate;
IElag = 0; % inter-epoch lag (how much does it overlap). e.g. 0.4 is a overlap of 40%!
Plength = Elength-(Elength*IElag);
    
Rdata = zeros(sum(grid.inside),floor(length(EEGdata.avg)/Plength));

cfg = [];
cfg.method = 'eloreta';
cfg.grid = grid;
cfg.headmodel = vol;

percentages = floor(length(EEGdata.avg)/Plength)/10:floor(length(EEGdata.avg)/Plength)/10:floor(length(EEGdata.avg)/Plength);
percentages2 = 10:10:100;
disp(['start calculating eloreta for ' num2str(floor(length(EEGdata.avg)/Plength)) ' epochs']);

EEGcopy = EEGdata;
tic;
i = 1;
while i < floor(length(EEGcopy.avg)/Plength)
    EEGdata.avg = double(EEGcopy.avg(:,(i-1)*Plength+1:(i-1)*Plength+Elength));
    EEGdata.var = double(EEGcopy.var(:,(i-1)*Plength+1:(i-1)*Plength+Elength));
    EEGdata.time = 1:Elength;
    
    source = ft_sourceanalysis(cfg, EEGdata);
    
    lead = source.avg.mom(source.inside);
    ii = 1;
    while ii < length(lead)+1
       Rdata(ii,i) = mean(abs(hilbert(sqrt(lead{ii}(1,:).^2 + lead{ii}(2,:).^2 + lead{ii}(3,:).^2))),2);
        ii = ii+1;
    end

    if ~isempty(find(percentages == i, 1))
        disp([num2str(percentages2(percentages == i)) ' % done']);
    end
    i = i+1;
end
toc;

copyR = Rdata;
overall_r = [];

addpath(path_fastICA)
addpath(path_ICASSO)

for ICs = 10:60
    Rdata = copyR;
    [ICAcalc] = icassoEst('both', Rdata, 10, 'g', 'tanh', 'approach', 'defl', 'lastEig', ICs, 'maxNumIterations',1000); % 'randinit' instead of 'both'?
    [rA] = icassoExp(ICAcalc);
    [Iq, A, W, S] = icassoResult(rA);
    
    r = [];
    p = [];

    for i = 1:size(A,2)

        Rdata = A(:,i); %abs? leave negative values?
        source.time = 1;
        source.avg.mom = [];
        source.avg.ori = [];

        source.avg.pow(source.inside) = Rdata;

        cfg            = [];
        cfg.spmversion = 'spm12';
        cfg.parameter = 'pow';
        sourceInt  = ft_sourceinterpolate(cfg, source , mri);

        % normalize to template
        sourceIntNorm = ft_volumenormalise(cfg, sourceInt);

        [r(i),p(i)] = corr(reshape(MRIIntNorm.anatomy,numel(MRIIntNorm.inside),1),reshape(sourceIntNorm.pow,numel(sourceIntNorm.pow),1));
    end
    
    [overall_r(ICs,1),overall_r(ICs,2)] = max(abs(r));
end

% cfg = [];
% cfg.spmversion = 'spm12';
% cfg.method         = 'surface';
% cfg.funparameter   = 'anatomy';
% cfg.maskparameter  = cfg.funparameter;
% cfg.funcolormap    = 'jet'; 
% cfg.opacitymap     = 'rampup';  
% cfg.projmethod     = 'nearest'; 
% cfg.surffile       = 'surface_white_both.mat';
% cfg.surfdownsample = 10; 
% ft_sourceplot(cfg, sourceIntNorm);
% view ([90 0])
% 
%     cfg              = [];
%     cfg.spmversion = 'spm12';
%     cfg.method       = 'slice';
%     cfg.funparameter = 'pow';
%     cfg.maskparameter = cfg.funparameter;
%     % cfg.funcolorlim   = [min(min(min(sourceIntNorm.pow))) max(max(max(sourceIntNorm.pow)))];
%     % cfg.opacitylim    = [min(min(min(sourceIntNorm.pow))) max(max(max(sourceIntNorm.pow)))]; 
%     cfg.opacitymap    = 'rampup'; 
%     ft_sourceplot(cfg,sourceIntNorm);
%% if continuous data: loop in 1s epochs through inverse model using eLORETA
% addpath('Z:\Matlab_Scripts\Fieldtrip\new_fieldtrip\external\openmeeg');

% cfg = [];
% cfg.method = 'eloreta';
% cfg.grid = grid;
% % cfg.eloreta.keepfilter              = 'yes';
% % cfg.eloreta.normalize               = 'yes';
% % cfg.eloreta.lambda                  = 0.05;
% % cfg.eloreta.projectnoise            = 'yes';
% cfg.headmodel = vol;
% 
% source = ft_sourceanalysis(cfg, EEGdata);
% 
% tic;
% lead = source.avg.mom(source.inside);
% Rdata = zeros(length(lead),length(lead{1}(1,:)),'single');
% 
% percentages = floor(length(lead)/10:length(lead)/10:length(lead));
% percentages2 = 10:10:100;
% disp(['start computing the norm for each source and time bin. This may take a while...']);
% 
% i = 1;
% while i < length(lead)+1
%     Rdata(i,:) = abs(hilbert(sqrt(lead{i}(1,:).^2 + lead{i}(2,:).^2 + lead{i}(3,:).^2))); % transform to hilbert envelope
%     if ~isempty(find(percentages == i))
%         disp([num2str(percentages2(percentages == i)) ' % done']);
%     end
%     i = i+1;
% end
% toc;
% 
% 
% timeVector = 1:size(Rdata,2);
% intervalVec = 1:srate:size(Rdata,2);
% Rsampled = [];
% 
% for i = 1:size(Rdata,1)
%     Rsampled(i,:) = interp1(timeVector,Rdata(i,:),intervalVec);
% end
% 
% Rdata = Rsampled;
% disp('Rdata is computed!');    
%     
%     
%     
% % if ndims(EEG.data)==3
% %     cfg.feedback = 'yes';  
% %     EpochLength = 1; %epoch length in seconds
% %     Elength = EpochLength*srate;
% %     IElag = 0; % inter-epoch lag (how much does it overlap). e.g. 0.4 is a overlap of 40%!
% %     Plength = srate-(srate*IElag);
% %     
% %     Rdata = zeros(sum(cfg.grid.inside),floor(length(EEGdata.avg)/Plength));
% %     lead = source.avg.mom(source.inside);
% %     
% %     percentages = floor(length(EEGdata.avg)/Plength)/10:floor(length(EEGdata.avg)/Plength)/10:floor(length(EEGdata.avg)/Plength);
% %     percentages2 = 10:10:100;
% %     disp(['start calculating norm for ' num2str(floor(length(EEGdata.avg)/Plength)) ' epochs']);
% %     tic;
% %     i = 1;
% %     while i < floor(length(EEGdata.avg)/Plength)
% %     
% %         for ii = 1:numel(lead)
% %             Rdata(ii,i) = mean(sqrt(lead{ii}(1,(i-1)*Plength+1:(i-1)*Plength+Elength).^2 + lead{ii}(2,(i-1)*Plength+1:(i-1)*Plength+Elength).^2 +lead{ii}(3,(i-1)*Plength+1:(i-1)*Plength+Elength).^2),2);
% %         end
% %     
% %         if ~isempty(find(percentages == i))
% %             disp([num2str(percentages2(percentages == i)) ' % done']);
% %         end
% %         i = i+1;
% %     end
% %     toc;
%     
%     
% % %     
%     lead = source.avg.mom(source.inside);
%     LFMx = [];
%     LFMy = [];
%     LFMz = [];
%     for i = 1:numel(lead)
%         LFMx(end+1,:) = lead{i}(1,:);
%         LFMy(end+1,:) = lead{i}(2,:);
%         LFMz(end+1,:) = lead{i}(3,:);
%     end
%     % 
%     Rdata = abs(hilbert(sqrt(LFMx.^2 + LFMy.^2 + LFMz.^2)));
%     copy = Rdata;
%     c_source = source;
%     disp('Rdata is computed!'); 
%     
%     
% % else
% %     cfg.feedback = 'no';
% %     EpochLength = 1; %epoch length in seconds
% %     Elength = EpochLength*srate;
% %     IElag = 0; % inter-epoch lag (how much does it overlap). e.g. 0.4 is a overlap of 40%!
% %     Plength = srate-(srate*IElag);
% % 
% %     Rdata = zeros(sum(cfg.grid.inside),floor(length(EEGdata.avg)/Plength));
% % 
% %     percentages = floor(length(EEGdata.avg)/Plength)/10:floor(length(EEGdata.avg)/Plength)/10:floor(length(EEGdata.avg)/Plength);
% %     percentages2 = 10:10:100;
% %     disp(['start calculating inverse model for ' num2str(floor(length(EEGdata.avg)/Plength)) ' epochs']);
% %     i = 1;
% %     while i < floor(length(EEGdata.avg)/Plength)
% %         IVdata = EEGdata;
% %         if i == 1
% %             IVdata.avg = IVdata.avg(:,1:Elength);
% %             IVdata.var = IVdata.var(:,1:Elength);
% %         else
% %             IVdata.avg = IVdata.avg(:,(i-1)*Plength+1:(i-1)*Plength+Elength);
% %             IVdata.var = IVdata.var(:,(i-1)*Plength+1:(i-1)*Plength+Elength);
% %         end
% % 
% %         IVdata.time = IVdata.time(1:Elength);
% %         source = ft_sourceanalysis(cfg, IVdata);
% % 
% %         Rdata(:,i) = source.avg.pow(grid.inside)'; % avg.pow is the power for each source!
% %         if ~isempty(find(percentages == i))
% %             disp([num2str(percentages2(percentages == i)) ' % done']);
% %         end
% %         i = i+1;
% %     end
% %     
% %     % normalize Rdata to grand mean across sources?????
% % %     Rdata = (Rdata-mean(mean(Rdata)))./mean(mean(Rdata));
% % end
% 
% 
% %% perform fastICA
% 
% path_fastICA = uigetdir([],'Give me the fastICA folder!');
% path_ICASSO = uigetdir([],'Give me the ICASSO folder!');
% addpath(path_fastICA)
% addpath(path_ICASSO)
% 
% % test = fastica(Rdata, 'interactivePCA','on');
% 
% [ICAcalc] = icassoEst('both', Rdata, 10, 'g', 'tanh', 'approach', 'defl', 'lastEig', 50, 'maxNumIterations',1000); % 'randinit' instead of 'both'?
% [rA] = icassoExp(ICAcalc);
% [Iq, A, W, S] = icassoResult(rA);
% 
% 
% %% correlate each IC with template DMN
% 
% mriTemplate = ft_read_mri('Z:\07_fNetworks_rest-state\data_attentional-load\DMN_template_melodic_IC_sum.nii');
% 
% ft_determine_coordsys(mriTemplate, 'interactive', 'no')
% 
% cfg = [];
% cfg.resolution = 1;
% cfg.dim = [256 256 256];
% mriT = ft_volumereslice(cfg,mriTemplate);
% 
% % cfg = [];
% % cfg.spmversion = 'spm12';
% % cfg.method = 'interactive';
% % cfg.coordsys = 'acpc';
% % cfg.unit = 'mm';
% % [mriT] = ft_volumerealign(cfg,mriTemplate);
% 
% cfg            = [];
% cfg.spmversion = 'spm12';
% % cfg.downsample = 2;
% cfg.parameter = 'anatomy';
% MRIInt  = ft_sourceinterpolate(cfg, mriT, mriT);
% 
% MRIIntNorm = ft_volumenormalise(cfg, MRIInt);
% 
% % plot on brain
% % cfg = [];
% % cfg.spmversion = 'spm12';
% % cfg.method         = 'surface';
% % cfg.funparameter   = 'anatomy';
% % cfg.maskparameter  = cfg.funparameter;
% % cfg.funcolormap    = 'jet'; 
% % cfg.opacitymap     = 'rampup';  
% % cfg.projmethod     = 'nearest'; 
% % cfg.surffile       = 'surface_white_both.mat';
% % cfg.surfdownsample = 10; 
% % cfg.funcolorlim    = [0 200];
% % cfg.opacitylim     = [0 200];
% % ft_sourceplot(cfg, MRIIntNorm);
% % view ([90 0])
% % 
% % cfg              = [];
% % cfg.spmversion = 'spm12';
% % cfg.method       = 'slice';
% % cfg.funparameter = 'anatomy';
% % cfg.maskparameter = cfg.funparameter;
% % cfg.funcolorlim    = [0 200];
% % cfg.opacitylim     = [0 200];
% % cfg.opacitymap    = 'rampup'; 
% % ft_sourceplot(cfg,MRIIntNorm);
% 
% 
% 
% r = [];
% p = [];
% 
% for i = 1:size(A,2)
%     
%     Rdata = abs(A(:,i)); %abs? leave negative values?
%     source.time = 1;
%     source.avg.mom = [];
%     source.avg.ori = [];
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
% 
% %% Interpolate leadfield and plot
% for icI = 28:28
%     Rdata = abs(A(:,icI)); %abs? leave negative values?
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
% 
%     % plot on brain
%     cfg = [];
%     cfg.spmversion = 'spm12';
%     cfg.method         = 'surface';
%     cfg.funparameter   = 'avg.pow';
%     cfg.maskparameter  = cfg.funparameter;
%     % cfg.funcolorlim    = [min(min(min(sourceIntNorm.pow))) max(max(max(sourceIntNorm.pow)))];
%     cfg.funcolormap    = 'jet';
%     % cfg.opacitylim     = [min(min(min(sourceIntNorm.pow))) max(max(max(sourceIntNorm.pow)))]; 
%     cfg.opacitymap     = 'rampup';  
%     cfg.projmethod     = 'nearest'; 
%     cfg.surffile       = 'surface_white_both.mat';
%     cfg.surfdownsample = 10; 
%     ft_sourceplot(cfg, sourceIntNorm);
%     view ([90 0])
% 
%     % plot on MRI
%     cfg              = [];
%     cfg.spmversion = 'spm12';
%     cfg.method       = 'slice';
%     cfg.funparameter = 'pow';
%     cfg.maskparameter = cfg.funparameter;
%     % cfg.funcolorlim   = [min(min(min(sourceIntNorm.pow))) max(max(max(sourceIntNorm.pow)))];
%     % cfg.opacitylim    = [min(min(min(sourceIntNorm.pow))) max(max(max(sourceIntNorm.pow)))]; 
%     cfg.opacitymap    = 'rampup'; 
%     ft_sourceplot(cfg,sourceIntNorm);
% end
% 
% %% plot several figures
% frame = 1000;  %in frames (ms)
% start_frame = 1000;
% end_frame = 2001;
% 
% % frame = 5;  %in frames (ms)
% % start_frame = 540;
% % end_frame = 546;
% 
% % frame = 1;  %in frames (ms)
% % start_frame = 602;
% % end_frame = 603;
% 
% for i = 1:floor((end_frame-start_frame)/frame)
%     Rdata = copy;
%     source = c_source;
%     DURING = mean(Rdata(:,start_frame+(i-1)*frame+1:start_frame+(i)*frame),2);
%     Rdata = DURING;
%     
%     if numel(Rdata) > sum(grid.inside)
%         source.avg.pow(grid.inside) = mean(Rdata,2);
%     else
%         source.avg.pow(source.inside) = Rdata;
%     end
% 
%     source.time = 1;
%     rmfield(source.avg,'mom');
%     rmfield(source.avg,'ori');
%     
%     cfg            = [];
% %     cfg.downsample = 2;
%     cfg.parameter = 'pow';
%     cfg.spmversion = 'spm12';
%     sourceInt  = ft_sourceinterpolate(cfg, source , mri);
% 
%     % normalize to template
%     sourceIntNorm = ft_volumenormalise(cfg, sourceInt);
% 
%     % plot on brain
%     cfg = [];
%     cfg.spmversion = 'spm12';
%     cfg.method         = 'surface';
%     cfg.funparameter   = 'avg.pow';
%     cfg.maskparameter  = cfg.funparameter;
% %     cfg.funcolorlim    = [0 1];
%     cfg.funcolormap    = 'jet';
% %     cfg.opacitylim     = [0 1]; 
%     cfg.opacitymap     = 'rampup';  
%     cfg.projmethod     = 'nearest'; 
%     cfg.surffile       = 'surface_white_both.mat';
%     cfg.surfdownsample = 10; 
%     ft_sourceplot(cfg, sourceIntNorm);
%     view ([90 0])
%     
%     % plot on MRI
%     cfg              = [];
%     cfg.spmversion = 'spm12';
%     cfg.method       = 'slice';
%     cfg.funparameter = 'pow';
%     cfg.maskparameter = cfg.funparameter;
% %     cfg.funcolorlim   = [0 1];
% %     cfg.opacitylim    = [0 1]; 
%     cfg.opacitymap    = 'rampup'; 
%     ft_sourceplot(cfg,sourceIntNorm);
% end
% 
% 
% %% write to nifti for MRIcron plotting
% % save_folder = uigetdir([],'where do you want to save nifti of sources?');
% 
% cfg = [];
% cfg.coordsys = 'mni';
% sourceMNI = ft_volumenormalise(cfg, sourceIntNorm);
% 
% cfg = [];
% cfg.filename  = [save_folder '\source_activity'];
% cfg.filetype  = 'nifti';
% cfg.parameter = 'pow';
% cfg.precision = 'single';
% ft_sourcewrite(cfg, sourceMNI);