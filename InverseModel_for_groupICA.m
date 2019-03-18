path_ft = uigetdir([],'Give me the Field Trip folder!');
path_fastICA = uigetdir([],'Give me the fastICA folder!');
path_ICASSO = uigetdir([],'Give me the ICASSO folder!');

% tmp = matlab.desktop.editor.getActive;
% cd(fileparts(tmp.Filename));

restoredefaultpath
addpath(path_ft)
ft_defaults

addpath('Z:\Matlab_Scripts\Fieldtrip\new_fieldtrip\external\eeglab')
addpath(genpath('Z:\Matlab_Scripts\GroupICATv4.0b\GroupICATv4.0b\icatb'))

subjIdx = {'1001', '1002', '1003', '1004', '1006', '1007', '1008', '1011'};

Data_path_general = 'Z:\07_fNetworks_rest-state\data_attentional-load\Processed\';
MRIdata_path_general = 'Z:\07_fNetworks_rest-state\data_attentional-load\MR\';

subjectDMN = [];
EstComp = [];
MaxCorr = [];
IdxCorr = [];

for i = 1:length(subjIdx)
    % load MRI
    MRI_data = 'IM_0001.nii';
    MRI_path = [MRIdata_path_general subjIdx{i} '\'];
    mri = ft_read_mri([MRI_path MRI_data]);
    
    % bring MRI from nifti to acpc. needed?
    cfg = [];
    cfg.spmversion = 'spm12';
    cfg.method = 'interactive';
    cfg.coordsys = 'acpc';
    cfg.unit = 'mm';
    [mri] = ft_volumerealign(cfg,mri);
%     mri = ft_volumereslice([], mri);

    % set file paths
    LF_Head_path = [Data_path_general subjIdx{i}];
    LF_data = [LF_Head_path,'\leadfield_12T_FEM_gray-only.mat'];
    Head_data = [LF_Head_path,'\headmodel_12T_FEM_prepared_sens_vol.mat'];
    EEG_path = [Data_path_general 's' subjIdx{i} '\auto_Rest\'];
    EEG_file = ['s' subjIdx{i} '_auto_Rest_ica_cleaned'];
    
    % load liedfield, headmodel, and eeg
    disp('Load MRI, leadfield, and headmodel. This may take a while');
    load(LF_data,'grid');
    disp('leadfield loaded')
    load(Head_data,'vol');
    disp('headmodel loaded')
    load([EEG_path EEG_file]);
    disp('EEG data loaded')
    
    % transform EEG to fieldtrip structure
    srate = EEG.srate;
    EEG.icachansind = 1:size(EEG.data,1);
    fieldbox = 'timelockanalysis';
    transform = 'none'; %or DIPTFIT transformation of channel locations
    EEGdata = eeglab2fieldtrip(EEG, fieldbox,transform);
    EEGdata.dimord = 'chan_time';
    EEGdata = rmfield(EEGdata,'fsample');
    
    EEGdata.avg = double(EEGdata.avg);
    EEGdata.var = double(EEGdata.var);
    
    % perform inverse model using eloreta
    cfg = [];
    cfg.method = 'eloreta';
    cfg.grid = grid;
    cfg.eloreta.keepfilter              = 'yes';
    cfg.eloreta.normalize               = 'yes';
    cfg.eloreta.lambda                  = 0.05;
    cfg.eloreta.projectnoise            = 'yes';
    cfg.headmodel = vol;
    source = ft_sourceanalysis(cfg, EEGdata);
    
    % perform source projection and downsample with anti-aliasing to 1s
    % epochs
    lead = source.avg.mom(source.inside);

    EpochLength = 1; %epoch length in seconds
    Elength = EpochLength*srate;
    IElag = 0; % inter-epoch lag (how much does it overlap). e.g. 0.4 is a overlap of 40%!
    Plength = Elength-(Elength*IElag);
    Rdata = zeros(sum(grid.inside),ceil(length(EEGdata.avg)/Plength));
    
    percentages = floor(length(lead)/10:length(lead)/10:length(lead));
    percentages2 = 10:10:100;
    disp(['start computing the source projection for ' num2str(sum(grid.inside)) ' sources. This may take a while...']);
    
    for ii = 1:sum(grid.inside)
        Stime = abs(hilbert(sqrt(double(lead{ii}(1,:)).^2 + double(lead{ii}(2,:)).^2 + double(lead{ii}(3,:)).^2)));
        Rdata(ii,:) = resample(Stime,1,srate);
        if ~isempty(find(percentages == ii))
            disp([num2str(percentages2(percentages == ii)) ' % done']);
        end
    end
 
    %some cleaning to save workspace
    clear lead EEG EEGdata Stime grid vol cfg
      
    
    % compute and save brain map for each time bin of the source-space
    source.time = 1;
    source.avg.mom = [];
    source.avg.ori = [];
    source.avg.pow(source.inside) = Rdata(:,1);
    
    cfg            = [];
    cfg.spmversion = 'spm12';
    cfg.parameter = 'avg.pow';
    cfg.downsample = 1;
    sourceInt  = ft_sourceinterpolate(cfg, source , mri);

    cfg            = [];
    cfg.nonlinear = 'no';
    cfg.spmversion = 'spm12';
    cfg.parameter = 'pow';
    cfg.downsample = 1;
    sourceIntNorm = ft_volumenormalise(cfg, sourceInt);
        
    files = [];
    mkdir([Data_path_general 'groupICA\' subjIdx{i} '\'])
    
    
    for ii = 1:size(Rdata,2)
        source.time = 1;
        source.avg.mom = [];
        source.avg.ori = [];
        source.avg.pow(source.inside) = Rdata(:,ii);
        
        cfg            = [];
        cfg.spmversion = 'spm12';
        cfg.parameter = 'avg.pow';
        cfg.downsample = 1;
        sourceInt  = ft_sourceinterpolate(cfg, source , mri);

        cfg            = [];
        cfg.nonlinear = 'no';
        cfg.spmversion = 'spm12';
        cfg.parameter = 'pow';
        cfg.downsample = 1;
        sourceIntNorm = ft_volumenormalise(cfg, sourceInt);
        
        cfg = [];
        cfg.filename  = [Data_path_general 'groupICA\' subjIdx{i} '\' subjIdx{i} '_' sprintf( '%06d', ii)];
        cfg.spmversion = 'spm12';
        cfg.filetype  = 'nifti';
        cfg.parameter = 'pow';
        cfg.precision = 'single';
        ft_sourcewrite(cfg, sourceIntNorm);
        
        files = [files; [Data_path_general 'groupICA\' subjIdx{i} '\' subjIdx{i} '_' sprintf( '%06d', ii) '.nii']];
    end
      
    
    %% from here the code is crap (icasso is not really working on reshaped 3d files but MDL is not computed in source space (nifti)
    %best result when performing group ICA in GIFT. Problem: loss of time
    %domain and only template can be used...
    
    % estimate number of ICs using MDL approach
    disp('start to estimate number of components by MDL. This may take a while...');
    [comp_est, mdl, aic, kic] = icatb_estimate_dimension(files);
    
    EstComp(i) = comp_est;
    
    addpath(path_fastICA)
    addpath(path_ICASSO)
    
    %perform sICA on single subject data
    [ICAcalc] = icassoEst('bootstrap', Rdata', 10, 'g', 'tanh', 'approach', 'defl', 'numOfIC', comp_est, 'maxNumIterations', 1000); % reconstruct first before ICA?
    sR = icassoCluster(ICAcalc, 'strategy', 'AL', 'simfcn', 'abscorr', 's2d', 'sim2dis', 'L', comp_est); 
    sR = icassoProjection(sR, 'cca', 's2d', 'sqrtsim2dis'); 
    [Iq2, A2, W2, S2] = icassoResult(sR, comp_est); 
    
    mriTemplate = ft_read_mri('Z:\07_fNetworks_rest-state\data_attentional-load\DMN_template_melodic_IC_sum.nii');
    
    % Template seems to be already in MNI space
    ft_determine_coordsys(mriTemplate)
    
    cfg = [];
    cfg.resolution = 1;
    cfg.dim = [256 256 184]; %is this still correct?
    mriT = ft_volumereslice(cfg,mriTemplate);
    

    % correlate each IC with the template and find the best fit (Perason's r)    
    copyR = Rdata;
    
    for ii = 1:size(S2,1)

        cfg            = [];
        cfg.spmversion = 'spm12';
        cfg.parameter = 'avg.pow';
        sourceInt  = ft_sourceinterpolate(cfg, source , mri);

        % normalize to template
        cfg            = [];
        cfg.nonlinear = 'no';
        cfg.spmversion = 'spm12';
        cfg.parameter = 'pow';
        sourceIntNorm = ft_volumenormalise(cfg, sourceInt);

        cfg              = [];
        cfg.method       = 'slice';
        cfg.funparameter = 'pow';
        ft_sourceplot(cfg,sourceIntNorm);
        
        cfg = [];
        cfg.spmversion = 'spm12';
        cfg.method         = 'surface';
        cfg.funparameter   = 'pow';
        cfg.maskparameter  = cfg.funparameter;
        cfg.funcolormap    = 'jet';
        cfg.opacitymap     = 'rampup';  
        cfg.projmethod     = 'nearest'; 
        cfg.surffile       = 'surface_white_both.mat';
        cfg.surfdownsample = 10; 
        ft_sourceplot(cfg, sourceIntNorm);
        view ([90 0])

        [r(ii),p(ii)] = corr(reshape(mriT.anatomy,numel(mriT.inside),1),reshape(sourceIntNorm.pow,numel(sourceIntNorm.pow),1), 'rows','complete');
    end

    [Vm,~,Im] = max(abs(r));
    
    MaxCorr(i) = Vm;
    IdxCorr(i) = Im;
    
    % plot and save the best fit (individual subjet's DMN)
    Rdata = S(Im,:);
    source.time = 1;
    source.avg.mom = [];
    source.avg.ori = [];

    source.avg.pow(source.inside) = Rdata';

    cfg            = [];
    cfg.spmversion = 'spm12';
    cfg.parameter = 'avg.pow';
    sourceInt  = ft_sourceinterpolate(cfg, source , mri);

    % normalize to template
        cfg            = [];
    cfg.nonlinear = 'no';
    cfg.spmversion = 'spm12';
    cfg.parameter = 'pow';
    sourceIntNorm = ft_volumenormalise(cfg, sourceInt);

    cfg              = [];
    cfg.spmversion = 'spm12';
    cfg.method       = 'slice';
    cfg.funparameter = 'pow';
    cfg.maskparameter = cfg.funparameter;
    cfg.opacitymap    = 'rampup'; 
    ft_sourceplot(cfg,sourceIntNorm);
    
    cfg = [];
    cfg.filename  = [Data_path_general 'reconstr_DMN\DMN_' subjIdx{i}];
    cfg.filetype  = 'nifti';
    cfg.parameter = 'pow';
    cfg.precision = 'single';
    ft_sourcewrite(cfg, sourceIntNorm);
    
    subjectDMN(:,:,:,i) = sourceIntNorm.pow;
    dummySource = sourceIntNorm;
end

% mean over subjects gives groupICA DMN
MsubjectDMN = mean(subjectDMN,4);

dummySource.pow = MsubjectDMN;

% plot and save groupICA DMN
cfg              = [];
cfg.spmversion = 'spm12';
cfg.method       = 'slice';
cfg.funparameter = 'anatomy';
cfg.maskparameter = cfg.funparameter;
cfg.opacitymap    = 'rampup'; 
ft_sourceplot(cfg,dummyDource);

cfg = [];
cfg.filename  = [Data_path_general 'reconstr_DMN\groupICA_DMN'];
cfg.filetype  = 'nifti';
cfg.parameter = 'pow';
cfg.precision = 'single';
ft_sourcewrite(cfg, sourceIntNorm);

logFile = struct;
logFile.ICest = EstComp;
logFile.MaxCorr = MaxCorr;
logFile.IdxCorr = IdxCorr;

save([Data_path_general 'OutputFile.mat'],'logFile','-v7.3');



%% test for validity
%- check if topoplots in EEG and source reconstructed look the same
%- check if topoplots of EEG-ICs and source ICs look the same

%- if yes: source reconstruction, ICA and source reconstruction works!
%   --> go with groupICA using GIFT

%- topoplots, timeline and plotting looks good!
%- ICA topoplots look similar but only 4 out of 5 had a counterpart
%
%- restriction to gray matter only (not basal ganglia etc.). this gives
%less dipoles and better reconstruction?


source.time = 1;
source.avg.mom = [];
source.avg.ori = [];

for k = 1:5
    source.avg.pow(source.inside) = S2(k,:);

    cfg            = [];
    cfg.spmversion = 'spm12';
    cfg.parameter = 'avg.pow';
    sourceInt  = ft_sourceinterpolate(cfg, source , mri);

    % normalize to template
    cfg            = [];
    cfg.nonlinear = 'no';
    cfg.spmversion = 'spm12';
    cfg.parameter = 'pow';
    sourceIntNorm = ft_volumenormalise(cfg, sourceInt);

    cfg              = [];
    cfg.spmversion = 'spm12';
    cfg.method       = 'slice';
    cfg.funparameter = 'pow';
    % cfg.funcolormap    = 'jet';
    cfg.opacitymap     = 'rampup';  
    cfg.projmethod     = 'nearest';
    ft_sourceplot(cfg,sourceIntNorm);

    cfg = [];
    cfg.spmversion = 'spm12';
    cfg.method         = 'surface';
    cfg.funparameter   = 'pow';
    cfg.maskparameter  = cfg.funparameter;
    cfg.funcolormap    = 'jet';
    cfg.opacitymap     = 'rampup';  
    cfg.projmethod     = 'nearest'; 
    cfg.surffile       = 'surface_white_both.mat';
    cfg.surfdownsample = 10; 
    ft_sourceplot(cfg, sourceIntNorm);
    view ([90 0])
end

