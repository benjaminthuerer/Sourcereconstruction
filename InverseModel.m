%% Constants and options
path_ft = uigetdir([],'Give me the Field Trip folder!');
[LF_Head_path] = uigetdir([],'Feed the leadfield and headmodel folder');
[EEG_file,EEG_path] = uigetfile('*.*','Feed me the EEG data');

LF_data = [LF_Head_path,'\leadfield.mat'];
Head_data = [LF_Head_path,'\headmodel.mat'];
MRI_data = [LF_Head_path,'\mri.mat'];

% Add path to field trip and set general settings


tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

addpath('Z:\Matlab_Scripts\eeglab13_5_4b\');

eeglab
close

restoredefaultpath
addpath(path_ft)
ft_defaults

%% load leadfield and EEG data
disp('Load MRI, leadfield, and headmodel. This may take a while');
load(LF_data,'grid');
load(Head_data,'vol');
load(MRI_data,'mri');
disp('Leadfield and headmodel ready to use');

EEG = pop_loadset('filename',EEG_file,'filepath',EEG_path);
EEG = eeg_checkset( EEG );

srate = EEG.srate;

    
label_cell = strfind({EEG.chanlocs.labels},'EOG');
for i = 1:length(label_cell)
    if ~isempty(label_cell{i} == 1)
        label_idx(i) = 1;
    else
        label_idx(i) = 0;
    end
end
label_idx = label_idx == 1;
    
EEG.data(label_idx,:) = [];
EEG.chanlocs(label_idx) = [];
EEG.nbchan = size(EEG.data,1);



%% transform data for fieldtrip processing

% just a hack for eeglab2fieldtrip
EEG.icachansind = randi(1,1,size(EEG.data,1));

fieldbox = 'timelockanalysis';
transform = 'none'; %or DIPTFIT transformation of channel locations
EEGdata = eeglab2fieldtrip(EEG, fieldbox,transform);

% EEGdata.elec = grid.cfg.elec; %elec positions transformed from leadfield??????

if ~isfield(EEGdata,'trial')
    EEGdata.dimord = 'chan_time';
else
    EEGdata.dimord = 'rpt_chan_time';
end

EEGdata = rmfield(EEGdata,'fsample');
EEGdata.label = EEGdata.label';

rmpath(genpath('Z:\Matlab_Scripts\eeglab13_5_4b\'));

%% if continuous data: loop in 1s epochs through inverse model using eLORETA

cfg = [];
cfg.method = 'eloreta';
% cfg.lambda = 3;
cfg.grid = grid;
cfg.headmodel = vol;
% cfg.eloreta.keepfilter              = 'yes';
% cfg.eloreta.normalize               = 'yes';
% cfg.eloreta.lambda                  = 0.05;
% cfg.eloreta.projectnoise            = 'yes';


if ndims(EEG.data)==3
    cfg.feedback = 'yes';
    Rdata = [];
    source = ft_sourceanalysis(cfg, EEGdata);
else
    cfg.feedback = 'no';
    EpochLength = 1; %epoch length in seconds
    Elength = EpochLength*srate;
    IElag = 0.5; % inter-epoch lag (how much does it overlap). e.g. 0.4 is a overlap of 40%!
    Plength = srate-(srate*IElag);

    Rdata = zeros(sum(cfg.grid.inside),floor(length(EEGdata.avg)/(Elength-Plength)));

    percentages = floor(length(EEGdata.avg)/(Elength-Plength))/10:floor(length(EEGdata.avg)/(Elength-Plength))/10:floor(length(EEGdata.avg)/(Elength-Plength));
    percentages2 = 10:10:100;
    disp(['start calculating inverse model for ' num2str(floor(length(EEGdata.avg)/(Elength-Plength))) ' epochs']);
    i = 1;
    while i < floor(length(EEGdata.avg)/(Elength-Plength))
        IVdata = EEGdata;
        if i == 1
            IVdata.avg = IVdata.avg(:,1:Elength);
            IVdata.var = IVdata.var(:,1:Elength);
        else
            IVdata.avg = IVdata.avg(:,(i-1)*(Elength-Plength)+1:(i-1)*(Elength-Plength)+Elength);
            IVdata.var = IVdata.var(:,(i-1)*(Elength-Plength)+1:(i-1)*(Elength-Plength)+Elength);
        end

        IVdata.time = IVdata.time(1:Elength);
        source = ft_sourceanalysis(cfg, IVdata);

        Rdata(:,i) = source.avg.pow(grid.inside)'; % avg.pow is the power for each source!
        if ~isempty(find(percentages == i))
            disp([num2str(percentages2(percentages == i)) ' % done']);
        end
        i = i+1;
    end
    
    % normalize Rdata to grand mean across sources?????
    Rdata = (Rdata-mean(mean(Rdata)))./mean(mean(Rdata));
end


lead = source.avg.mom(source.inside);
LFMx = [];
LFMy = [];
LFMz = [];
for i = 1:numel(lead)
    LFMx(end+1,:) = lead{i}(1,:);
    LFMy(end+1,:) = lead{i}(2,:);
    LFMz(end+1,:) = lead{i}(3,:);
end
% 
Rdata = sqrt(LFMx.^2 + LFMy.^2 + LFMz.^2);
copy = Rdata;
c_source = source;

disp('Rdata is computed!');

%%
Rdata = Rdata_pos + (Rdata_neg.*-1);

Rdata = (Rdata-mean(Rdata(:,1:500),2))./mean(Rdata(:,1:500),2); %normalize each source to a baseline
% Rdata = (Rdata-mean(Rdata(:,:),1))./mean(Rdata(:,:),1); %normalize each time bin to the mean source activity
copy = Rdata;
c_source = source;



% 
% % PRE = mean(Rdata(:,1:500),2);
% 
% DURING = mean(Rdata(:,800:820),2);
% % 
% Rdata = DURING;
%% perform fastICA

path_fastICA = uigetdir([],'Give me the fastICA folder!');
addpath(path_fastICA)

icaResult = fastica(copy, 'numOfIC', 50); %numOfIC change to individual threshold using PCA or minimum discription lenght (MDL) criterion


%% Interpolate leadfield and plot

if numel(Rdata) > sum(grid.inside)
    source.avg.pow(grid.inside) = mean(Rdata,2);
else
    source.avg.pow(source.inside) = Rdata;
end

cfg            = [];
cfg.downsample = 2;
cfg.parameter = 'pow';
sourceInt  = ft_sourceinterpolate(cfg, source , mri);

% normalize to template
sourceIntNorm = ft_volumenormalise(cfg, sourceInt);

% plot on brain
cfg = [];
cfg.method         = 'surface';
cfg.funparameter   = 'avg.pow';
cfg.maskparameter  = cfg.funparameter;
cfg.funcolorlim    = [min(min(min(sourceIntNorm.pow))) max(max(max(sourceIntNorm.pow)))];
cfg.funcolormap    = 'jet';
cfg.opacitylim     = [min(min(min(sourceIntNorm.pow))) max(max(max(sourceIntNorm.pow)))]; 
cfg.opacitymap     = 'rampup';  
cfg.projmethod     = 'nearest'; 
cfg.surffile       = 'surface_white_both.mat';
cfg.surfdownsample = 10; 
ft_sourceplot(cfg, sourceIntNorm);
view ([90 0])

% plot on MRI
cfg              = [];
cfg.method       = 'slice';
cfg.funparameter = 'pow';
cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim   = [min(min(min(sourceIntNorm.pow))) max(max(max(sourceIntNorm.pow)))];
cfg.opacitylim    = [min(min(min(sourceIntNorm.pow))) max(max(max(sourceIntNorm.pow)))]; 
cfg.opacitymap    = 'rampup'; 
ft_sourceplot(cfg,sourceIntNorm);


%% plot several figures
% frame = 100;  %in frames (ms)
% start_frame = 1000;
% end_frame = 1101;

frame = 5;  %in frames (ms)
start_frame = 540;
end_frame = 546;

% frame = 1;  %in frames (ms)
% start_frame = 602;
% end_frame = 603;

for i = 1:floor((end_frame-start_frame)/frame)
    Rdata = copy;
    source = c_source;
    DURING = mean(Rdata(:,start_frame+(i-1)*frame+1:start_frame+(i)*frame),2);
    Rdata = DURING;
    
    if numel(Rdata) > sum(grid.inside)
        source.avg.pow(grid.inside) = mean(Rdata,2);
    else
        source.avg.pow(source.inside) = Rdata;
    end

    cfg            = [];
%     cfg.downsample = 2;
    cfg.parameter = 'pow';
    sourceInt  = ft_sourceinterpolate(cfg, source , mri);

    % normalize to template
    sourceIntNorm = ft_volumenormalise(cfg, sourceInt);

    % plot on brain
    cfg = [];
    cfg.method         = 'surface';
    cfg.funparameter   = 'avg.pow';
    cfg.maskparameter  = cfg.funparameter;
%     cfg.funcolorlim    = [0 1];
    cfg.funcolormap    = 'jet';
%     cfg.opacitylim     = [0 1]; 
    cfg.opacitymap     = 'rampup';  
    cfg.projmethod     = 'nearest'; 
    cfg.surffile       = 'surface_white_both.mat';
    cfg.surfdownsample = 10; 
    ft_sourceplot(cfg, sourceIntNorm);
    view ([90 0])
    
    % plot on MRI
    cfg              = [];
    cfg.method       = 'slice';
    cfg.funparameter = 'pow';
    cfg.maskparameter = cfg.funparameter;
%     cfg.funcolorlim   = [0 1];
%     cfg.opacitylim    = [0 1]; 
    cfg.opacitymap    = 'rampup'; 
    ft_sourceplot(cfg,sourceIntNorm);
end