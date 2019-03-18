% create the MIDA model with 12 tissue types

%change path to MIDA-surface files
cd('Z:\07_fNetworks_rest-state\MIDA_head-model\MIDAv1.0\MIDA_v1.0\MIDA_v1_surfaces'); 

Tissues = {'Skin','Eyes','Muscle','Fat','SpongyBone','CompactBone','Gray','CerebellarGray','White','CerebellarWhite','CSF','BrainStem'};

%Classify surface files to tissues
Skin = {'Epidermis_Dermis.stl'};
Eyes = {'Eye Aqueous.stl','Eye Cornea.stl','Eye Lens.stl','Eye Retina_Choroid_Sclera.stl','Eye Vitreous .stl'};
Muscle = {'Muscle (General).stl'};
Fat = {'Adipose Tissue.stl','Subcutaneous Adipose Tissue.stl'};
SpongyBone = {'Skull Diploe.stl'};
CompactBone = {'Skull Inner Table.stl','Skull Outer Table .stl'};
Gray = {'Brain Gray Matter.stl'};
CerebellarGray = {'Cerebellum Gray Matter.stl'};
White = {'Brain White Matter.stl'};
CerebellarWhite = {'Cerebellum White Matter.stl'};
CSF = {'CSF General.stl'};
BrainStem = {'Brainstem Medulla.stl','Brainstem Midbrain.stl','Brainstem Pons.stl'};

for i = 1:length(Tissues)
    Segment = ft_read_headshape(eval(Tissues{i}),'concatenate','no');
    assignin('base', Tissues{i}, Segment);
end

% concatinate meshes Eyes
newEyes = struct;

newEyes.pos = vertcat(Eyes(1).pos, Eyes(2).pos);
faceUpdate = Eyes(2).tri+max(max(Eyes(1).tri));
newEyes.tri = vertcat(Eyes(1).tri,faceUpdate);

newEyes.pos = vertcat(newEyes.pos, Eyes(3).pos);
faceUpdate = Eyes(3).tri+max(max(newEyes.tri));
newEyes.tri = vertcat(newEyes.tri,faceUpdate);

newEyes.pos = vertcat(newEyes.pos, Eyes(4).pos);
faceUpdate = Eyes(4).tri+max(max(newEyes.tri));
newEyes.tri = vertcat(newEyes.tri,faceUpdate);

newEyes.pos = vertcat(newEyes.pos, Eyes(5).pos);
faceUpdate = Eyes(5).tri+max(max(newEyes.tri));
newEyes.tri = vertcat(newEyes.tri,faceUpdate);

Eyes = newEyes;
Eyes.unit = 'mm';
clear newEyes

% concatinate meshes Fat
newFat = struct;

newFat.pos = vertcat(Fat(1).pos, Fat(2).pos);
faceUpdate = Fat(2).tri+max(max(Fat(1).tri));
newFat.tri = vertcat(Fat(1).tri,faceUpdate);

Fat = newFat;
Fat.unit = 'mm';
clear newFat

% concatinate meshes compactBone
newBone = struct;

newBone.pos = vertcat(CompactBone(1).pos, CompactBone(2).pos);
faceUpdate = CompactBone(2).tri+max(max(CompactBone(1).tri));
newBone.tri = vertcat(CompactBone(1).tri,faceUpdate);

CompactBone = newBone;
CompactBone.unit = 'mm';
clear newBone

% concatinate meshes Brainstem
newStem = struct;

newStem.pos = vertcat(BrainStem(1).pos, BrainStem(2).pos);
faceUpdate = BrainStem(2).tri+max(max(BrainStem(1).tri));
newStem.tri = vertcat(BrainStem(1).tri,faceUpdate);

newStem.pos = vertcat(newStem.pos, BrainStem(3).pos);
faceUpdate = BrainStem(3).tri+max(max(newStem.tri));
newStem.tri = vertcat(newStem.tri,faceUpdate);

BrainStem = newStem;
BrainStem.unit = 'mm';
clear newStem


%% downsample
addpath([path_ft '\external\iso2mesh']);
[Skin.pos, Skin.tri] = meshresample(Skin.pos, Skin.tri,1000/size(Skin.pos,1));
[CSF.pos, CSF.tri] = meshresample(CSF.pos, CSF.tri,100000/size(CSF.pos,1));

Seg.dim = size(MIDA.anatomy);
% Seg.coordsys = 'ctf';
Seg.transform = MIDA.transform;
Seg.unit = 'mm';

% 
% 
% %plot meshes
% colorV = rand(length(Tissues),3);
% 
% figure;
% hold on
% for i = 1:length(Tissues)
%     ft_plot_mesh(eval(Tissues{i}),'edgecolor','none', 'facecolor',colorV(i,:), 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
% end


% headmodel
cfg        = [];
cfg.method ='openmeeg';
cfg.conductivity = [0.33];   % order follows mesh.tissyelabel
vol        = ft_prepare_headmodel(cfg, Seg);    