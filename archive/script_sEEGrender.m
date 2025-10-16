
%% example script to render sEEG electrode positions 

% Dora Hermes, Mayo Clinic 2025

%% Set data path

localDataPath = 'where I put the iEEG BIDS data';


subjects = {'01'};

ss = 1;
sub_label = subjects{ss};

% get elecmatrix
electrodes_tsv_name = fullfile(localDataPath,['sub-' sub_label],'ses-ieeg01','ieeg',...
    ['sub-' sub_label '_ses-ieeg01_electrodes.tsv']);
loc_info = readtable(electrodes_tsv_name, 'FileType', 'text', 'Delimiter', '\t');

% load surface should be the same hemi/space as electrodes
hemi = 'r';
g = gifti(fullfile(localDataPath,'derivatives','freesurfer',['sub-' sub_label],['pial.' upper(hemi) '.surf.gii']));
            
figure, hold on

% render
tH = ieeg_RenderGifti(g); 
tH.FaceAlpha = 0.5;

% add electrode positions
xyz = [loc_info.x loc_info.y loc_info.z];
s = scatter3(xyz(:,1),xyz(:,2),xyz(:,3),'o','Filled');
s.SizeData = 10;
s.MarkerFaceColor = [.9 .5 .5];
s.MarkerEdgeColor = [0 0 0];

% adjust light
ieeg_viewLight(90,0)



