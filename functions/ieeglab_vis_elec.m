function ieeglab_vis_elec(EEG)

% Requirements
%   - Vistasoft --> should be in iEEGLAB > external (as submodule)
% addpath(genpath(fullfile(plugin_path, 'external')));

figure('color','w')
hold on
try icadefs; set(gcf, 'color', BACKCOLOR); catch; end  % eeglab color

sub_files = { dir(fullfile(EEG.filepath)).name }';
surf_files = sub_files(contains(sub_files, 'pial'));
if isempty(surf_files)
    error("Sorry, no pial surface file in this subject's folder. Aborting visualization. If this file does not exist, you need to generate it with Freesurfer using the subject's MRI.")
end

% Left hemisphere
% idx = strcmpi({EEG.chanlocs.hemisphere}, 'L');
idx = [EEG.chanlocs.X] < 0;
if sum(idx) > 1
    surfFile = surf_files{contains(surf_files, 'L')};
    % g = gifti(fullfile(EEG.filepath,'pial.L.surf.gii' ));
    % g = gifti(fullfile(EEG.filepath,'pial_desc-qsiprep.L.surf.gii' ));    
    g = gifti(fullfile(EEG.filepath, surfFile));
    tH = ieeg_RenderGifti(g);  % render
    tH.FaceAlpha = 0.2;
    s = scatter3([EEG.chanlocs(idx).X], [EEG.chanlocs(idx).Y], [EEG.chanlocs(idx).Z], 'o', 'Filled');
    s.SizeData = 10;
    s.MarkerFaceColor = [.9 .5 .5];
    s.MarkerEdgeColor = [0 0 0];
    ieeg_viewLight(90,0)
else
    disp("No electrodes in left hemisphere to plot")
end

% Right hemisphere
% idx = strcmpi({EEG.chanlocs.hemisphere}, 'R');
idx = [EEG.chanlocs.X] > 0;
if sum(idx) > 1
    surfFile = surf_files{contains(surf_files, 'R')};
    g = gifti(fullfile(EEG.filepath, surfFile));

    tH = ieeg_RenderGifti(g);  % render
    tH.FaceAlpha = 0.2;
    s = scatter3([EEG.chanlocs(idx).X], [EEG.chanlocs(idx).Y], [EEG.chanlocs(idx).Z], 'o', 'Filled');
    s.SizeData = 10;
    s.MarkerFaceColor = [.9 .5 .5];
    s.MarkerEdgeColor = [0 0 0];
    ieeg_viewLight(90,0)  % 180,90 (from the top with the front downward);  90,0 (from the right)
else
    disp("No electrodes in right hemisphere to plot")
end
