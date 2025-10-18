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
idx_elecs_hemi = [EEG.chanlocs.X] < 0;
idx_surf_hemi = contains(surf_files, 'L');
if sum(idx_elecs_hemi) > 1 && any(idx_surf_hemi)
    surfFile = surf_files{idx_surf_hemi};
    % g = gifti(fullfile(EEG.filepath,'pial.L.surf.gii' ));
    % g = gifti(fullfile(EEG.filepath,'pial_desc-qsiprep.L.surf.gii' ));    
    g = gifti(fullfile(EEG.filepath, surfFile));
    tH = ieeg_RenderGifti(g);  % render
    tH.FaceAlpha = 0.2;
    s = scatter3([EEG.chanlocs(idx_elecs_hemi).X], [EEG.chanlocs(idx_elecs_hemi).Y], [EEG.chanlocs(idx_elecs_hemi).Z], 'o', 'Filled');
    s.SizeData = 10;
    s.MarkerFaceColor = [.9 .5 .5];
    s.MarkerEdgeColor = [0 0 0];
    ieeg_viewLight(90,0)
else
    disp("No electrodes in left hemisphere to plot")
end

% Right hemisphere
% idx = strcmpi({EEG.chanlocs.hemisphere}, 'R');
idx_elecs_hemi = [EEG.chanlocs.X] > 0;
if sum(idx_elecs_hemi) > 1 && any(idx_surf_hemi)
    surfFile = surf_files{idx_surf_hemi};
    g = gifti(fullfile(EEG.filepath, surfFile));

    tH = ieeg_RenderGifti(g);  % render
    tH.FaceAlpha = 0.2;
    s = scatter3([EEG.chanlocs(idx_elecs_hemi).X], [EEG.chanlocs(idx_elecs_hemi).Y], [EEG.chanlocs(idx_elecs_hemi).Z], 'o', 'Filled');
    s.SizeData = 10;
    s.MarkerFaceColor = [.9 .5 .5];
    s.MarkerEdgeColor = [0 0 0];
    ieeg_viewLight(90,0)  % 180,90 (from the top with the front downward);  90,0 (from the right)
else
    disp("No electrodes in right hemisphere to plot")
end
