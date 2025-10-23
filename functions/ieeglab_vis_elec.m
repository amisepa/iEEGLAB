function ieeglab_vis_elec(EEG)
% ieeglab_vis_elec
% - If subject pial *.gii surfaces exist in EEG.filepath, uses your original code (unchanged).
% - Else: uses high-res cortex from 'cortex.mat' (variable cortex_highres),
%         auto-aligns electrodes (orientation + rigid PCA + uniform scale),
%         and plots a single rotatable 3D.
% - If cortex.mat is missing, falls back to dipfit/standard_BEM (smooth).

figure('color','w'); hold on
try icadefs; set(gcf, 'color', BACKCOLOR); catch, end  % eeglab color

sub_files  = { dir(fullfile(EEG.filepath)).name }';
surf_files = sub_files(contains(sub_files, 'pial')); % 'pial' or 'white' work

% Subject pial-surface plotting from vistasoft
if ~isempty(surf_files)        
        
    try icadefs; set(gcf, 'color', BACKCOLOR); catch; end  % eeglab color

    % Left hemisphere
    idx_elecs_hemi = [EEG.chanlocs.X] < 0;
    idx_surf_hemi  = contains(surf_files, 'L');
    if sum(idx_elecs_hemi) > 1 && any(idx_surf_hemi)
        surfFile = surf_files{idx_surf_hemi};
        g = gifti(fullfile(EEG.filepath, surfFile));
        tH = ieeg_RenderGifti(g);  % your renderer
        tH.FaceAlpha = 0.1;
        s = scatter3([EEG.chanlocs(idx_elecs_hemi).X], [EEG.chanlocs(idx_elecs_hemi).Y], [EEG.chanlocs(idx_elecs_hemi).Z], 'o', 'Filled');
        s.SizeData = 10; 
        s.MarkerFaceColor = [.9 .5 .5]; 
        s.MarkerEdgeColor = [0 0 0];
        ieeg_viewLight(90,0)
    else
        disp("No electrodes in left hemisphere to plot")
    end

    % Right hemisphere
    idx_elecs_hemi = [EEG.chanlocs.X] > 0;
    idx_surf_hemi  = contains(surf_files, 'R');
    if sum(idx_elecs_hemi) > 1 && any(idx_surf_hemi)
        surfFile = surf_files{idx_surf_hemi};
        g = gifti(fullfile(EEG.filepath, surfFile));
        tH = ieeg_RenderGifti(g);
        tH.FaceAlpha = 0.1;
        s = scatter3([EEG.chanlocs(idx_elecs_hemi).X], [EEG.chanlocs(idx_elecs_hemi).Y], [EEG.chanlocs(idx_elecs_hemi).Z], 'o', 'Filled');
        s.SizeData = 10; 
        s.MarkerFaceColor = [.9 .5 .5]; 
        s.MarkerEdgeColor = [0 0 0];
        ieeg_viewLight(90,0)  % 180,90 (top/front); 90,0 (right)
    else
        disp("No electrodes in right hemisphere to plot")
    end

    axis equal off 
    % vis3d tight
    title("Visualization (using subject's pial surface file from Freesurfer)",'Interpreter','none');


else
    % Fallback if no Freesurfer file is available: dipfit standard_BEM (smooth)
    dipfit_root = fileparts(which('dipfitdefs'));
    assert(~isempty(dipfit_root), 'dipfit not found on path. Enable dipfit in EEGLAB.');
    std_dir = fullfile(dipfit_root, 'standard_BEM');
    f_sccn  = fullfile(std_dir, 'standard_vol_SCCN.mat');
    f_std   = fullfile(std_dir, 'standard_vol.mat');

    if exist(f_sccn,'file'), S = load(f_sccn); else, S = load(f_std); end
    assert(isfield(S,'vol'), 'standard_BEM volume not found in dipfit.');

    bnd = S.vol.bnd; if numel(bnd) > 1, bnd = bnd(end); end
    if isfield(bnd,'pnt'), Vb = double(bnd.pnt); else, Vb = double(bnd.pos); end
    if isfield(bnd,'tri'), Fb = double(bnd.tri); else, Fb = double(bnd.face); end

    patch('Faces',Fb,'Vertices',Vb, ...
        'FaceColor',[0.75 0.80 0.90],'EdgeColor','none','FaceAlpha',0.12);

    E = [[EEG.chanlocs.X]' [EEG.chanlocs.Y]' [EEG.chanlocs.Z]'];
    if median(abs(E(:))) < 2, E = E*1000; end  % units only
    s = scatter3(E(:,1), E(:,2), E(:,3), 'o', 'filled');
    s.SizeData = 10; s.MarkerFaceColor = [.9 .5 .5]; s.MarkerEdgeColor = [0 0 0];

    try icadefs; set(gcf, 'color', BACKCOLOR); catch; end  % eeglab color
    axis equal off 
    % vis3d tight
    camlight headlight; camlight right
    lighting gouraud; material dull
    view([-135 20]);
    title('Visualization using standard BEM template (fallback when no Freesurfer pial surface files is detected)','Interpreter','none');

end


