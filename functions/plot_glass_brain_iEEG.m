function h = plot_glass_brain_iEEG(elec_xyz, varargin)
% plot_glass_brain_iEEG  Minimal-dependency "glass brain" + iEEG electrodes
%
% INPUTS
%   elec_xyz   [N x 3] double, electrode XYZ coordinates (e.g., MNI mm)
%
% NAME-VALUE OPTIONS (all optional)
%   'Labels'       - {N x 1} cellstr electrode labels (default: {})
%   'MarkerSize'   - scalar, electrode marker size (default: 28)
%   'MarkerEdge'   - color, marker edge color (default: [0 0 0])
%   'MarkerFace'   - color, marker face color (default: [0.85 0.1 0.1])
%   'BrainAlpha'   - scalar [0..1], glass brain transparency (default: 0.08)
%   'BrainColor'   - color, glass brain color (default: [0.6 0.7 0.9])
%   'Center'       - [1x3] translation of brain shell (default: [0 0 0])
%   'Radii'        - [1x3] radii of ellipsoidal brain shell, mm (default: [70 95 75])
%   'SuperExp'     - scalar >=1, superquadric exponent (1.6 gives "brain-ier") (default: 1.6)
%   'ShowAxes'     - logical, show axes (default: false)
%   'ShowGrid'     - logical, add a faint 3D grid (default: false)
%   'LabelFontSize'- scalar, text size for labels (default: 10)
%
% OUTPUT
%   h  struct of handles (patch, scatter, text, axes)
%
% NOTES
%   1) Coordinates should be in a consistent space (e.g., MNI, RAS of your template).
%   2) This draws an analytic "brain-like" superquadric shell—fast and dependency-free.
%   3) If you have SPM/FreeSurfer, see the OPTIONAL block at the end for a real cortical mesh.

% 1. Parse inputs
p = inputParser;
p.addRequired('elec_xyz', @(x) ismatrix(x) && size(x,2)==3);
p.addParameter('Labels', {}, @(x) iscellstr(x) || isempty(x));
p.addParameter('MarkerSize', 28, @(x) isscalar(x) && x>0);
p.addParameter('MarkerEdge', [0 0 0], @(x) isnumeric(x) && numel(x)==3);
p.addParameter('MarkerFace', [0.85 0.1 0.1], @(x) isnumeric(x) && numel(x)==3);
p.addParameter('BrainAlpha', 0.08, @(x) isscalar(x) && x>=0 && x<=1);
p.addParameter('BrainColor', [0.6 0.7 0.9], @(x) isnumeric(x) && numel(x)==3);
p.addParameter('Center', [0 0 0], @(x) isnumeric(x) && numel(x)==3);
p.addParameter('Radii', [70 95 75], @(x) isnumeric(x) && numel(x)==3);   % ~MNI-ish head/brain scale
p.addParameter('SuperExp', 1.6, @(x) isscalar(x) && x>=1);
p.addParameter('ShowAxes', false, @(x) islogical(x) || ismember(x,[0 1]));
p.addParameter('ShowGrid', false, @(x) islogical(x) || ismember(x,[0 1]));
p.addParameter('LabelFontSize', 10, @(x) isscalar(x) && x>0);
p.parse(elec_xyz, varargin{:});
opt = p.Results;

% 2. Make or reuse figure/axes
ax = gca;
fig = ancestor(ax,'figure');
set(fig, 'Color', 'w');

% 3. Generate a smooth "superquadric" brain-like shell (dependency-free)
%    Equation: |x/a|^e + |y/b|^e + |z/c|^e = 1  (e ~1.6-2 gives a rounded brain-ish shape)
a = opt.Radii(1); b = opt.Radii(2); c = opt.Radii(3);
e = opt.SuperExp;

% Create a uniform grid in normalized coords, then map to mm
n = 80; % resolution of the shell surface
u = linspace(-1, 1, n);
[vx, vy, vz] = ndgrid(u,u,u);

F = (abs(vx).^e) + (abs(vy).^e) + (abs(vz).^e) - 1;   % implicit function = 0 on the surface
% Use isosurface on a *slightly negative* iso to keep shape a bit "puffed"
iso = -0.02;
S = isosurface(vx*a + opt.Center(1), vy*b + opt.Center(2), vz*c + opt.Center(3), F, iso);

% 4. Plot the translucent shell ("glass brain")
hold(ax, 'on');
hp = patch(ax, S);
hp.FaceColor = opt.BrainColor;
hp.EdgeColor = 'none';
hp.FaceAlpha = opt.BrainAlpha;
daspect(ax, [1 1 1]);
axis(ax, 'vis3d');
if ~opt.ShowAxes, axis(ax, 'off'); end
if opt.ShowGrid, grid(ax, 'on'); else, grid(ax, 'off'); end

% 5. Nice lights/materials
camlight(ax, 'headlight');
camlight(ax, 'right');
lighting(ax, 'gouraud');
material(hp, 'dull');

% 6. Plot electrodes
hs = scatter3(ax, elec_xyz(:,1), elec_xyz(:,2), elec_xyz(:,3), ...
    opt.MarkerSize, 'o', ...
    'MarkerEdgeColor', opt.MarkerEdge, ...
    'MarkerFaceColor', opt.MarkerFace, ...
    'LineWidth', 0.8);

% 7. Optional labels
ht = gobjects(0);
if ~isempty(opt.Labels)
    assert(numel(opt.Labels)==size(elec_xyz,1), 'Labels must match number of rows in elec_xyz.');
    ht = gobjects(size(elec_xyz,1),1);
    for k = 1:size(elec_xyz,1)
        ht(k) = text(elec_xyz(k,1), elec_xyz(k,2), elec_xyz(k,3), ...
            ['  ' opt.Labels{k}], 'FontSize', opt.LabelFontSize, ...
            'Color', [0.1 0.1 0.1], 'Interpreter','none');
    end
end

% 8. Camera/view presets that look "glass brain"-ish
view(ax, [-135 20]); % oblique
axis(ax,'tight');

% 9. Output handles
h = struct('axes', ax, 'brain', hp, 'elec', hs, 'labels', ht);

end
