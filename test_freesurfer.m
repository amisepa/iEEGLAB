%% FreeSurfer run from MATLAB on macOS (M1)
% Subject: UCI29
% Data folder:
%   /Users/cedriccannard/Documents/MATLAB/eeglab/plugins/iEEGLAB/tutorial/dataset5
%
% This script
% 1. Sets environment variables for FreeSurfer inside MATLAB
% 2. Converts T1 NIfTI to 001.mgz in $SUBJECTS_DIR/<SUB>/mri/orig
% 3. Kicks off recon-all (autorecon1 first)

% 0. User inputs
clear; clc
SUB = 'UCI29';

FREESURFER_HOME = '/Applications/freesurfer/8.1.0';
SUBJECTS_DIR    = '/Users/cedriccannard/Documents/MATLAB/eeglab/plugins/iEEGLAB/tutorial';  % parent folder only
FS_LICENSE      = '/Users/cedriccannard/Documents/license.txt';

% 1. Env for this MATLAB session
setenv('FREESURFER_HOME', FREESURFER_HOME);
setenv('SUBJECTS_DIR', SUBJECTS_DIR);
setenv('FS_LICENSE', FS_LICENSE);
setenv('FSF_OUTPUT_FORMAT','nii.gz');
setenv('PATH', [fullfile(FREESURFER_HOME,'bin') ':' fullfile(FREESURFER_HOME,'fsfast','bin') ':' getenv('PATH')]);

% 2. Create subject directories
orig_dir = fullfile(SUBJECTS_DIR, SUB, 'mri', 'orig');
if ~exist(orig_dir, 'dir'); mkdir(orig_dir); end

% 3. Convert T1 to 001.mgz
t1_nii = fullfile(SUBJECTS_DIR, SUB, 'SubjectUCI29_MR_acpc.nii');
assert(isfile(t1_nii), 'T1 NIfTI not found');
t1_mgz = fullfile(orig_dir, '001.mgz');
if ~isfile(t1_mgz)
    cmd = sprintf('mri_convert "%s" "%s"', t1_nii, t1_mgz);
    disp(cmd);  assert(system(cmd)==0, 'mri_convert failed');
end

% 4. CT to MGZ (Optional)
ct_nii = fullfile(SUBJECTS_DIR, SUB, 'SubjectUCI29_CT_acpc_f.nii');
if isfile(ct_nii)
    ct_mgz = fullfile(SUBJECTS_DIR, SUB, 'mri', 'CT.mgz');
    if ~isfile(ct_mgz)
        cmd = sprintf('mri_convert "%s" "%s"', ct_nii, ct_mgz);
        disp(cmd);  assert(system(cmd)==0, 'mri_convert CT failed');
    end
end

%%

SUB = 'UCI29';

% env
setenv('FREESURFER_HOME', '/Applications/freesurfer/8.1.0');
setenv('FREESURFER',      '/Applications/freesurfer/8.1.0');  % required by tcsh scripts
setenv('SUBJECTS_DIR',    '/Users/cedriccannard/Documents/MATLAB/eeglab/plugins/iEEGLAB/tutorial');
setenv('FS_LICENSE',      '/Users/cedriccannard/Documents/license.txt');
setenv('FSF_OUTPUT_FORMAT','nii.gz');
setenv('PATH', [fullfile(getenv('FREESURFER_HOME'),'bin') ':' fullfile(getenv('FREESURFER_HOME'),'fsfast','bin') ':' getenv('PATH')]);

setenv('OMP_NUM_THREADS','1');
setenv('ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS','1');
setenv('MKL_NUM_THREADS','1');
setenv('OPENBLAS_NUM_THREADS','1');

% run recon via tcsh and csh setup
cmd = sprintf('/bin/tcsh -c "source %s/SetUpFreeSurfer.csh; recon-all -sd \\"%s\\" -s \\"%s\\" -autorecon1 -nogcareg"', ...
              getenv('FREESURFER_HOME'), getenv('SUBJECTS_DIR'), SUB);
disp(cmd);
st = system(cmd);
if st ~= 0
    error('recon-all failed. See log at %s', fullfile(getenv('SUBJECTS_DIR'), SUB, 'scripts', 'recon-all.log'));
end


%% 5. Start recon-all (short pass first)
% 1) -sd points to parent folder
% 2) -s is only the subject ID

cmd = sprintf('/bin/tcsh -c "source %s/SetUpFreeSurfer.csh; recon-all -sd \\"%s\\" -s \\"%s\\" -autorecon1 -nogcareg"', ...
              getenv('FREESURFER_HOME'), getenv('SUBJECTS_DIR'), SUB);
disp(cmd);
st = system(cmd);
if st ~= 0
    error('recon-all failed. See log at %s', fullfile(getenv('SUBJECTS_DIR'), SUB, 'scripts', 'recon-all.log'));
end


%% 6. Full pipeline
% When ready for the full run, uncomment below. This can take many hours.

cmd = sprintf('recon-all -s "%s" -all', SUB);
fprintf('Running: %s\n', cmd);
st = system(cmd);
if st~=0
    error('recon-all -all failed. See log at %s', fullfile(SUBJECTS_DIR,SUB,'scripts','recon-all.log'));
end

%% 6. Quick sanity checks
% Show version and where binaries come from
system('recon-all --version');
system('which recon-all');
system(sprintf('mri_info "%s" | head -n 10', t1_mgz));

fprintf('\nDone. Subject folder: %s\n', fullfile(SUBJECTS_DIR,SUB));
