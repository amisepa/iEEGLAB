# iEEGLAB

EEGLAB plugin for analyzing intracranial EEG (iEEG) data. Supports both stereoEEG (sEEG) and eCoG data. 
The plugin supports continuous iEEG data applications:

	- Epilepsy research
	
	- Clinical monitoring

Although it is mainly designed for event-related applications:

	- Stimulus-/Response-locked
	
	- Cortico-Cortical Evoked Potentials (CCEP; e.g., single pulse stimulation experiments)

A lot of the code and algorithms implemented were adapted from work by Dora Hermes and the Multimodal Neuroimaging Lab (https://github.com/MultimodalNeuroimagingLab). 
Please cite the following references when using this plugin: 
[to update]



## Requirements

- MATLAB

- EEGLAB (installed and path added to MATLAB)
  
- Vistasoft (for visualizations; cloned/downloaded and path added to MATLAB): https://github.com/vistalab/vistasoft.git
  
- Data importation plugins (depends on the data format; e.g., .mefd, .edf, .vhdr, etc.)
  
- Cartesian (XYZ) electrode locations (in .tsv or .csv file)
  
- Events either directly in the data or in a .tsv or .csv file (for event-related applications)


### FreeSurfer: only if you need to generate new brain surface files for your subject from MRI data

1) Install XQuartz: https://www.xquartz.org/

2) Download FreeSurfer from the official site: https://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall
3) MAC: Right-click the .pkg file in Downloads → Open
--> If macOS Gatekeeper blocks the FreeSurfer installer:
	2.1) go to: Apple Menu → System Settings → Privacy & Security
	2.2) Scroll down until you see a message like: “‘freesurfer-macOS-darwin_arm64-8.1.0.pkg’ was blocked because it is not from an identified developer.”
	2.3) Click “Allow Anyway”.
	2.4) Control-click the .pkg file again and choose Open → then click Open in the new dialog. It will now install normally.

4) once FreeSurfer is installed, you’ll want your terminal to automatically load its environment variables (like FREESURFER_HOME, SUBJECTS_DIR, etc.) every time you open a new Terminal window.

