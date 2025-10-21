# iEEGLAB

EEGLAB plugin for analyzing intracranial EEG (iEEG) data. Supports both stereoEEG (sEEG) and eCoG data. 
The plugin supports continuous iEEG data applications:

- Epilepsy research
	
- Clinical monitoring


Although it is mainly designed for event-related applications:

- Stimulus-/Response-locked
	
- Cortico-Cortical Evoked Potentials (CCEP; e.g., single pulse stimulation experiments)


A lot of the code and algorithms implemented in this plugin were adapted from work by Dora Hermes and the Multimodal Neuroimaging Lab (https://github.com/MultimodalNeuroimagingLab). 
Please cite the following references when using this plugin: 

- Valencia, G. O.et al., (2023). Signatures of electrical stimulation driven network interactions in the human limbic system. Journal of Neuroscience, 43(39), 6697-6711. https://pubmed.ncbi.nlm.nih.gov/37620159/

- Huang, H., Valencia, G. O., Gregg, N. M., Osman, G. M., Montoya, M. N., Worrell, G. A., ... & Hermes, D. (2024). CARLA: Adjusted common average referencing for cortico-cortical evoked potential data. Journal of neuroscience methods, 407, 110153. https://pubmed.ncbi.nlm.nih.gov/38710234/

- Miller, K. J., Müller, K. R., Valencia, G. O., Huang, H., Gregg, N. M., Worrell, G. A., & Hermes, D. (2023). Canonical Response Parameterization: Quantifying the structure of responses to single-pulse intracranial electrical brain stimulation. PLoS computational biology, 19(5), e1011105. https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011105


## Requirements

- MATLAB installed with license

- EEGLAB Toolbox installed and path added to MATLAB: https://github.com/sccn/eeglab?tab=readme-ov-file#installingcloning

- iEEGLAB plugin

	option 1: download + unzip (or clone if you use git version control; ideal to pull up to date fixes and improvements before new releases are made) this iEEGLAB repository into your local eeglab > plugins folder

	option 2 (via EEGLAB user interface): open EEGLAB > File > Manage extensions > type "iEEGLAB" in the search bar > Install

Note: eeglab will add the paths to the plugin and functions automatically when you launch EEGLAB. Otherwise, you can also clone/place the iEEGLAB folder somewhere else on your computer and manually add the path to the folder in MATLAB (if preferred). 

- Vistasoft for visualizations: clone/download repository (https://github.com/vistalab/vistasoft.git) and add the path in MATLAB manually (see above)
  
- Data importation plugins (depends on the data format; e.g., .mefd, .edf, .vhdr, etc.; See tutorial for more details). 
  
- Have 3D cartesian (XYZ) electrode locations for each file you wish to analyze (in .tsv/.csv file)
  
- Have events either directly in the data or in a .tsv/.csv file (for event-related applications)


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


