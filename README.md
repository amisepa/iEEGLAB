# iEEGLAB
EEGLAB plugin for analyzing intracranial EEG (iEEG) data.

Adapted from work by Dora Hermes. 
Please cite the following references when using this plugin:




## Requirements

- MATLAB
- EEGLAB
- Mac Os: Xcode




MACOS:
for reading MEFD files, MATLAB needs to build the MEX file 
for that you need to:
1) install xcode
	Open Mac App Store → Search for Xcode → Install.

	Switch to full Xcode developer directory:
	sudo xcode-select --switch /Applications/Xcode.app/Contents/Developer

	Accept the license:
	sudo xcodebuild -license
	press (to display license): enter 
	type: agree

3) Verify MATLAB detects the compiler (back in matlab)
	mex -setup C
	You should now see an option like: Xcode with Clang

## FreeSurfer

1) Install XQuartz: https://www.xquartz.org/

2) Download FreeSurfer from the official site: https://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall
3) MAC: Right-click the .pkg file in Downloads → Open
--> If macOS Gatekeeper blocks the FreeSurfer installer:
	2.1) go to: Apple Menu → System Settings → Privacy & Security
	2.2) Scroll down until you see a message like: “‘freesurfer-macOS-darwin_arm64-8.1.0.pkg’ was blocked because it is not from an identified developer.”
	2.3) Click “Allow Anyway”.
	2.4) Control-click the .pkg file again and choose Open → then click Open in the new dialog. It will now install normally.

4) once FreeSurfer is installed, you’ll want your terminal to automatically load its environment variables (like FREESURFER_HOME, SUBJECTS_DIR, etc.) every time you open a new Terminal window.

