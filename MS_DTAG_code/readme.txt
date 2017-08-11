MATLAB code to accompany manuscript titled:

"Characterizing Chilean blue whale vocalizations with DTAGs: a test of using tag accelerometers for caller identification" (2016-2017).

Authors:

Mark R. Saddler, Alessandro Bocconcelli, Leigh S. Hickmott, Gustavo Chiang, Rafaela Landea-Briones, Paulina A. Bahamonde, Gloria Howes, Paolo S. Segre, Laela S. Sayigh

Contents:

(1) dtagmark.m
(2) dtagsignal.m
(3) decimate_for_correlation.m
(4) dtagcorrelate.m
(5) dtagcorrelate_accelerations.m


------------------------------------------------------------------------------
I. dtagmark
------------------------------------------------------------------------------
Function generates graphic user interface for auditing DTAG data for low-frequency blue whale calls, especially for auditing synchronous hydrophone and accelerometer recordings. No parameters are required and running the function in MATLAB will open the GUI. The data directory must be specified in the code (line 14) and the code is designed to expect the following directory structure (dataDir = '...\data_bm\'):

    data_bm\
        audio\
            deployment folders containing .wav and .xml DTAG data files
        audit\
            .xlsx spreadsheets for each deployment (audit logs) following the format specified in 'BLANK_AUDIT_LOG.xlsx'
        prh\
            prh.mat files for each deployment (DTAG accelerometer data)

The GUI has buttons that prompt the user to load the audio data folder, audit log, and prh.mat file for the desired deployment. The accelerometer data can be bandpass filtered if appropriate values for the low and high frequency cutoffs are provided. If only a low frequency cutoff is provided, the accelerometer data will be high-pass filtered. To add an entry to the audit log, right-click on the audio spectrogram. This will generate cross-hairs which can be used to select the start and end of an event in the audio spectrogram (use left-click to select desired times). Once valid start and end times have been selected, the GUI will prompt the user to provide a call type label and a comment. Once an entry has been successfully to the audit log, it will display in the listbox in the side panel. Clicking on an entry in the listbox will pull it up in the audio and accelerometer records.


------------------------------------------------------------------------------
II. dtagsignal
------------------------------------------------------------------------------
Function generates graphic user interface for reviewing calls marked using dtagmark interface. The GUI is especially optimized for reviewing low-frequency blue whale calls with synchronous accelerometer data. The GUI allows for signals to be precisely selected and saved in a .mat file for further analysis and parameter calculations. The data directory must be specified in the code (line 12) and the output filename for the generated .mat file can be specified in line 37. The duration of noise samples collected pre- and post-signal can be specified in line 38. The GUI will open when the function is run. Buttons in the side panel prompt the user to load the deployment data (select the audio directory containing .wav and .xml files for a single deployment. NOTE: the prh.mat file for the deployment is expected to have the same name as the audio directory followed by 'prh.mat') and the audio log .xlsx file. Navigate to the events roughly marked using dtagmark by clicking on the log entries in the upper listbox. When the audio and accelerometer data are displayed for a log entry, specific calls can be selected from the window by right-clicking the audio spectrogram or waveform. Use the mouse's left-button to drag a rectangle around the desired signal and press 'Add signal to list' once happy with the selection. This will add the signal to the lower listbox. Note that the output .mat file, containing a cell array named 'savedSignals', will not be saved automatically. Press 'Save signals: ...' button after adding each signal to avoid losing progress if the GUI crashes. The output variable 'savedSignals' is an N-by-2 cell array where N is the number of signals selected. Each cell contains a data structure corresponding to a single call with fields for the signal and time vectors as well as metadata (call type, comment, source .wav filename, sample rates, etc.). Signals in column 1 of 'savedSignals' contain audio data and signals in column 2 of 'savedSignals' contain accelerometer data.


------------------------------------------------------------------------------
III. decimate_for_correlation
------------------------------------------------------------------------------
Simple script for decimating the acoustic data to the same sample rate as the accelerometer data, which is necessary for cross-correlation analyses between acoustic and accelerometer data. Simply load-in the 'savedSignals' variable saved using dtagsignal and run the decimate_for_correlation script. The script will add 'y_dec',  'A_dec', and 'fs_accel' fields to the audio signal column in 'savedSignals', containing the resampled signals as well as the new sampling rate. These fields are used by dtagcorrelate and dtagcorrelate_accelerations.


------------------------------------------------------------------------------
IV. dtagcorrelate
------------------------------------------------------------------------------
Function to analyze low-frequency baleen whale calls and use cross-correlation methods for phase analysis of accelerometer and acoustic signal. This function takes in the 'savedSignals' variable as a parameter and returns the 'signalParams' output data structure, which contains fields for a large variety of signal parameters calculated (see comments in script). This function also cross-correlates synchronous audio and accelerometer signals and saves the peak coefficients and lags in the 'signalParams' output data structure. If the variable input argument "save" is provided, then summary figures for each signal in 'savedSignals' will be saved in the specified output directory (line 22), which include spectrograms, power spectra, and cross-correlation plots. This function calls a number of auxiliary functions in this repository to generate and modify plots for the figures:
    MS_annotate_axes.m
    MS_powerspectrum_noise_correction.m
    MS_spectrogram.m
    MS_waveform.m


------------------------------------------------------------------------------
V. dtagcorrelate_accelerations
------------------------------------------------------------------------------
Function to calculate cross-correlations between all possible pairs of acceleration signals (inter-call cross-correlations). This function takes in the 'savedSignals' variable as a parameter and will output data structures containing matrices of the peak cross-correlation coefficients and lags from cross-correlations. The returned data structures have fields corresponding to each of the DTAG accelerometer's axes (X, Y, Z). The returned 'output' variable contains a convenient N-by-5 matrix with rows containing the peak coefficient values from all unique cross-correlation pairs (autocorrelations excluded) in the format: [index1, index2, coefX, coefY, coefZ]. The 'output' variable is useful for identifying groups of calls that exhibit high intra-correlation (one easy way is to paste the entire matrix into excel and sort rows by coefficient values). NOTE: this function can be easily adapted to cross-correlate noise or audio signals as well as accelerometer signals between calls.



MS 2017.08.10
