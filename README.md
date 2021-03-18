# PhysioKit

## Installation

1. install pipenv (https://pypi.org/project/pipenv/)

2. Download the source code of PhysioKit 

3. Change current path to the location of PhysioKit 

4. Run the following command: 

> pipenv install

## Support data format

1. BIOPAC (.acq file) 

2. BrainVision (.vhdr, .vmrk, .eeg)

## Function

### get_rpeaks.py [-h] [-i INPUT_FILE] [-o OUTPUT_DIR] [-c CHANNEL_NAME] [-p]

Process PPG file to get R peaks.

#### Arguments:

##### -h, --help
Show this help message and exit
  
##### -i INPUT_FILE, --input_file INPUT_FILE
Path to PPG file. (The extension of INPUT_FILE must be .acq or .vhdr)
                        
##### -o OUTPUT_DIR, --output_dir OUTPUT_DIR
Directory to store output files. default: ./output
                        
##### -c CHANNEL_NAME, --channel_name CHANNEL_NAME
Name of PPG channel. default: PPG
                        
##### -p, --plot
Option to plot time series of heart rate and detected R peaks
  
#### Outputs:

##### 1. preprocessed signals:

[OUTPUT_DIR]/[CHANNEL_NAME]/signal/[CHANNEL_NAME]_ preproc_[INPUT_FILE name]

1. csv file: raw ppg signals, filtered ppg signals, heart rate, and detected r peaks
2. txt file: the indices (start from 0) and time (unit:second) of every detected r peaks
3. pkl file: including the information in csv and txt files


##### 2. figures (only if set up the argument -p):

1. [CHANNEL_NAME]/fig/hr/[INPUT_FILE name]_hr.png: time series of heart rate

2. [CHANNEL_NAME]/fig/rpeak/[INPUT_FILE name]/rpeak_and_hr_sec_[START TIME].png: the detected r peaks (red dots) overlaid on the filtered ppg signals, which are segmented by a 90-second interval. 

#### Example:

> python get_rpeaks.py -i subj1_run1.vhdr  

Compute R peaks from the channel named 'PPG'. The R peak detection results are saved in ./output

> python get_rpeaks.py -i subj1_run1.vhdr -o 'preprocessed_file/subj1' -c 'ppg1' -p       

Compute R peaks from the channel named 'ppg1'. The R peak detection results and the figures are saved in ./preprocessed_file/subj1. 
