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

#### get_rpeaks.py [-h] [-i INPUT_FILE] [-o OUTPUT_DIR] [-c CHANNEL_NAME] [-p]

Process PPG file to get R peaks.

##### optional arguments:

  -h, --help                                      show this help message and exit
  
  -i INPUT_FILE, --input_file INPUT_FILE          path to PPG file.
                        
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR          directory to store output files. default: ./output
                        
  -c CHANNEL_NAME, --channel_name CHANNEL_NAME    name of PPG channel. default: PPG
                        
  -p, --plot                                      option to plot time series of heart rate and detected R peaks
  
  
##### Example:

> python get_rpeaks.py -i subj1_run1.vhdr  

Compute R peaks from the channel named 'PPG'. The R peaks detection results are saved in ./output

> python get_rpeaks.py -i subj1_run1.vhdr -o 'preprocessed_file/subj1' -c 'ppg' -p       

Compute R peaks from the channel named 'ppg'. The R peaks detection results and the figures are saved in ./reprocessed_file/subj1. 

   
