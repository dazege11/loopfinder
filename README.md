# Loopfinder
 
Code for the paper "LoopFinder: Sequence-Driven Prioritization of Enzyme Loop-Replacement Candidates".
 
The script 'Run_Prediction.py' was used to run the prediction in paper for the activity of Loop Replacement protein.

## Installation

You can installing the dependency through following command:

```
pip3 install -r requirements.txt
```

## Example Usage

```
python3.12 Run_Prediction.py --Input_Sequence input.tsv --Critical_Loop_Region 215-221 --Model_type DAAO --Blosum_Loop 10 --Output_File Predicted.csv
```
* `--Input_Sequence` -- Protein sequence input file.
* `--Critical_Loop_Region` -- Critical Loop Region Position in Protein sequence.
* `--Model_type` -- Select Predict Model for Prediction
* `--Blosum_Loop` -- Top-ranked loop sequences with the highest BLOSUM scores
* `--Output_File` -- Prediction Output File location.

## Input Structure

Input File Structure Should be Like:

```
MHSQKRVVVLGSGVIGLSSALILARKGYSVHILARDLPEDVSSQTFASPWAGAVWTPQMTLTDGPRQAKWEESTFKKWVELVPTGHAMWLKGTRRFAQNEDGLLGHWYKDITPNYRPLPSSECPPGAIGVTYDTLSVHAPKYCQYLARELQKLGATFERRTVTSLEQAFDGADLVVNATGLGAKSIAGIDDQAAEPVRGQTVLVKSPCKRCTSDSSDPASPAYIIPRPGGEVICGGTYGVGDWDLSVNPETVQRILKHCLRLDPTISSDGTIEGIEVLRHNVGLRPARRGGPRVEAERIVLPLDRTKSPLSLGRGSARAAKEKEVTLVHAYGFSSAGYQQSWGAAEDVAQLVDEAFQRYHGAARESKL
```


## Data Availablity

The raw data from the analyses in the paper can be found here: https://drive.google.com/drive/folders/1JEhmPBjn1A7GYO0D2YcNL9hkIFMGT-Dh?usp=sharing


