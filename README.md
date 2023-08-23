# CIN_analysis

## Data preparation
annotate_cytoband_data.ipynb 
Chromosome_data.ipynb
These scripts were used to annotate the data for the CIN computation and are only needed if the files 'cytoband_UCSC_annotated.csv' and 'chromosome_data.csv' not exist, because these are needed for the computation.



## Computation CIN

### Classes
Patient.py
Chromosome.py
ChrArm.py
Cytoband.py

### Compute the CIN
The functions are all in mai_cin.ipynb (wanted to move the function to the cin_functions.ipynb but these are older ones)
At the end of the Notebook one can find the main function to let the computation run (cells containing the functions need to be run first)

## HRR cytobands
cind_hrd_cytoband.ipynb
This can be used to reduce the size of the cytoband data as it extracts all cytobands that have a HRR related gene on them.
It also creates new Patients as well and both are stored with the addition of 'HRR' for example there 'general_cin_hrr' would be general cin on the sample level only taking into account HRR cytobands.

## Analysis CIN
Plotting and statistical analysis are done in analysis_CIN.R. The first function are for loading and there are functions for the heatmap, the are small parts for boxplots and functions to write the analysis txt files (will be updated with more documentation at the end of October)
