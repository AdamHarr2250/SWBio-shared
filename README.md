# Interrogating data from a small molecule screen for antibacterial activity against methicillin-resistant _Staphylococcus aureus_

## Introduction to dataset and code

This script ingests a csv file containing chemical information about a library of covalent fragments from the chemical supplier, Enamine, and biological antibacterial data recorded by our lab. The script processes the data and outputs distribution plots of the biological responses, dynamically and statistically determines a threshold for activity and outputs a visual grid of putative active chemical structures with unique ID and corresponding biological response data.

The script in this repository aims to:

1) **Format** data to enable correlation and distribution analysis
2) **Calculate** properties from the measured replicate outputs (average and range of replicates) and the performs a chemical substructure search (electrophilic warhead)
3) **Highlight putative 'hits'** based on the distribution of biological effects
4) Interrogate the **pairwise correlation and distribution of different parameters** including chemical properties and biological effects in the population of fragments categorised as active compared to the entire screened library

# Instructions

<br> 
1) Before running the code, please ensure you have first downloaded the correct modules for working with this data. Copy and paste the following lines into your command prompt.
<br>python3 -m pip install pandas
<br>python3 -m pip install seaborn
<br>python3 -m pip install matplotlib
<br>python3 -m pip install numpy
<br>python3 -m pip install rdkit 
If this does not work for RDkit, please find documentation on other ways to install the module at this [web address](https://www.rdkit.org/docs/Install.html)
<br>
<br>
2) Please find hyperlinks to the script, data file and GitHub repository in the attached .docx report and here:
   <br> GitHub repository: [web address](https://github.com/AdamHarr2250/SWBio-shared/tree/main)
   <br> Script: [web address](https://github.com/AdamHarr2250/SWBio-shared/blob/main/MRSAscreen_plots.py)
   <br> Data: [web address](https://raw.githubusercontent.com/AdamHarr2250/SWBio-shared/refs/heads/main/SFFdata.csv) (Note: The script contains a relative data reference and will directly import data from the GitHub repository without any need to manually download the data).
<br>
<br>
3) When the script is run, it will first prompt the user to input a folder name which will be created or found (if it already exists) in the current working directory. All plots and figures will be saved to this file path. Plots and figures will also be opened as the script runs for immediate manual inspection.
