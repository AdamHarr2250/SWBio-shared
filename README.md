# Interrogating data from a small molecule screen for antibacterial activity against methicillin-resistant _Staphylococcus aureus_

## Introduction to dataset and code

This script, 'smallmol_screen.py' ingests a csv file containing chemical information about a library of covalent fragments from the chemical supplier, Enamine, and biological antibacterial data recorded by our lab against methicillin-resistant Staphylococcus aureus (MRSA). The script processes the data and outputs distribution plots of the biological responses, dynamically and statistically determines a threshold for inhibitory activity that is significantly different from inactivity and outputs a visual grid of putative active chemical structures with unique ID and corresponding biological response data. All plots and figures are saved to a folder in the current working directory.

The script in this repository aims to:

1) **Format** data to enable correlation and distribution analysis
2) **Calculate** properties from the measured replicate outputs (average and range of replicates) and the performs a chemical substructure search (electrophilic warhead)
3) **Highlight putative 'hits'** based on the distribution of biological effects
4) Interrogate the **pairwise correlation and distribution of different parameters** including chemical properties and biological effects in the population of fragments categorised as active compared to the entire screened library

# Instructions

<br> 
1) Before running the code, please ensure you have first downloaded the correct modules for working with this data. Copy and paste the following lines into your command prompt.
<br>python -m pip install pandas
<br>python -m pip install seaborn
<br>python -m pip install matplotlib
<br>python -m pip install numpy
<br>python -m pip install rdkit 
If this does not work for RDkit, please find documentation on other ways to install the module at this web address: https://www.rdkit.org/docs/Install.html
<br>
<br>
2) Please find hyperlinks to the script, data file and GitHub repository in the attached .docx report and here:
   <br> GitHub repository: https://github.com/AdamHarr2250/SWBio-shared/tree/main
   <br> Script: https://github.com/AdamHarr2250/SWBio-shared/blob/main/MRSAscreen_plots.py
   <br> Data: https://raw.githubusercontent.com/AdamHarr2250/SWBio-shared/refs/heads/main/SFFdata.csv (Note: The script contains a relative data reference and will directly import data from the GitHub repository without any need to manually download the data).
<br>
<br>
3) Ensure your workig directory is pointing to the same file path as where your script has been saved. This can be changed in the terminal: 'chdir("{filepath to folder containing script}")'
4) When the script is run using 'python smallmol_screen.py', it will first prompt the user to input a folder name which will be created (or found if it already exists) in the current working directory. All plots and figures will be saved to this file path. Plots and figures will also be opened as the script runs for immediate manual inspection.

## Outputs
1) reghist_plot.png: 2 by 1 grid containing a regression plot of two replicates of inhibition of methicillin-resistant MRSA growth in pane 1 and overlaid histograms of each replicate to compare distribution of biological response data.
2) avhist_plot.png: a histogram of the distribution of the average MRSA growth inhibition calculated from the two replicate values for every fragment tested.
3) Values are printed: mean and standard deviation of average inhibition of MRSA growth, and threshold calculated by mean + 3* stdev. Theory is included in report.
4) hits_figure.png: a grid of fragments classed as displaying significant average inhibition values.
5) corr_matrix: a heatmap of pearson correlation coefficient values between biolgical summary values (mean and range MRSA growth inhibition) against the chemical property values. This aims to identify if any single chemical property is predictive of biological output in the overall data and the active fragment set.
