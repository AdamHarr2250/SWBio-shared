# Interrogating data from a small molecule screen for antibacterial activity against methicillin-resistant _Staphylococcus aureus_

## Introduction to dataset and code

The script in this repository aims to:

1) **Format** data to enable correlation and distribution analysis
2) **Calculate** properties from the measured replicate outputs (average and range of replicates) and the chemical structure (electrophilic warhead)
3) **Highlight putative 'hits'** based on the distribution of biological effects
4) Interrogate the **pairwise correlation and distribution of different parameters** including chemical properties and biological effects in the population of hits compared to the entire screened library

# Instructions

1) Before running the code, please ensure you have first downloaded the correct modules for working with this data:
<br>python3 -m pip install pandas
<br>python3 -m pip install seaborn
<br>python3 -m pip install matplotlib
<br>python3 -m pip install numpy
<br>python3 -m pip install rdkit # If this does not work, please find documentation on other ways to install the module at this [web address](https://www.rdkit.org/docs/Install.html)
3) 
