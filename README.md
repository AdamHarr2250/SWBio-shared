# Interrogating data from a small molecule screen for antibacterial activity against methicillin-resistant _Staphylococcus aureus_

## Introduction to dataset and code

Infections arising from pathogenic bacteria are typically combatted with antibiotics but resistace to these drugs is on the rise. Multidrug-resistant bacterial species such as methicillin-resistant _Staphylococcus aureus_ (MRSA) exhibit resistance to most of the drugs currently used making them difficult to treat and leading to greater disease burden and, in many cases, death. Resistance mechanisms arise through random mutagenesis that confer specific resistance to a class of related antibiotics thereby improving the fitness of these bacteria when subjected to antibiotics such as oxacillin. Unfortunately, there have been very few novel chemotypes developed as new antibiotics in the last 40 years meaning that there are often very few alternative treatments when resistance arises. It is therefore imperative to identify new chemistry that can be used effectively as antibiotic agents. Herein, growth inhibition data from an initial screen of a compound library against live culture of MRSA is to be interrogated to highlight potential chemical starting points for antibiotic development. The compound library consisted of 960 diverse low molecular weight 'fragments', each containing a stable electrophilic warhead capable of forming irreversible covalent bonds with nucleophilic side chains of specific amino acids under physiological conditions. Each fragment was tested at 100 uM final assay concentration and replicates were performed in singlicate in independent experiments. All raw data was normalised to on-plate controls and a reference antibiotic, oxacillin, was used to confirm consistent pharmacological performance between plates. A range of growth inhibition values was measured which must be deconvoluted to identify a statistically-based threshold, apply this to the screen to label actives, and then look to see if there are any described chemical properties that are more likely to be found in an active fragment in this assay.

## Aims

The script in this repository aims to:

1) **Format** data to enable correlation and distribution analysis
2) **Calculate** properties from the measured replicate outputs (average and range of replicates) and the chemical structure (electrophilic warhead, substructure searches)
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
