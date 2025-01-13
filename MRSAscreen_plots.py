# Import all modules required for working with this data: 
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import rdkit

# Import the data from the csv file uploaded to GitHub (direct import)
df = pd.read_csv("https://raw.githubusercontent.com/AdamHarr2250/SWBio-shared/refs/heads/main/SFFdata.csv")

# Specify a folder in the current directory to save plots into. If folder doesn't already exist, one is made in current working directory. Folder name is user-specified to enable users to prevent overwriting.
plot_folder = input("What would you like the folder to be called for your plots?")
if not os.path.exists(plot_folder):
    os.makedirs(plot_folder)

# Create a grid to show the initial descriptive plots together
fig, axes = plt.subplots(1, 2, figsize=(15, 5))

# Plot data into each of the grid subplots
# First, a regression plot with x and y as columns "% inh rep1" and "% inh rep2"; biological replicate measurements of cell growth after incubation. 
sns.regplot(ax = axes[0], x = df["% inh rep1"], y = df["% inh rep2"], data = df, line_kws={"color": "blue", "lw": 0.5}, label = "SuFEx fragments")

# Label the title and axes
axes[0].set_title("Scatter Plot with Line of Best Fit")
axes[0].set_xlabel("Inhibition of MRSA cell growth n=1 (%)")
axes[0].set_ylabel("Inhibition of MRSA cell growth n=2 (%)")

# Dynamically define the minimum and maximum values
min_val = min(df['% inh rep1'].min(), df['% inh rep2'].min())
max_val = max(df['% inh rep1'].max(), df['% inh rep2'].max())

# Add a line of identity between the minimum and maximum values for visual comparison to the regression line
axes[0].plot([min_val,max_val], [min_val,max_val], color='red', linestyle='--', linewidth=.5, label='Identity Line')
axes[0].legend()


# The distribution of the two datasets can also be shown overlaid to confirm similarity.
sns.histplot(ax = axes[1], data = df, x = "% inh rep1", color = "blue", binwidth = 2, label = "Repeat 1", legend = True, binrange = [min_val, max_val], kde = True)
sns.histplot(ax = axes[1], data = df, x = "% inh rep2", color = "red",  binwidth = 2, label = "Repeat 2", legend = True, binrange = [min_val, max_val], kde = True)

# Define the title and x axis label for the overlaid histogram
axes[1].set_title("Histogram of replicate MRSA growth screens")
axes[1].set_xlabel("Relative inhibition of MRSA cell growth (%)")
axes[1].legend()

plt.tight_layout()

# Save the plot and open it
reghistplot_filename = os.path.join(plot_folder, "reghist_plot.png")
plt.savefig(reghistplot_filename)
plt.show()

# The replicated data is expected to be colinear and can therefore be summarized by taking an average, "% inh average (n=2)" and range, "% inh range (n=2)".
df["% inh average (n=2)"] = (df["% inh rep1"] + df["% inh rep2"])/2
df["% inh range (n=2)"] = abs(df["% inh rep1"] - df["% inh rep2"])

# The distribution of average growth inhibition between the two replicates can be shown to be similar to each replicate
sns.histplot(data = df, x = "% inh average (n=2)", color = "purple",binwidth = 2, label = "Average MRSA growth inhibition across two replicates (%)", legend = True, binrange = [min_val, max_val], kde = True)

# Set title and xlabel of average inhibition of MRSA growth
plt.title("Histogram of average inhibition of MRSA growth")
plt.xlabel("Average MRSA growth inhibition across two replicates (%)")

# Save to same folder and show
avhistplot_filename = os.path.join(plot_folder, "avhist_plot.png")
plt.savefig(avhistplot_filename)
plt.show()


# It might also be of interest to know which electrophilic warhead is contained within each fragment to establish if there is a strong preference for one or the other.
from rdkit import Chem
from rdkit.Chem import Draw

electrophile1 = "S(=O)(=O)(-F)"
electrophile2 = "S(=O)(=O)(-F)-O"

substruc1 = Chem.MolFromSmarts(electrophile1)
substruc2 = Chem.MolFromSmarts(electrophile2)

smiles_list = df["SMILES"]

# all fragments will contain O=S(=O)F but a subselection will contain O=S(=O)(O)F. Here, we only search for the presence of O=S(=O)(O)F and add a True/False column for whether each fragment contains the O=S(=O)(O)F substructure. If False, it can be assumed to contain O=S(=O)F instead. Applies a lambda function to dataframe column "SMILES" to iterate through each row, perform substructure search and return boolean. Note that the 'Chem.MolFromSmiles(x)' is important to transform the SMILES code into an rdkit Mol data format.
df["electrophile2"] = df["SMILES"].apply(lambda x: Chem.MolFromSmiles(x).HasSubstructMatch(substruc2)) 

# Isolate putative active fragments based on average % inhibition
# First, summary statistics must be collected on the % inh average (n=2) column
mean = df["% inh average (n=2)"].mean()
stdev = df["% inh average (n=2)"].std()
print(f"Average of % inhibition: {mean:.2f} %\nStandard deviation of % inhibition: {stdev:.2f} %") #.2f limits the float values to 2 decimal points

# We can then find the threshold value for a compound to be considered a 'hit'. See methods for explanation of the below calculation.
signifthreshold = mean + 3* stdev
print(f"Significance threshold for anti-MRSA activity: {signifthreshold:.2f} %")

# Add column for activity cataegory: 'active' or 'inactive'
average = [df["% inh average (n=2)"]]

def activitycat(average):
    if average > signifthreshold:
        return "active"
    else:
        return "inactive"

df["activity"] = df["% inh average (n=2)"].apply(activitycat)

# Create a new subsection of the dataframe containing only active hits
hits = df[df["activity"] == "active"]

# This dataframe subsection can then be reordered in order of potency
hits = hits.sort_values("% inh average (n=2)", ascending = 0)

# Extract each parameter for hit visualisation into a list
hit_list_smiles = hits["SMILES"].tolist() # Creates a list object containing the SMILES codes for the hits
hit_cat_no = hits["Catalog ID"].tolist() # Creates a list object containing the unique identifier for the hits
hit_list_avpercentinh = hits["% inh average (n=2)"].tolist() # Creates a list object containing the average inhibition of MRSA cell growth
trunc_avpercentinh = [float(f"{x:.2f}") for x in hit_list_avpercentinh] # Truncates each average percentage inhibition to 2 decimal places for clarity
hit_list_inhrange = hits["% inh range (n=2)"].tolist() # Creates a list object containing the percentage range for the hits
trunc_inhrange = [float(f"{x:.2f}") for x in hit_list_inhrange] # Truncates each percentage inhibition range to 2 decimal places for clarity

# Converts SMILES to rdkit SMILES format for visualisation
hit_smiles = [Chem.MolFromSmiles(smiles) for smiles in hit_list_smiles]

# Create a legend format that combines the above lists
hit_legends = [f"{cat_no}\nAverage inhibition: {inh} %\nRange inhibition: {inhrange} %" for cat_no, inh, inhrange in zip(hit_cat_no, trunc_avpercentinh, trunc_inhrange)]

# Display the hits as 2D structures with legends containing unique identifiers and average and range of MRSA growth inhibition
hit_structures = Draw.MolsToGridImage(
    hit_smiles, 
    molsPerRow = 5, 
    subImgSize = (400, 400),
    legends = hit_legends) # creates a grid view of the structures of the chemical hits

# Save hits figure to same directory as before and show
hits_figure_filename = os.path.join(plot_folder, "hits_figure.png")
hit_structures.save(hits_figure_filename)
hit_structures.show()


# Set up a correlation matrix of biological measurements, % inh average (n=2) and % inh range (n=2), against chemical parameters for all data and actives only to identify whether there is any chemical feature that correlates with activity.
# Remove unwanted columns:
# "Structure [idcode], "Formula", and "Chemical name" contains the same information as SMILES and are not easily read by the modules to be used in the current script
# "Plate_ID", "Well", "V, ÂµL", "PO", Conc, mM", "Purity" are all values associated with the physical storage of the compounds and are not variables associated with the final results from the experiment
# "Stereochem.data" and "Geometric.isomer" may be important but are not relevant to a majority of the compounds so have been omitted for now from further analysis
df = df.drop(columns = ["Structure [idcode]","Plate_ID","Well","V, ÂµL","PO","Conc, mM","Formula","Purity",'Stereochem.data', 'Geometric.isomer',"Chemical name"])

# Some of the remaining columns contain non-numerical data types due to commas replacing decimal points. 
# Replace all commas in ClogP,logS and TPSA columns with decimal points:
columnstomodify = ["ClogP", "logS", "TPSA"]

df[columnstomodify] = df[columnstomodify].apply(lambda col: col.str.replace(',', '.').astype(float))

# Catalog ID is an identifier string and should not be considered in any correlation. Other identifying factors also don't make sense to include: "Structure No". "MW_salt" is irrelevent since the salt ion is not expected to be the active species. % inh rep1 and % inh rep 2 are known to be colinear and are now represented by average % inh (n=2) and % inh range (n=2)
# An intermediate dataframe, "selected_columns" is created containing only the factors we are interested in correlating.
dfcolumns = df.columns.tolist()
itemstoremove = ["Structure No","Catalog ID","SMILES","MW_salt","% inh rep1", "% inh rep2", "activity"]
selected_columns = [item for item in dfcolumns if item not in itemstoremove]

# 2 subsets of the selected_columns list are then created to enable filtering of correlation data in the matrix
biologicalcols = ["% inh average (n=2)", "% inh range (n=2)"]
correlating_columns = [item for item in selected_columns if item not in biologicalcols]

# Calculate the correlation between each pair of variables in each subset of data including all data and actives only
corr_matrix = df[selected_columns].corr()
corr_matrix_actives = df[selected_columns][df["% inh average (n=2)"] > signifthreshold].corr()
corr_matrix_inactives = df[selected_columns][df["% inh average (n=2)"] < signifthreshold].corr()

# Specifies a section of the matrix
corr_section = corr_matrix.loc[biologicalcols, correlating_columns]
corr_section_actives = corr_matrix_actives.loc[biologicalcols, correlating_columns]
corr_section_inactives = corr_matrix_inactives.loc[biologicalcols, correlating_columns]

# Create a grid to show correlation heatmaps concurrently
fig, axes = plt.subplots(2, 1, figsize=(10, 5))

sns.heatmap(corr_section, annot=True, cmap="coolwarm", fmt=".1f", linewidths=0.5, vmin=-1, vmax=1, ax=axes[0])
sns.heatmap(corr_section_actives, annot=True, cmap="coolwarm", fmt=".1f", linewidths=0.5, vmin=-1, vmax=1, ax=axes[1])
axes[0].set_title("Correlation Matrix for all activities")
axes[1].set_title("Correlation Matrix for actives")

plt.tight_layout()

# Save correlation matrix to previously specified folder and show
corr_matrix_filename = os.path.join(plot_folder, "corr_matrix.png")
plt.savefig(corr_matrix_filename)
plt.show()