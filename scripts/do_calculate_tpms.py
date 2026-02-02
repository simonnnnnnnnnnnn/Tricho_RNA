import pandas as pd
import sys
import os

filename = sys.argv[1]
file = pd.read_csv(filename, sep="\t", comment = "#")

#file = pd.read_csv("~/Documents/uni/semester8/praxis/rna/fikt_A_marked_counts.txt", sep="\t")

# check
#print(file.columns)
"""
colname_list = list(file.columns)
temp_name = colname_list[-1]
final_name = temp_name.split("_")[0] + "_" + temp_name.split("_")[1] + "TPMs.csv" # get the name for the resulting file
"""
# better with os
variant_name = os.path.basename(filename).replace("_marked_counts.txt", "")
final_name = f"{variant_name}_TPMs.csv"

# plan: calculate RPK, then use that to calculate TPM from it

# rename last column (featureCounts names that after the filename --> a uniform name is desired here)
file.columns = ['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length', 'Counts'] # --> "Counts" introduced in place of the name of the bamfile



# divisions by zero should be avoided, therefore drop all rows where length = 0
file = file[file['Length'] > 0]

# now calculate the RPKs
file['RPK'] = (file['Counts'] / (file['Length'] / 1000)) # here reads are divided by kilobase... so RPK :D

# a scaling factor needs to be used in order to calculate the TPMs
# --> it is the sum of all RPKs divided by 1e6 --> per sample scaling factor (because RPKs may differ between samples)
factor = (file['RPK'].sum() / 1e6)

# now the TPMs can finally be calculated
file['TPM'] = (file['RPK'] / factor)

# now save it all to a new file --> no replacing in place
file.to_csv(final_name, sep="\t", index=False)

