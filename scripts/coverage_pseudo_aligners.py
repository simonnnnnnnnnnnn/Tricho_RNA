import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import seaborn as sns
import json


info = sys.argv[1]
files = sys.argv[2:]

# sort the files
"""
kallisto = []
salmon = []

for i in range(len(files)):
    if "abundance" in files [i]:
        kallisto.append(files[i])
    elif "quant" in files[i]:
        salmon.append(files[i])
    else:
        print("unexpected file detected")
"""
openinfo = pd.read_csv(info, sep=",", header=0)
namepattern = openinfo['name']
print(namepattern)



def get_salmon_rate(infile):
    with open(infile, "r") as file:
        data = json.load(file)
        print(data)
    return data.get('percent_mapped', 0)

def get_kallisto_rate(infile):
    with open(infile, "r") as file:
        data = json.load(file)
        print(data)
    return data.get('p_pseudoaligned', 0)

# read and plot both files
def plot_coverage(pseudofiles, name):
    kallisto = None
    salmon = None
    for elements in pseudofiles:
        if "meta" in elements:
            salmon = get_salmon_rate(elements)
        elif "run" in elements:
            kallisto = get_kallisto_rate(elements)

    # lets plot
    print(kallisto)
    print(salmon)
    map_input =  []
    map_input.append({
        'Tool': 'Kallisto',
        'Mapping Rate': kallisto
    })
    map_input.append({
        'Tool': 'Salmon',
        'Mapping Rate': salmon
    })

    myframe = pd.DataFrame(map_input)
    print(myframe.head())


    sns.set_style("whitegrid")
    plt.figure(figsize=(10,6))
    sns.barplot(x='Tool', y='Mapping Rate', data=myframe)
    plt.ylabel('Mapping Rate in %')
    plt.title(f'Mapping Rate {name}')
    plt.tight_layout()
    plt.savefig(f'comparison mapping rate {name}.png')





for i in range(len(namepattern)): # iterate through all samples
    filelist = []
    for j in range(len(files)): # get the files belonging to each sample --> saved in filelist
        if namepattern[i] in files[j]:
            filelist.append(files[j])
    # additionally append the namepattern, so the plotting function can name the plot
    #filelist.append(namepattern[i])

    # now filelist should have 2 files: one quant and one abundance
    plot_coverage(filelist, namepattern[i]) # logic must be copied to the end...


