import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys


# read files as CL arguments
# the idea is to wrap everything in functions
# then the list of files can be worked through
# although with nextflow i could just call this process for all files
# functions may be better if i want all files at once to compare or bundle them together
# --> so i would get an average for all files, not sure of that even necessary or makes a lot of sense

aligner = sys.argv[1]
file = sys.argv[2] # the rest of the args
name = sys.argv[3]

def read_file(infile):
    stats = {} # the general infos are stored in a dict --> easier to compare later on with key:value storage
    mq = [] # here the mapping quality info gets saved

    with open(infile, "r") as file:
        for line in file:
            if line.startswith("SN"):
                elements = line.strip().split("\t") # should be tab separated
                key = elements[1].replace(':', '') # get the second elements (--> first one would be "SN"), those describe the stats
                stats[key] = elements[2] # insert the value
            elif line.startswith('MAPQ'):
                elements = line.strip().split("\t")
                mq.append({'mq': elements[1], 'counts': elements[2]})

    return stats, pd.DataFrame(mq) # returning the df is a bit more organized than just to pass on a list


def plot_it_all(current_file, outname):
    stats, mq = read_file(current_file)

    figure, axes = plt.subplots(1,2, figsize=(18,6))
    figure.suptitle(f'Mapping Quality {outname}', fontsize=16)

    # first of the three parts
    """
    axes[0].bar(mq['mq'], mq['counts'], alpha=0.7, color = 'magenta')
    axes[0].set_xlabel('mapping quality (phred)')
    axes[0].set_ylabel('reads')
    axes[0].set_title('Mapping Quality')
    axes[0].set_yscale('log') # otherwise that would look odd
    """

    # MAPQ plot
    num_slices = len(mq['counts'])
    colormap = matplotlib.colormaps['plasma'].resampled(num_slices)
    colors = colormap(np.arange(num_slices))
    axes[0].pie(mq['counts'], labels=mq['mq'], colors=colors, autopct='%1.1f%%', startangle=90)
    axes[0].set_title("Mapping Quality (MAPQ) Distribution")

    # now a pie chart for the mapped reads ... not super necessary, but looks nice
    """
    mapped_reads = stats.get('reads mapped', 0)
    total_reads = stats.get('raw total sequences', 1) # thats why python is better than R: there are proper dicts, just so good
    not_mapped = int(total_reads) - int(mapped_reads)
    labels = ['mapped', 'unmapped']
    amounts = [mapped_reads, not_mapped]
    colors = ['green', 'red']
    axes[1].pie(amounts, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90)
    axes[1].set_title("mapping rate")
    """

    # last but not least: the stats
    my_stats = {
        'Total Reads': stats.get('raw total sequences', 0),
        'Mapped Reads': stats.get('reads mapped', 0),
        'Properly Paired': stats.get('reads properly paired', 0),
        'Average Quality': stats.get('average quality', 0),
        'Average Length': stats.get('average length', 0)
    }

    names = list(my_stats.keys())
    vals = []#list(my_stats.values())
    for key, value in my_stats.items():
        vals.append(float(value))

    num_bars = len(vals)
    barcolors = matplotlib.colormaps['plasma'].resampled(num_bars)(np.arange(num_bars))
    axes[1].barh(names, vals, color = barcolors)
    axes[1].set_xlabel('count')
    axes[1].set_title('stats')

    for x,y in enumerate(vals):
        axes[1].text(y + max(vals)*0.01, x, f'{y:,.0f}', va='center', fontsize=10)

    plt.tight_layout()
    return figure # i really hope that all works

def plot_avg(*all_files):
    pass # here all files are read first --> make averages
    # maybe just handle all the averaging here, then pass it on to plotting function...


# now call it all
temp = plot_it_all(file, name)
temp.savefig(f'{aligner}_{name}_mq.png')