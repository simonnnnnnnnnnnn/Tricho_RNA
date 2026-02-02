import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

aligner = sys.argv[1]
file = sys.argv[2]
filename = sys.argv[3]

# optionally filelist and namelist

def genome_coverage(infile, outname):
    data = pd.read_csv(infile, sep="\t", names = ['chr', 'start', 'end', 'coverage'])

    # setup
    figure, axes = plt.subplots(2,2, figsize=(15,10))
    figure.suptitle(f'Genome Coverage {outname}', fontsize=16)

    # so this is the coverage distribution
    axes[0,0].hist(data['coverage'], bins=50, alpha=0.7, color='skyblue') # thats the nicer version of blue :D
    axes[0,0].axvline(data['coverage'].mean(), color='magenta', linestyle='--', label=f'mean: {data["coverage"].mean():.1f}')
    # now the labels
    axes[0,0].set_xlabel('coverage depth')
    axes[0,0].set_ylabel('frequency')
    axes[0,0].set_title('coverage distribution')
    axes[0,0].set_yscale('log') # get everything in the same range
    axes[0,0].legend()

    # next diagram: coverage and position
    selected_data = data.sample(n=min(10000, len(data)))
    selected_data['position'] = (selected_data['start'] + selected_data['end']) / 2

    for chromosome_name in selected_data['chr'].unique()[:10]: # first 10 chromosomes are shown, i fear that using everything will be too cluttered, on the other hand trichoplay only has 10..
        chrom_data = selected_data[selected_data['chr'] == chromosome_name]
        axes[0,1].plot(chrom_data['position'], chrom_data['coverage'], alpha=0.7, label=chromosome_name)

        axes[0,1].set_xlabel('genomic position')
        axes[0,1].set_ylabel('coverage depth')
        axes[0,1].set_title('coverage along the genome')
        axes[0,1].legend()

    # now some boxplots --> coverage
    chrom_coverage = []
    chrom_names = [] # must not be confused with chrom_name....maybe i should use better names
    for chrom_name in data['chr'].unique()[:10]:
        chr_data = data[data['chr'] == chrom_name]['coverage']# i really hope that works
        chrom_coverage.append(chr_data)
        chrom_names.append(chrom_name)

    axes[1,0].boxplot(chrom_coverage, labels=chrom_names)#plural
    axes[1,0].set_xlabel('chromosome')
    axes[1,0].set_ylabel('coverage depth')
    axes[1,0].set_title('coverage by chromosome')
    axes[1,0].tick_params(axis='x', rotation=45)# axes != axis

    # lets see if there are regions without coverage... unlikely but still
    no_cover = data[data['coverage'] == 0]
    no_cover['length']  = no_cover['end'] - no_cover['start']
    if len(no_cover) > 0:
        axes[1,1].hist(no_cover['length'], bins=30, alpha=0.7, color='red')
        axes[1,1].set_xlabel('length zero coverage regions')
        axes[1,1].set_ylabel('frequency')
        axes[1,1].set_title('length zero coverage regions')
        axes[1,1].set_xscale('log')
    else:
        # well aeverythin is covered
        axes[1,1].text(0.5,0.5, 'everything is covered', ha='center', va='center', transform=axes[1,1].transAxes)
        axes[1,1].set_title('zero coverage regions')
    
    plt.tight_layout()
    return figure


# now the actual logic starting everything

temp = genome_coverage(file, filename)
temp.savefig(f'{aligner}_{filename}_genome_coverage.png')