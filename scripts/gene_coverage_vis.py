import matplotlib.pyplot as plt
import sys
import numpy as np
import pandas as pd

# as with the other vis scripts: files come in via CL, aligner as well

aligner = sys.argv[1]
file = sys.argv[2] # the rest of the args
name = sys.argv[3]

# alternatively filelist and namelist

def plot_gene_cov(infile, outname):

    current = pd.read_csv(infile, sep="\t", names = ['chr', 'source', 'feature', 'start', 'end', 'score',
                                                'strand', 'frame', 'attributes', 'bases_covered',
                                                'total_bases', 'fraction_covered', 'mean_coverage'])
    
    # make e column with the gene ids
    current['gene_id'] = current['attributes'].str.extract(r'gene_id "([^"]+)"')# regex....useful but the code looks like someone rolled his head over the keyboard

    # setuo, 1
    
    figure, axes = plt.subplots(2,2, figsize=(15,10))
    figure.suptitle(f'Gene Coverage {outname}', fontsize=16)
    """
    # mean coverage compared to gene length
    axes[0,0].scatter(current['total_bases'], current['mean_coverage'])
    axes[0,0].set_xscale('log')
    axes[0,0].set_yscale('log')
    axes[0,0].set_xlabel('gene length in bp')
    axes[0,0].set_ylabel('mean coverage')
    axes[0,0].set_title('coverage vs gene length')

    # trendline, because why not, gotta learn that anyway...
    log_len = np.log10(current['total_bases'])
    log_cov = np.log10(current['mean_coverage'].replace(0, np.nan))

    my_mask = ~(np.isnan(log_len) | np.isnan(log_cov))
    if my_mask.sum() > 1:
        z = np.polyfit(log_len[my_mask], log_cov[my_mask], 1)
        p = np.poly1d(z)
        axes[0,0].plot(current['total_bases'], 10**p(log_len), "r--", alpha=0.8)"""
#--------------------------------------------------------------------
    # hexbin plot
    prep = ((current['total_bases'] > 0) & (current['mean_coverage'] > 0) & (current['total_bases'].notna()) & (current['mean_coverage'].notna()))
    x = current['total_bases'][prep]
    y = current['mean_coverage'][prep]
    hexb = axes[0,0].hexbin(x, y, gridsize=50, xscale='log', yscale='log', cmap='viridis', mincnt=1)
    axes[0,0].set_ylabel('mean coverage')
    axes[0,0].set_xlabel('gene length in bp')
    axes[0,0].set_title('coverage vs gene length')
    cb = plt.colorbar(hexb, ax=axes[0,0])
    cb.set_label('count in bin')

    # trendline for hexbin plot
    verts = hexb.get_offsets() # get center and counts of all the bins
    hexcounts = hexb.get_array()
    valid_mask = (verts[:,0] > 0) & (verts[:,1] > 0)
    valid_verts = verts[valid_mask]
    valid_hexcounts = hexcounts[valid_mask]
    hex_log_x = np.log10(verts[:,0])
    hex_log_y = np.log10(verts[:,1])# log transform + bin centers

    # now kikck out all the bins with NaN / inf --> get rid of bad entries that will only cause problems
    """
    hexmask = ~(np.isnan(hex_log_x) | np.isnan(hex_log_y) | np.isinf(hex_log_x) | np.isinf(hex_log_y))
    hex_log_x = hex_log_x[hexmask]
    hex_log_y = hex_log_y[hexmask]
    weights = hexcounts[hexmask] # improve"""
    weights = valid_hexcounts

    # fit lin reg
    z = np.polyfit(hex_log_x, hex_log_y, 1, w=weights)
    p = np.poly1d(z)

    # finally fit the line into the plot:
    x_vals = np.logspace(np.log10(current['total_bases'].min()), np.log10(current['total_bases'].max()), 500)
    axes[0,0].plot(x_vals, 10**p(np.log10(x_vals)), "r--", alpha=0.8, label='weighted trendline')
    axes[0,0].legend()
#---------------------------------------------------------------------
    # completeness
    axes[0,1].hist(current['fraction_covered'], bins=30, alpha=0.7, color='green')
    axes[0,1].axvline(current['fraction_covered'].mean(), color='red', linestyle='--', label=f'mean: {current["fraction_covered"].mean():.2f}')
    axes[0,1].set_xlabel('fraction of gene covered')
    axes[0,1].set_ylabel('number of genes')
    axes[0,1].set_title('gene coverage completeness')
    axes[0,1].legend()

    # best covered genes
    top = current.drop_duplicates(subset='gene_id').nlargest(20, 'mean_coverage')
    y_pos = np.arange(len(top))
    axes[1,0].barh(y_pos, top['mean_coverage'], color='teal')#color = magenta
    axes[1,0].set_yticks(y_pos)
    axes[1,0].set_yticklabels(top['gene_id'], fontsize=8)
    axes[1,0].set_xlabel('mean coverage')
    axes[1,0].set_ylabel('20 best covered genes')

    # dependence coverage to chromosome position --> that could be interesting
    current['midpoint'] = (current['start'] + current['end']) / 2
    for a,b in enumerate(current['chr'].unique()[:10]): #get the first 10 chromosomes again
        chrom_genes = current[current['chr'] == b]
        axes[1,1].scatter(chrom_genes['midpoint'], chrom_genes['mean_coverage'], alpha=0.6, label=b, s=15)

    axes[1,1].set_xlabel('position on chromosome')
    axes[1,1].set_ylabel('mean coverage')
    axes[1,1].set_title('gene coverage along chromosomes')
    axes[1,1].legend()
    axes[1,1].set_yscale('log')

    plt.tight_layout()
    return figure

def make_average(*infiles):
    pass



# now finally call everything and run it
temp = plot_gene_cov(file, name)
temp.savefig(f'{aligner}_{name}_gene_coverage.png')
