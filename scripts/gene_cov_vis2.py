import matplotlib.pyplot as plt
import sys
import numpy as np
import pandas as pd
import seaborn as sns

plt.style.use('default')
sns.set_palette("husl")
# as with the other vis scripts: files come in via CL, aligner as well

aligner = sys.argv[1]
file = sys.argv[2] # the rest of the args
name = sys.argv[3]

# prepare: this is the method on which all other rely
def load(infile):
    frame = pd.read_csv(infile, sep='\t', header=None, names=['scaffold', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes', 'covered_bases', 'feature_length', 'total_bases', 'coverage_fraction'])

    # now i want to get some additional columns that will be useful later on
    frame['gene_id'] = frame['attributes'].str.extract(r'gene_id "([^"]+)"')
    frame['transcript_id'] = frame['attributes'].str.extract(r'transcript_id "([^"]+)"')
    frame['feature_size'] = frame['end'] - frame['start'] + 1
    frame['coverage_depth'] = frame['covered_bases'] / frame['feature_length']

    return frame

# ok we start with a coverage distribution plot with 4 subplots

def coverage_distribution(df):
    figure, axes = plt.subplots(2,2, figsize=(15,10))

    # coverage fraction --> histogram
    axes[0,0].hist(df['coverage_fraction'], bins=50, alpha=0.7, color='skyblue', edgecolor='black')
    axes[0,0].set_xlabel('coverage fraction')
    axes[0,0].set_ylabel('frequency')
    axes[0,0].set_title('Coverage Fraction Distributions')
    # trednline is always nice
    axes[0,0].axvline(df['coverage_fraction'].mean(), color='red', linestyle='--', label=f'Mean:{df["coverage_fraction"].mean():.1f}')
    axes[0,0].legend()

    # now lets look at the feature type --> the coverage by feature type
    featcov = df.groupby('feature')['coverage_fraction'].mean().reset_index()
    bars = axes[0,1].bar(featcov['feature'], featcov['coverage_fraction'], color='lightcoral', alpha=0.8)
    axes[0,1].set_xlabel('Feature Type')
    axes[0,1].set_ylabel('Mean coverage fraction')
    axes[0,1].set_title('Average coverageby feature type')
    axes[0,1].tick_params(axis='x', rotation=45)

    # that is all noce, but the bars should have some labels --> values, so that it is easier to read
    for bar in bars:
        h = bar.get_height()
        axes[0,1].text(bar.get_x() + bar.get_width()/2., h+0.01, f'{h:.1f}', ha='center', va='bottom')

    # now one of the most interesting ones: feature length vs coverage
    axes[1,0].scatter(df['feature_length'], df['coverage_fraction'], alpha=0.6, s=30, c=df['feature'].astype('category').cat.codes)
    axes[1,0].set_xlabel('feature length in bp (logscaled)')
    axes[1,0].set_ylabel('coverage fraction')
    axes[1,0].set_title('feature length vs coverage')
    axes[1,0].set_xscale('log')

    # coverage depth distribution
    sns.boxplot(data=df, x='feature', y='coverage_fraction', ax=axes[1,1])
    axes[1,1].set_xlabel('Feature Type')
    axes[1,1].set_ylabel('coverage fraction')
    axes[1,1].set_title('ccoverage depth by feature type')
    axes[1,1].tick_params(axis='x', rotation=45)

    plt.tight_layout()
    return figure
    #plt.savefig(f'{name} coverage plots.png')

def coverage_heatmap(df, outname, top=20,):
    

    gene_counts = df['gene_id'].value_counts().head(top)
    top_genes = gene_counts.index.tolist()
    heat_data = df[df['gene_id'].isin(top_genes)].pivot_table(values='coverage_fraction', index='gene_id', columns='feature', aggfunc='mean')

    plt.figure(figsize=(10,8))
    sns.heatmap(heat_data, annot=True, fmt='.3f', cmap='RdYlBu_r', cbar_kws={'label': 'Coverge Fraction'})
    plt.title(f'coverage heatmap of top {top} genes by feature count')
    plt.xlabel('feature type')
    plt.ylabel('gene id')
    plt.xticks(rotation=45)
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.savefig(f'coverage heatmap {outname}')



frame = load(file)
plot1 = coverage_distribution(frame)
plot1.savefig(f'coverage_distribution_{name}.png')
coverage_heatmap(frame, name)