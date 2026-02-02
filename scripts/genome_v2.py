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

data = pd.read_csv(file, sep="\t", names = ['chr', 'start', 'end', 'coverage'])

# some additional columns that will be helpful later on

data['length'] = data['end'] - data['start']
data['midpoint'] = (data['start'] + data['end'])/2

figure, axes = plt.subplots(2,2, figsize=(15,10))
figure.suptitle('Genome coverage analysis', fontsize=16, fontweight='bold')

# first the coverage along the whoile geome
# this will be a step plot
ax1 = axes[0,0]
positions = []
cover = []

for sth,row in data.iterrows():
    positions.extend([row['start'], row['end']])
    cover.extend([row['coverage'], row['coverage']])

ax1.plot(positions, cover, linewidth=2, color='seagreen')
ax1.fill_between(positions, cover, alpha=0.3, color='seagreen')
ax1.set_xlabel('genomic postion in bp')
ax1.set_ylabel('coverage depth')
ax1.set_title('coverage along genome')
ax1.grid(True, alpha=0.3)


# now a coverage distribution histogram, same idea as woth the gene coverage hist
ax2 = axes[0,1]
sns.histplot(data=data, x='coverage', bins=range(1, data['coverage'].max()+2), kde=True, ax=ax2, color='coral')
ax2.set_xlabel('coverage depth')
ax2.set_ylabel('frequency')
ax2.set_title('coverage distribution')
ax2.grid(True, alpha=0.3)

# something i havent done before: plot interval lenght against a coverage scatter plot
ax3 = axes[1,0]
scatter = ax3.scatter(data['length'], data['coverage'], c=data['coverage'], cmap='viridis', alpha=0.7, s=50)
ax3.set_xlabel('interval length in bp')
ax3.set_ylabel('coverage depth')
ax3.set_title('interval length vs coverage')
ax3.grid(True, alpha=0.3)
plt.colorbar(scatter, ax=ax3, label='coverage')


# last but not least: a boxplot showing the coverage

ax4 = axes[1,1]
data['coverage_bin'] = pd.cut(data['coverage'], bins=[0,2,5,10, float('inf')], labels=['low (1-2)', 'medium (3-5)', 'high (6-10)', 'very high (>10)'])
sns.boxplot(data=data, x='coverage_bin', y='length', ax=ax4, palette='Set2')
ax4.set_xlabel('coverage categories')
ax4.set_ylabel('interval length in bp')
ax4.set_title('interval length by coverage category')
ax4.tick_params(axis='x', rotation=45)

plt.tight_layout()
plt.savefig(f'genome_coverage_{name}_{aligner}.png')

