import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt
import seaborn as sns

info = sys.argv[1]
files = sys.argv[2:] # all tpm files

# plan:
# group files by mapfile[name] --> compare same samples
# then match each of the 4 files per group (4 aligners so 4 files)
# then compare tpm vals per gene

mapfile = pd.read_csv(info, sep=",", header=0)
namepatterns = mapfile['name']
"""
for i in range(len(namepatterns)):
    filelist = []
    for j in range(len(files)):
        if namepatterns[i] in files[j]:
            filelist.append(files[j])# sort by name --> therefore making sure same condition and replicate are used
    # attach namepattern as well for naming
    filelist.append(namepatterns[i])

    # now plot the differences in TPM values between all 4 aligners
    plot_tpm_diffs(filelist, namepatterns[i])
"""

def make_tpm_dict_kallisto(infile):
    temp = pd.read_csv(infile, sep="\t", header=0)
    temp = temp.rename(columns={"target_id": "id", "tpm": "TPM"})
    res_dict = dict(zip(temp['id'], temp['TPM']))
    return res_dict

def make_tpm_dict_salmon(infile):
    temp = pd.read_csv(infile, sep="\t", header=0)
    temp = temp.rename(columns={"Name": "id"})
    res_dict = dict(zip(temp['id'], temp['TPM']))
    return res_dict

def make_tpm_dict_bowtie2(infile):
    temp = pd.read_csv(infile, sep="\t", header=0)
    temp = temp.rename(columns={"Geneid": "id"})
    res_dict = dict(zip(temp['id'], temp['TPM']))
    return res_dict

def make_tpm_dict_star(infile):
    temp = pd.read_csv(infile, sep="\t", header=0)
    temp = temp.rename(columns={"Geneid": "id"})
    res_dict = dict(zip(temp['id'], temp['TPM']))
    return res_dict

def plot_tpm_diffs(tpm_files, name):
    # map to aligners
    for i in range(len(tpm_files)):
        if "abundance" in tpm_files[i]:
            kallisto = tpm_files[i]
        elif "quant" in tpm_files[i]:
            salmon = tpm_files[i]
        elif "Aligned" in tpm_files[i]:
            star = tpm_files[i]
        else:
            bowtie2 = tpm_files[i] # that could probably be a tad more efficient, but with the astronomuc runtimes for the alignments, saving a second here really doesnt matter...

    # get dicts amtching ids to tpm vals
    kallisto_tpms = make_tpm_dict_kallisto(kallisto)
    salmon_tpms = make_tpm_dict_salmon(salmon)
    bowtie2_tpms = make_tpm_dict_bowtie2(bowtie2)
    star_tpms = make_tpm_dict_star(star)

    # make dfs from the dicts
    tpm_dicts = {'kallisto': kallisto_tpms,
                 'salmon': salmon_tpms,
                 'bowtie2': bowtie2_tpms,
                 'star': star_tpms}
    
    tpm_frame = pd.DataFrame(tpm_dicts)
    tpm_frame.index.name = 'id'
    tpm_frame.reset_index(inplace=True)

    # get the average tpm for all
    tpm_frame['avg'] = tpm_frame[['kallisto', 'salmon', 'bowtie2', 'star']].mean(axis=1)

    # now a nice line graph/diagram
    longframe = pd.melt(tpm_frame, id_vars='id', value_vars=['kallisto', 'salmon', 'bowtie2', 'star', 'avg'], var_name='aligner', value_name='TPM')

    plt.figure(figsize=(12,6))
    sns.lineplot(data=longframe, x='id', y='TPM', hue='aligner')

    plt.title(f'TPM comparison between aligner, {name}')
    plt.xlabel('Gene ID')
    plt.ylabel('TPM Values')
    plt.xticks(rotation=45)
    plt.grid(True)
    plt.tight_layout()
    plt.legend(title='aligner')
    plt.savefig(f'TPM comparison {name}.png')


for i in range(len(namepatterns)):
    filelist = []
    for j in range(len(files)):
        if namepatterns[i] in files[j]:
            filelist.append(files[j])# sort by name --> therefore making sure same condition and replicate are used
    # attach namepattern as well for naming
    #filelist.append(namepatterns[i])

    # now plot the differences in TPM values between all 4 aligners
    plot_tpm_diffs(filelist, namepatterns[i])


