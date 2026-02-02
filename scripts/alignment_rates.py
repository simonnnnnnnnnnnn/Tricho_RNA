import matplotlib.pyplot as plt
import json
import sys
import pandas as pd
import seaborn as sns

samples = int(sys.argv[1])
files = sys.argv[2:]
kallisto_files = files[0:samples]
salmon_files = files[samples:2*samples]
bowtie_files = files[2*samples:3*samples]
star_files = files[3*samples:4*samples]

kallisto_list = []
salmon_list = []
bowtie_list = []
star_list = []

def read_kallisto(infile):
    with open(infile, "r") as file:
        data = json.load(file)
        #print(data)
    return data.get('p_pseudoaligned', 0)

def read_salmon(infile):
    with open(infile, "r") as file:
        data = json.load(file)
        #print(data)
    return data.get('percent_mapped', 0)

def read_picard(infile):
    with open(infile, "r") as file:
        for line in file:
            if line.startswith("#") or not line.strip():
                continue # blank lines and comments are not super interesting :D

            if line.startswith("CATEGORY"):
                headers = line.strip().split("\t")
                continue

            if line.startswith("PAIR"):
                values = line.strip().split("\t")
                break # found what i need, looping further is no use, that may not be clean...but it does the job

    for i in range(len(headers)):
        if headers[i] == "PCT_PF_READS_ALIGNED":
            return (float(values[i]) * 100)# get whatever number is there, written as a fraction so times 100
        
    raise ValueError("no percentage found")# thats better than returning something odd and using this as an error signal

for i in range(len(kallisto_files)):
    temp = read_kallisto(kallisto_files[i])
    kallisto_list.append(temp)

for j in range(len(salmon_files)):
    temp = read_salmon(salmon_files[i])
    salmon_list.append(temp)

for k in range(len(bowtie_files)):
    temp = read_picard(bowtie_files[i])
    bowtie_list.append(temp)

for l in range(len(star_files)):
    temp = read_picard(star_files[i])
    star_list.append(temp)


kallisto_avg = sum(kallisto_list) / len(kallisto_list)
salmon_avg = sum(salmon_list) / len(salmon_list)
bowtie_avg = sum(bowtie_list) / len(bowtie_list)
star_avg = sum(star_list) / len(star_list)

# and now plot it all --> now i have the mean alignment rate for all (pseudo)aligners
# make a df so seaborn can easily read that
proto = []

proto.append({'Aligner': 'Kallisto', 'Mapping Rate': kallisto_avg})
proto.append({'Aligner': 'Salmon', 'Mapping Rate': salmon_avg})
proto.append({'Aligner': 'Bowtie2', 'Mapping Rate': bowtie_avg})
proto.append({'Aligner': 'STAR', 'Mapping Rate': star_avg})

frame = pd.DataFrame(proto)

sns.set_style('whitegrid')
plt.figure(figsize=(10,6))
sns.barplot(x='Aligner', y='Mapping Rate', data=frame)
plt.ylabel('Mapping Rate in %')
plt.title('Average Mapping Rates')
plt.tight_layout()
plt.savefig('average_mapping_rates.png')