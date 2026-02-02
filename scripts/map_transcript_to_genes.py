import pandas as pd
import sys


filepath = sys.argv[1]
with open(filepath, "r") as file:
    good_lines = []
    for line in file:
        if line.startswith("#"):
            pass
        else:
            good_lines.append(line.strip().split("\t"))# these are the good lines i want to keep

tx2gene = []

# one thing to note: the third field here, so with index 2 has the info if its an mRNA or something else
# the last field (index = 8), the one with the description is split to extract the gene id
for element in good_lines:
    if element[2] == "mRNA":
        temp = element[8].split(";")
        transcript = None
        gene = None
        for entries in temp:
            if "ID=transcript" in entries:
                transcript = entries.split(":")[1]
            if "Parent=gene" in entries:
                gene = entries.split(":")[1]
        
        if transcript and gene:#only put stuff into tx2gene when both entries exist --> no defective entries
            tx2gene.append([transcript, gene])
            #key = kv_pair.split("=")[0]
            #values = kv_pair.split("=")[1]

# now pandas is finally used
frame_tx2gene = pd.DataFrame(tx2gene, columns=["transcriptID", "geneID"])
frame_tx2gene.to_csv("tx2gene.csv", index = False)