#!/usr/bin/python

import re
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint as IP
import random
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Levenshtein import distance as levenshtein_distance

dist_file = ("RAxML_distances.T4_diag")
ind_file = ("Leiden_1993_final_SW_het_data.tsv")
fasta_file = ("Haplotypes_MACSE_HVR_AA.fas")
fasta_file_full = ("superalgn_aa.fasta")


dist_file_handle = open(dist_file)
ind_file_handle = open(ind_file)
ind_file_handle2 = open(ind_file)

dist_hash = dict()
seqNoHVR_hash = dict()
seq_hash = dict()
len_hash = dict()
HVR_hash = dict()
Iso_hash = dict()

mer_dict=dict()

#https://www.pnas.org/doi/10.1073/pnas.2018234118
kd = { 'A': 1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C': 2.5,
       'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I': 4.5,
       'L': 3.8,'K':-3.9,'M': 1.9,'F': 2.8,'P':-1.6,
       'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V': 4.2 }

## for HVR len diff etc
for record in SeqIO.parse(fasta_file, "fasta"):
    AA_noGap = re.sub('-', '', str(record.seq))
    #print(record.id)
    len_hash[record.id] = len(AA_noGap)
    HVR_hash[record.id]= AA_noGap
    values = []
    for residue in AA_noGap:
        values.append(kd[residue])

    #print(AA_noGap.count('N'))
    Iso_hash[record.id] = sum(values)
    #print(sum(values))

## for levenshtein dist
for record in SeqIO.parse(fasta_file_full, "fasta"):
    AA_HVR = re.sub(r'[0-9]+', '', str(record.seq))
    #print(AA_HVR)
    seqNoHVR_hash[record.id] = SeqRecord(Seq(AA_HVR), id=record.id)
    seq_hash[record.id] = str(record.seq)
    length= len(str(record.seq))
    #values = []
    #for residue in AA_noGap:
    #    values.append(kd[residue])

    #print(AA_noGap.count('N'))
    #Iso_hash[record.id] = sum(values)
    #print(sum(values))


for line in dist_file_handle:
    line_strip = line.rstrip()
    values = line_strip.split()
    dist_hash[tuple(values[0:2])] = values[2]
    dist_hash[tuple(values[2:0])] = values[2]
    #print(tuple(values[0:2]),values[2])


female_list = list()
male_list = list()
male_pairs = list()

male_dat = list()
female_dat = list()


out_file_handle = open("Leiden_1993_final_dists.tsv",'w')

for line in ind_file_handle:
    line_strip = line.rstrip()
    values = line_strip.split()
    #print(values)
    if values[1] == 'male':

        HVR_len_diff = abs(len_hash[values[2]]-len_hash[values[3]])

        HVR_N_diff = abs(Iso_hash[values[2]] - Iso_hash[values[3]])

        align = levenshtein_distance(seqNoHVR_hash[values[2]].seq,seqNoHVR_hash[values[3]].seq)

        HVR_dist = levenshtein_distance(HVR_hash[values[2]],HVR_hash[values[3]])
        #print(seqNoHVR_hash[values[2]].seq,seqNoHVR_hash[values[3]].seq)
        #print(align)
        #print("\t".join(values)+"\t"+dist_hash[tuple(values[2:4])]+"\n")
        out_file_handle.write("\t".join(values)+"\t"+dist_hash[tuple(values[2:4])] +"\t"+ str(HVR_len_diff) +"\t"+ str(HVR_N_diff)+"\t"+ str(align) + "\t" + str(HVR_dist) +"\n")
        male_list.append(values[2])
        male_list.append(values[3])
        male_pairs.append(tuple(values[2:4]))
        male_dat.append(dist_hash[tuple(values[2:4])])
    elif values[1] == 'fem':
        #print(values[2])
        female_list.append(values[2])

#print(male_pairs)
window = 1
for pos in range(0,length-window):
    mer_dict[pos] = [0,0]
    for pair in male_pairs:
     #print(seq_hash[values[2]][pos:pos+window],seq_hash[values[3]][pos:pos+window])
        if seq_hash[pair[0]][pos:pos+window] == seq_hash[pair[1]][pos:pos+window]:
            #print("Hom")
            mer_dict[pos][0]+=1
        else:
            mer_dict[pos][1]+=1
            #print("Het")


for i in mer_dict:
    print(i, str(mer_dict[i][0]), str(mer_dict[i][1]),sep=",")




#print(male_list)
i=0
while i <= 1000:
#    for keys in male_list:
        #print(keys,dist_hash[keys])
        pseudo_male = random.sample(male_list,2)
        while pseudo_male[0] == pseudo_male[1]:
            pseudo_male = random.sample(male_list, 2)
        female_dat.append(dist_hash[tuple(pseudo_male)])

        HVR_len_diff = abs(len_hash[pseudo_male[0]] - len_hash[pseudo_male[1]])
        HVR_N_diff = abs(Iso_hash[pseudo_male[0]] - Iso_hash[pseudo_male[1]])
        align = levenshtein_distance(seqNoHVR_hash[pseudo_male[0]].seq, seqNoHVR_hash[pseudo_male[1]].seq)
        HVR_dist = levenshtein_distance(HVR_hash[pseudo_male[0]], HVR_hash[pseudo_male[1]])

        out_file_handle.write("dummy"+str(i)+"\tpseu\t" + pseudo_male[0] + "\t" + pseudo_male[1] + "\t" + dist_hash[tuple(pseudo_male)] +"\t"+ str(HVR_len_diff) +"\t"+ str(HVR_N_diff)+ "\t" + str(align) + "\t" + str(HVR_dist) + "\n")
        i+=1
        #print(pseudo_male, dist_hash[tuple(pseudo_male)])



#print(male_dat)
#sns.kdeplot(male_dat)
#sns.kdeplot(female_dat, hist=False)


#plt.show()
