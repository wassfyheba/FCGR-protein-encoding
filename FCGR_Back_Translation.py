# -*- coding: utf-8 -*-
"""
Created on Tue May  2 15:55:17 2023

@author: D E L L
"""


"""code to remove header and reverse translation for multipe protein sequences"""


from itertools import product
import collections
from collections import OrderedDict
from matplotlib import pyplot as plt
from matplotlib import cm
import pylab
import math
import os
import re
import pickle 

table = {
        'ATT':'I', 'ATG':'M','ACT':'T','AAC':'N',
        'AAG':'K', 'CTA':'L', 'CCA':'P','CAC':'H',
        'CAG':'Q','CGA':'R', 'GTG':'V', 'GCT':'A',
        'GAC':'D','GAG':'E','GGT':'G', 'TCA':'S',
        'TTC':'F', 'TAC':'Y','TGC':'C', 'TGG':'W'}

BASE_FILE = "D:/phd thesis/mydata/modified GDS/gds_dataset"
protein_file="{0}/ClassA_Cannabinoid_Cannabinoid.txt".format(BASE_FILE)



ALL_DATA = []
with open(protein_file, 'r') as fin:
    DATA = ''
    for line in fin:
        if line.startswith('>'):
            ALL_DATA.append(DATA)
            DATA = ''
        else:
            DATA += line.strip()
ALL_DATA.append(DATA)
ROW_DATA=ALL_DATA[1:]


ambig=['X','x','Z','z','B','b','U','u','_','O','o','J','j',' ']

deleted_ambig=[]
for v,line in enumerate(ROW_DATA):
    for resid in line:
        if resid in ambig:
            del(ROW_DATA[v])
            deleted_ambig.append(v)
        # if resid=='X' or resid=='Z':
        #     del(ROW_DATA[v])
        #     deleted_ambig.append(v)
        # if resid=='B' or resid=='U':
        #     del(ROW_DATA[v])
        #     deleted_ambig.append(v)
    
# remone sequences longer than 10000

# deleted_long=[]
# filterd_seq=[]
# for num, line in enumerate(ROW_DATA):
#     if len(line)<=1000:
#         filterd_seq.append(line)
#     else:
#         deleted_long.append(num)
        
# Back translation from protein to DNA

degenDict = dict()
for v,l in table.items():
    if l in degenDict:
        degenDict[l].append(v)
    else:
        degenDict[l]=[]
        degenDict[l].append(v)
        
#del(ROW_DATA[69])
#del(ROW_DATA[69])

# del(ROW_DATA[63])
# del(ROW_DATA[171])
#del(ROW_DATA[388])
# del(ROW_DATA[388])
# del(ROW_DATA[52])
# del(ROW_DATA[1])


# del(ROW_DATA[171])
# del(ROW_DATA[202])
# del(ROW_DATA[202])
# del(ROW_DATA[480])
# del(ROW_DATA[618])
# del(ROW_DATA[671])
# del(ROW_DATA[807])


reverse_prot=[]

for line in ROW_DATA:
    nucs = [degenDict[resid] for resid in line]
    for degenNuc in product(*nucs):
        data=''.join(degenNuc)
        reverse_prot.append(data)
        
# # print(reverse_prot)

# # """FCGR CODE"""

def count_kmers(sequence, k):
    d = collections.defaultdict(int)
    for i in range(len(sequence)-(k-1)):
        d[sequence[i:i+k]] +=1
    return d
    
def probabilities(kmer_count, k):
    probabilities = collections.defaultdict(float)
    N = len(data)
    for key, value in kmer_count.items():
        probabilities[key] = float(value) / (N - k + 1)
    return probabilities

def chaos_game_representation(probabilities, k):
    array_size = int(math.sqrt(4**k))
    chaos = []
    for i in range(array_size):
        chaos.append([0]*array_size)
 
    maxx = array_size
    maxy = array_size
    posx = 1
    posy = 1
    
    for key, value in probabilities.items():
        for char in key:
            if char == "G":
                posx=int(posx + maxx / 2)
            elif char == "C":
                posy=int(posy + maxy / 2)
            elif char == "T":
                posx =int(posx + maxx / 2)
                posy=int(posy + maxy / 2)
            maxx = maxx / 2
            maxy /= 2
        chaos[posy-1][posx-1] = value
        maxx = array_size
        maxy = array_size
        posx = 1
        posy = 1
    # print(chaos)
    return chaos


for i, line in enumerate(reverse_prot):
    f3= (count_kmers(line, 9))
    f3_prob = probabilities(f3, 9)
    chaos_k3=(chaos_game_representation(f3_prob, 9))
    # pylab.title('FCGR for sequence no.'+str (i))
    pylab.imshow(chaos_k3,interpolation='nearest',cmap=cm.gray_r)
    plt.axis('off')
    pylab.show()
    
