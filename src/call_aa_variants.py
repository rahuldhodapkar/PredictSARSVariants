#!/usr/bin/env python
#
# Load generated sequences and call amino acid changes from reference
#

import os
import sys
from functools import reduce
from Bio import Align
import numpy as np
import pandas as pd
import plotnine as pn
import seaborn as sns
import matplotlib.pyplot as plt

################################################################################
## Build Output Scaffolding
################################################################################

os.makedirs('./calc', exist_ok=True)
os.makedirs('./fig', exist_ok=True)

################################################################################
## Load Reference Sequence
################################################################################

with open('./data/sars_surface_glycoprotein.fasta') as f:
    metadata = f.readline().rstrip()
    refseq = ''.join([line.replace('\n', '') for line in f])

################################################################################
## Align Generated Sequences
################################################################################

aligner = Align.PairwiseAligner()
aligner.substitution_matrix = Align.substitution_matrices.load('BLOSUM62')
aligner.open_gap_score = -2

REFSEQ_START = 426
REFSEQ_STOP = 515

LINEAR_PATH = ((0,0), (REFSEQ_STOP-REFSEQ_START, REFSEQ_STOP-REFSEQ_START))

def check_alignment_substitutions(x):
    subs = filter(lambda x: x != '', 
        map(lambda i: (''.join([x.target[i], str(REFSEQ_START + i + 1), x.query[i]])
            if x.target[i] != x.query[i]
            else ''), range(len(x.query))))
    return list(subs)


all_substitutions = set()

with open('./calc/generated_rbd_small.txt') as f:
    for line in f:
        s = line.rstrip()
        alignments = aligner.align(refseq[REFSEQ_START:REFSEQ_STOP], s)
        if (len(alignments) == 1) and alignments[0].path == LINEAR_PATH:
            all_substitutions.update(
                check_alignment_substitutions(alignments[0]))


# check all alignments from training set

all_substitutions_train = set()

with open('./data/sars_spike_train.fasta') as f:
    s = ''
    for line in f:
        if line.startswith('>') and s != '':
            try:
                alignments = aligner.align(
                    refseq[REFSEQ_START:REFSEQ_STOP], s[REFSEQ_START:REFSEQ_STOP])
            except ValueError:
                s = ''
                continue
            if (len(alignments) == 1) and alignments[0].path == LINEAR_PATH:
                all_substitutions_train.update(
                    check_alignment_substitutions(alignments[0]))
            s = ''
            continue
        s = s + line.rstrip()



# check all alignments from test set

all_substitutions_test = set()

with open('./data/sars_spike_test.fasta') as f:
    s = ''
    for line in f:
        if line.startswith('>') and s != '':
            try:
                alignments = aligner.align(
                    refseq[REFSEQ_START:REFSEQ_STOP], s[REFSEQ_START:REFSEQ_STOP])
            except ValueError:
                s = ''
                continue
            if (len(alignments) == 1) and alignments[0].path == LINEAR_PATH:
                all_substitutions_test.update(
                    check_alignment_substitutions(alignments[0]))
            s = ''
            continue
        s = s + line.rstrip()

print('===== Variants in the Training Set =====')
print(all_substitutions_train)

print('===== Variants in the Test Set =====')
print(all_substitutions_test)

print('===== Variants in the Predicted Set =====')
print(all_substitutions)

print('===== New Variants in the Test Set that were predicted =====')
new_subs = list(all_substitutions_test.difference(all_substitutions_train))
predicted = list(map(lambda x: x in all_substitutions, new_subs))
np.array(new_subs)[predicted]

print('===== All Variants in the Test or Train Set that were predicted =====')
new_subs = list(all_substitutions_test.union(all_substitutions_train))
predicted = list(map(lambda x: x in all_substitutions, new_subs))
np.array(new_subs)[predicted]

print('===== Novel Variants Predicted =====')
# novel variants predicted by the algorithm
all_substitutions.difference(
    all_substitutions_test
).difference(
    all_substitutions_train
)

################################################################################
## Generate Plots
################################################################################

####
# Variant location along primary amino acid sequence density plot
####

subs_dict = {
    'train': all_substitutions_train,
    'test': all_substitutions_test,
#    'only_test': all_substitutions_test.difference(all_substitutions_train),
#    'only_train': all_substitutions_train.difference(all_substitutions_test),
    'predict': all_substitutions
}

df_dicts = []
for k,v in subs_dict.items():
    df_dicts += [{'group':k, 'value':int(x[1:-1])} for x in v]

plot_df = pd.DataFrame(df_dicts)

sns.displot(data=plot_df,
            x='value',
            hue='group',
            kind='kde',
            common_norm=False,
            fill=True,
            palette=sns.color_palette('colorblind'),
            bw_adjust=0.2)


####
# Variant Pie Chart
####


substitutions_df = pd.DataFrame({
    'type': 
        ['PREDICTED'] * len(all_substitutions_test.intersection(
                            all_substitutions)) +
        ['NOT_PREDICTED'] * len(all_substitutions_test.difference(
                            all_substitutions))
})
data = substitutions_df.groupby("type")['type'].count()
data.plot.pie(autopct="%.1f%%", explode=[0.05]*2, colors=sns.color_palette('colorblind'))



substitutions_df = pd.DataFrame({
    'type': 
        ['PREDICTED'] * len(all_substitutions_train.intersection(
                            all_substitutions)) +
        ['NOT_PREDICTED'] * len(all_substitutions_train.difference(
                            all_substitutions))
})
data = substitutions_df.groupby("type")['type'].count()
data.plot.pie(autopct="%.1f%%", explode=[0.05]*2, colors=sns.color_palette('colorblind'))



substitutions_df = pd.DataFrame({
    'type': 
        ['IN_TEST_ONLY'] * len(all_substitutions.intersection(
                            all_substitutions_test.difference(
                            all_substitutions_train))) +
        ['IN_TRAIN_ONLY'] * len(all_substitutions.intersection(
                            all_substitutions_train.difference(
                            all_substitutions_test))) +
        ['IN_TRAIN_AND_TEST'] * len(all_substitutions.intersection(
                            all_substitutions_train.intersection(
                            all_substitutions_test))) +
        ['NOVEL'] * len(all_substitutions.difference(
                            all_substitutions_train.union(
                            all_substitutions_test)))
})
data = substitutions_df.groupby("type")['type'].count()
data.plot.pie(autopct="%.1f%%", explode=[0.05]*4,
             colors=sns.color_palette('colorblind'))

################################################################################
## Generate BLOSUM plots
################################################################################

blosum = Align.substitution_matrices.load("BLOSUM80")

samples = [np.random.choice(
    np.array(range(len(blosum62.alphabet)), dtype=int), 2, replace=False)
    for x in range(len(all_substitutions))]


subs_dict = {
    'train': all_substitutions_train,
    'test': all_substitutions_test,
    'predict': all_substitutions,
    'random': [''.join(map(lambda x: blosum.alphabet[x], x)) for x in samples]
}

df_dicts = []
for k,v in subs_dict.items():
    df_dicts += [{'group':k, 'value':blosum[x[0],x[-1]]} for x in v]

plot_df = pd.DataFrame(df_dicts)

sns.boxplot(data=plot_df,
            x='group',
            y='value',
            palette=sns.color_palette('colorblind'),
            showfliers=False)



