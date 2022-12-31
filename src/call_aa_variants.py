#!/usr/bin/env python
#
# Load generated sequences and call amino acid changes from reference
#

import transformers as tfs
import os
import sys
from functools import reduce
from Bio import Align
import numpy as np

################################################################################
## Build Output Scaffolding
################################################################################

os.makedirs('./calc', exist_ok=True)

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

with open('./calc/generated_rbd.txt') as f:
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

print('===== Variants in the Test Set that were predicted =====')
new_subs = list(all_substitutions_test.difference(all_substitutions_train))
predicted = list(map(lambda x: x in all_substitutions, new_subs))
np.array(new_subs)[predicted]

print('===== Novel Variants Predicted =====')

# novel variants predicted by the algorithm
all_substitutions.difference(
    all_substitutions_test
).difference(
    all_substitutions_train
)
