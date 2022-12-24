#!/usr/bin/env python
#
# Clean fasta files and replace with appropriate end of text token
# '<|endoftext|>'
#

import os
from tqdm import tqdm
import random

random.seed(42)

################################################################################
## Build Output Scaffolding
################################################################################

os.makedirs('./calc', exist_ok=True)

################################################################################
## Load Data
################################################################################

sets = [
    ('./data/sars_spike_train.fasta', './calc/train.txt', 10000),
    ('./data/sars_spike_test.fasta', './calc/test.txt', 1000)
]

for in_fn, out_fn, n in sets:
    # first count and index the entries in the file
    n_entries = 0
    with open(in_fn) as infile:
        for line in tqdm(infile):
            if line.startswith('>'):
                n_entries += 1

    # then randomly select the number of entries to print
    valid_ixs = set(random.sample(range(n_entries), n))

    n_proc = -1                 # required for zero-indexing of valid_ixs
    with open(out_fn, 'w') as outfile:
        with open(in_fn) as infile:
            for line in tqdm(infile):
                if line.startswith('>'):
                    n_proc += 1
                    line = '<|endoftext|>\n'
                if n_proc in valid_ixs:
                    print(line, end='', file=outfile)

print('All done!')