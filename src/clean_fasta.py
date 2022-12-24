#!/usr/bin/env python
#
# Clean fasta files and replace with appropriate end of text token
# '<|endoftext|>'
#

import os
from tqdm import tqdm

################################################################################
## Build Output Scaffolding
################################################################################

os.makedirs('./calc', exist_ok=True)

################################################################################
## Load Data
################################################################################

sets = [
    ('./data/sars_spike_train.fasta', './calc/train.txt', 1000),
    ('./data/sars_spike_test.fasta', './calc/test.txt', 100)
]

for in_fn, out_fn, n in sets:
    n_proc = 0
    with open(out_fn, 'w') as outfile:
        with open(in_fn) as infile:
            for line in tqdm(infile):
                if line.startswith('>'):
                    n_proc += 1
                    line = '<|endoftext|>\n'
                if n_proc > n:
                    break
                print(line, end='', file=outfile)


print('All done!')