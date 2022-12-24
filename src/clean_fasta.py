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

# transfer text file and process
with open('./calc/sars_spike_test.fasta', 'w') as outfile:
    with open('./data/sars_spike_test.fasta') as infile:
        for line in tqdm(infile):
            if line.startswith('>'):
                print('<|endoftext|>', file=outfile)
            else:
                print(line, end='', file=outfile)


# transfer text file and process
with open('./calc/sars_spike_train.fasta', 'w') as outfile:
    with open('./data/sars_spike_train.fasta') as infile:
        for line in tqdm(infile):
            if line.startswith('>'):
                print('<|endoftext|>', file=outfile)
            else:
                print(line, end='', file=outfile)


print('All done!')