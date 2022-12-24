#!/usr/bin/env python
#
# Use a pretrained large language model to:
#
#   (1) generate predicted protein sequences for the SARS-CoV2 
#       receptor binding domain (RBD)
#
#   (2) perform unsupervised clustering of predicted sequences to obtain
#       homology groups.
#
#   (3) compare homology groups to RBD sequences derived from after the time
#       that the pretraining data set was constructed.
#

import transformers as tfs
import jellyfish as jf
import os
import time

################################################################################
## Build Output Scaffolding
################################################################################

os.makedirs('./calc', exist_ok=True)

################################################################################
## Load Data
################################################################################

protgpt2 = tfs.pipeline('text-generation', model="nferruz/ProtGPT2")

################################################################################
## (1) Generate Predicted Sequences
################################################################################

prompt = '<|endoftext|>'.join([
    'MGMTA', 'MGMTAACC',
    'MGMT', 'MGMTACC',
    'MGM'
])

start = time.time()
sequences = protgpt2('M',
    max_length=100, do_sample=True,
    top_k=950, repetition_penalty=1.2,
    num_return_sequences=10, eos_token_id=0)
end = time.time()
print(end - start)

sequences[0]

################################################################################
## (2) Cluster Predicted Sequences
################################################################################




################################################################################
## (3) Compare to Ground Truth Sequences
################################################################################


print('All done!')