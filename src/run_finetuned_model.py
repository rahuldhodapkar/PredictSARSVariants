#!/usr/bin/env python
#
# Load trained model using `run_clm.py` script.
#

import transformers as tfs
import os
import sys
from tqdm import tqdm

################################################################################
## Build Output Scaffolding
################################################################################

os.makedirs('./calc', exist_ok=True)

################################################################################
## Load Model
################################################################################

generator = tfs.pipeline(
    task='text-generation',
    model=sys.argv[1],
    tokenizer='nferruz/ProtGPT2'
)

################################################################################
## Run Model
################################################################################

ntd_seed = '<|endoftext|>'
rbd_seed = 'DDFTGCVIAW'
rbd_length = 69
rbd_pad_len = 5

n_iters = 1000

with open('./calc/generated_rbd.txt', 'w') as f:
    for i in tqdm(range(n_iters)):
        sequences = generator(
            rbd_seed,
            max_length=len(rbd_seed) + rbd_length + rbd_pad_len,
            do_sample=True,
            top_k=950,
            repetition_penalty=1.2,
            num_return_sequences=100,
            eos_token_id=0)

        for seq in sequences:
            print(seq['generated_text'][:90].replace('\n',''), file=f)

print('All done!')

