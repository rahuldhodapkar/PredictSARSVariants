# Predict SARS-Cov2 Variants

This repository contains an implementation of using pretrained language models
to predict the evolution of the SARS-CoV2 spike protein as variants develop.

## Models
The pretrained ProtGPT2 model was used for this purpose - to predict variations
in the RBD of the SARS-CoV2 spike (S) protein.  

## Data
SARS-CoV2 spike protein sequences were obtained from the NIH Sars-CoV2 Data Hub
accessible at 

    https://www.ncbi.nlm.nih.gov/labs/virus/vssi/

Note that the reference sequence for the surface glycoprotein can be found at:

    https://www.ncbi.nlm.nih.gov/protein/1791269090

## Evaluation
Distance between predicted sequences and true variant spike proteins were
computed by using the string edit distance between the amino acid sequences,
as compared to chance. 

As the loaded ProtGPT2 model was pretrained on the
UniRef50 (version 2021_04) dataset, it cannot have contained sequencing
data that was generated after that date.  Evaluations will be conducted using
SARS-CoV2 sequences generated on or after May 2021.

### MutaBind2
The MutaBind2 tool was used to evaluate the effect of single-residue amino
acid substitutions on viral fitness. (https://lilab.jysw.suda.edu.cn/research/mutabind2/)

Note the following relevant PDB codes:

    7DF4 (Spike-ACE2 complex Cryo-EM solved structure)

Binding of RBD to therapeutic antibodies:

    8D8Q (Tixagevimab + cilgavimab (Evushield) used for COVID PrEP)

Some manual entry of chain identities were required to define
which particular amino acid chains should be included in the binding
affinity simulations.

### SAAMBE-3D
The SAAMBE-3D tool was also set up to evaluate fitness changes after amino
acid point substitution, **however this implementation is unstable and should
not be used.** (http://compbio.clemson.edu/saambe_webserver/)

A standalone `python2` implementation was provided by the SAAMBE authors
and was used for rapid evaluation of the effect of point substitutions
on binding affinities/free energy.

After copying the appropriate files into the unzipped `standaloneCode`
folder, the following invocations can be used to call SAAMBE-3D as of Jan 2023.
For the most up-to-date documentation regarding the tool, please contact
the original authors.

Several notes on the installation requirements: (1) The SAAMBE tool requires
a deprecated version of python, which is not being actively maintained.  This
puts the tool's installation at risk (worth keeping in mind when setting up).
(2) For my setup using anaconda with `python 2.7.18`, `conda-forge` installation
of the `xgboost` dependency was required to prevent segmentation fault on
running the tool.  Further, conda was unable to appropriately set up the
path search with the correct `CONDA_PREFIX`, so this needed to be done manually.

After these additional setup steps, SAAMBE was able to predict changes to the
binding affinities of the relevant complexes with the following commands:

    export PATH=$CONDA_PREFIX/bin:$PATH
    python Mutation_pred.py -i 7df4.pdb -f saambe3d_predicted_substitutions.txt -d 1 >saambe3d_predicted_output.txt
    python Mutation_pred.py -i 8d8q.pdb -f saambe3d_evushield.txt -d 1 >saambe3d_evushield_output.txt


## Orchestration
Code to run these experiments is organized using a simple `Makefile`.
All software is designed to be run within a `conda` environment, requirements
file to be `pip` installed once environment is initialized.

    conda create -n covid_pred python=3.9

Then:

    conda activate covid_pred
    pip install -r requirements.txt

If you are running on the yale Farnam HPC, you will need to additionally run:

    conda install pytorch==1.11.0 torchvision==0.12.0 torchaudio==0.11.0 cudatoolkit=11.3 -c pytorch

And if you would like to finetune the `ProtGPT2` model, you will need to install
huggingface Tranformers from source with:

    pip install git+https://github.com/huggingface/transformers
