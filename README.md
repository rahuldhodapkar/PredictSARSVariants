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

## Evaluation
Distance between predicted sequences and true variant spike proteins were
computed by using the string edit distance between the amino acid sequences,
as compared to chance. 

As the loaded ProtGPT2 model was pretrained on the
UniRef50 (version 2021_04) dataset, it cannot have contained sequencing
data that was generated after that date.  Evaluations will be conducted using
SARS-CoV2 sequences generated on or after May 2021.

## Orchestration
Code to run these experiments is organized using a simple `Makefile`.
All software is designed to be run within a `conda` environment, requirements
file to be `pip` installed once environment is initialized.

    conda create -n covid_pred python=3.9

Then:

    conda activate covid_pred
    pip install -r requirements.txt