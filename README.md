# ACDC-NN

ACDC-NN is a novel antisymmetric neural network to predict proteins free energy changes upon point variations along the amino acid sequence.
The ACDC-NN model was built so that it can be used to make predictions in two different ways: 

1. when both the wild-type and variant structure are available, these are respectively used as direct and inverse inputs so that the network can provide a prediction that, by construction, is perfectly antisymmetric; 

2. when only the wild-type structure is available, as usual, the input for the inverse mutation is created starting from the direct one by inverting the variation encoding but preserving the same structure. 

For further information about the ACDC-NN architecture and properties, please see the related paper (LINK)

## About this repository

Here you can find the instructions to easily install ACDC-NN on your computer using pip (see commands below).
In this version, ACDC-NN was trained using all datasets available in the literature without correcting for sequence similarity.
In case you want to replicate our paper results you will find a jupyter notebook inside the 'results_replication' folder.
There ACDC-NN was trained using a 10-fold cross-validation taking into account sequence similarity to avoid overfitting.

## Installation

We recommend using pip:
```
pip install acdc-nn
```

Requirements:
<table>
  <tr><th>Requirement</th><th>Minimum tested version</th></tr>
  <tr><td>python</td><td>3.6</td></tr>
  <tr><td>tensorflow</td><td>2.3.1</td></tr>
  <tr><td>Biopython</td><td>1.78</td></tr>
  <tr><td>numpy</td><td>1.19.5</td></tr>
  <tr><td>pandas</td><td>1.1.5</td></tr>
  <tr><td>silence_tensorflow</td><td>1.1.1</td></tr>
</table>

## Usage
To predict the change of the folding free energy (DDG) due to a point mutation in a protein sequence, ACDC-NN needs both evolutionary and structural information about the protein itself. The structural information is from a PDB file. The evolutionary information is from a profile file, simple tab-separated table of the frequencies of each residue in each position in homologous proteins. Positive DDG are destabilizing.

When information is available only for the wild-type protein, the predictor can be run as:
```
acdc-nn single MUT PROF PDB CHAIN
```
where MUT is the point mutation, PROF and PDB are the paths to the profile and PDB files, and CHAIN is the PDB chain where the mutation occurs. MUT is in the form XNY where X is the wild-type residue, N is the position of the mutation, and Y is the mutated residue. X and Y are given as a one-letter amino acid code and N is 1-based and referred to the PDB numbering of the relevant chain, and not the position in the sequence. Both PDB and profile files are automatically decompressed when they have a ".gz" extension.

When information is available also for the mutated protein, a better prediction can be got as:
```
> acdc-nn inverse MUT WT-PROF WT-PDB WT-CHAIN INV_MUT MT-PROF MT-PDB MT-CHAIN 
```

To predict more than a few mutations, we provide a batch mode:
```
> acdc-nn batch MUT-TABLE
```
where MUT-TABLE is the path to a table with a row for each mutation to be predicted. For mutations where only the wild-type protein data is available, the row format is:
```
MUT PROF PDB CHAIN
```
For mutations where also the mutated protein data is available, the row format is:
```
MUT WT-PROF WT-PDB WT-CHAIN INV_MUT MT-PROF MT-PDB MT-CHAIN
```

## Examples
Single mutation:
```
> acdc-nn single Q104H tests/profiles/2ocjA.prof tests/structures/2ocj.pdb A
0.15008962
```
Single mutation with the structure of the mutated protein
```
> acdc-nn inverse V51I tests/profiles/1bsaA.prof tests/structures/1bsa.pdb A I51V tests/profiles/1bniA.prof tests/structures/1bni.pdb A 
0.48577148
> acdc-nn inverse I51V tests/profiles/1bniA.prof tests/structures/1bni.pdb A V51I tests/profiles/1bsaA.prof tests/structures/1bsa.pdb A
-0.48577148
```

NB: In the above example we have specifically chosen two homologous proteins that have similar structure.

## Cose ancora da fare:
- scegliere una licenza: GPL / MIT?
- spiegare come produrre i profili, magari mettere uno script
- controllare che gli esempio coincidano con tests/test.sh
