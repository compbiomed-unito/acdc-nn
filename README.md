# ACDC-NN

ACDC-NN is a novel antysymmetric neural network to predict proteins free energy changes upon point variations along the aminoacid sequence.
The  ACDC-NN  model  was  built  so that it can be used to make predictions in two different ways: 

1) when both the wild-type and variant structure are available, these are respectively used as direct and inverse  inputs  so  that  the  network  can  provide a  prediction  that,  by  construction,  is  perfectly antisymmetric ; 

2) when only the wild-type structure is available, as usual,  the  inverse  input  is  created  starting  from the direct one by inverting the variation encoding but preserving the same structure. 

For further information about the ACDC-NN architecture and properties, please see the related paper (LINK)

## About this repository

Here you can find the istructions to easily install ACDC-NN on your computer using pip (see commands below).
In this version ACDC-NN was trained using all datasets available in the literature without correcting for sequence similarity.
In case you want to replicate our paper results you will find a jupyter notebook inside the 'results_replication' folder.
There ACDC-NN was trained using a 10-fold cross validation taking into account sequence similarity to avoid overfitting.


## Installation
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
ACDC-NN expects the usage of pdb numbering, not the sequence one.

## Examples
Single mutation:
```
> acdc-nn single Q104H ./profiles/2ocjA.prof ./pdb/2ocj.pdb A
0.15008962
```
Single mutation with the structure of the mutated protein
```
> acdc-nn inverse V51I ./profiles/1bsaA.prof ./pdb/1bsa.pdb A I51V ./profiles/1bniA.prof ./pdb/1bni.pdb A 
0.48577148
> acdc-nn inverse I51V ./profiles/1bniA.prof ./pdb/1bni.pdb A V51I ./profiles/1bsaA.prof ./pdb/1bsa.pdb A
-0.48577148
```

NB: In the above example we have specifically chosen two homologous proteins that have similar structure.

## Cose ancora da fare:
- scegliere una licenza: GPL / MIT?
- spiegare come produrre i profili, magari mettere uno script
- sistemare il readme

- controllare che gli esempio coincidano con tests/test.sh
