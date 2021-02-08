# ACDC NN

Breve descrizione, referenza, modifica di Silvia.

## Installation
```
pip install acdc-nn
```

Requirements:
	- python 3.6.12
	- tensorflow 2.3.1
	- Biopython 1.78
	- numpy 
	- silence_tensorflow 


## Usage


## Examples
Single mutation:
```
> acdc-nn Q104H ./profiles/2ocjA.prof ./pdb/2ocj.pdb A
0.1735
```
Single mutation with the structure of the mutated protein
```
> acdc-nn V51I ./profiles/1bsaA.prof ./pdb/1bsa.pdb A --inverse I51V ./profiles/1bniA.prof ./pdb/1bni.pdb A 
0.6103
```

NB: In the above example we have specifically chosen two homologous proteins that have similar structure.

## Cose ancora da fare:
- scegliere una licenza: GPL / MIT?
- spiegare come produrre i profili, magari mettere uno script
- sistemare il readme

- controllare che gli esempio coincidano con tests/test.sh
