#!/bin/bash -x

acdc-nn single Q104H profiles/2ocjA.prof structures/2ocj.pdb A
acdc-nn single V51I profiles/1bsaA.prof structures/1bsa.pdb A
acdc-nn inverse V51I profiles/1bsaA.prof structures/1bsa.pdb A I51V profiles/1bniA.prof structures/1bni.pdb A
acdc-nn batch batch_test.tsv

