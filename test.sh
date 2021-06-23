#!/bin/bash -x

acdc-nn seq Q104H tests/profiles/2ocjA.prof.gz
acdc-nn struct Q104H tests/profiles/2ocjA.prof.gz tests/structures/2ocj.pdb.gz A
acdc-nn struct V51I tests/profiles/1bsaA.prof.gz tests/structures/1bsa.pdb.gz A
acdc-nn istruct V51I tests/profiles/1bsaA.prof.gz tests/structures/1bsa.pdb.gz A I51V tests/profiles/1bniA.prof.gz tests/structures/1bni.pdb.gz A
acdc-nn batch tests/batch_test.tsv

