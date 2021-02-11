#!/bin/bash -x

acdc-nn single Q104H tests/profiles/2ocjA.prof.gz tests/structures/2ocj.pdb.gz A
acdc-nn single V51I tests/profiles/1bsaA.prof.gz tests/structures/1bsa.pdb.gz A
acdc-nn inverse V51I tests/profiles/1bsaA.prof.gz tests/structures/1bsa.pdb.gz A I51V tests/profiles/1bniA.prof.gz tests/structures/1bni.pdb.gz A
acdc-nn batch tests/batch_test.tsv

