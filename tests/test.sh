#!/bin/bash

python3 -m acdc-nn Q104H tests/profiles/2ocjA.prof tests/structures/2ocj.pdb A
acdc-nn Q104H  profiles/2ocjA.prof structures/2ocj.pdb  A
echo "Expected: 0.1735"
acdc-nn V51I  profiles/1bsaA.prof structures/1bsa.pdb A  I51V profiles/1bniA.prof structures/1bni.pdb A
echo "Expected: 0.6103"

