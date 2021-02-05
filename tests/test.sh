#!/bin/bash

echo acdc-nn Q104H profiles/2ocjA.prof structures/2ocj.pdb A
acdc-nn Q104H profiles/2ocjA.prof structures/2ocj.pdb A
echo acdc-nn V51I profiles/1bsaA.prof structures/1bsa.pdb A
acdc-nn V51I profiles/1bsaA.prof structures/1bsa.pdb A
echo acdc-nn V51I profiles/1bsaA.prof structures/1bsa.pdb A --inverse I51V profiles/1bniA.prof structures/1bni.pdb A
acdc-nn V51I profiles/1bsaA.prof structures/1bsa.pdb A --inverse I51V profiles/1bniA.prof structures/1bni.pdb A

