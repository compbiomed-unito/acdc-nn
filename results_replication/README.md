# Content of this folder
Directories:

  -Ssym_cv.mut: contains the cross-validation sets with the mutations of Ssym divided by homology with the training sets;
  
  -pdbs: the structure files of the proteins mutated in Ssym;
  
  -profiles: the profiles files of the proteins mutated in Ssym;
  
  -weights_cv: the 10 different weights to use for each fold of the cross-validation.

Notebook:
ACDC_NN-reproduce_results.ipynb: notebook with all the code necessary to replicate the results. 

The reader that would like to replicate the results we obtained, could simply run the notebook. The usage of ACDC-NN presented in the notebook is slightly different from what one would do if it used it as a predictor. We decided to show in great detail how the method works so that the user could have a greater comprehension.
