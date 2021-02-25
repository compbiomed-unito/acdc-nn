## Description
In this directory, the user can find two notebooks to replicate the paper results on the Ssym dataset and on p53 and myoglobin datasets.
The usage of ACDC-NN presented in the notebook is slightly different from the predictor one. Here, we decided to show in great detail how the method works so that the user could have a deeper comprehension.

### Directories:

  -Ssym_cv.mut: contains the cross-validation sets with the mutations of Ssym considering the sequence similarities with the training sets;
  
  -pdbs: the proteins structure files mutated in Ssym, p53 and myoglobin; moreover it contains the structures created with MODELLER to predict the results on the reverse p53 and myoglobin mutations; 
  
  -profiles: the proteins profiles files mutated in Ssym, p53 and myoglobin;
  
  -weights_cv: the 10 fold cross-validation weights.

### Files:
  -1bz6A_TS_0.mut: set of the myoglobin mutations considering the sequence similarities with the training sets; these mutations can only be tested with the learning set 0;
  
  -2ocjA_TS_5.mut: set of the p53 mutations considering the sequence similarities with the training sets; these mutations can only be tested with the learning set 5;

### Notebooks:
  -ACDC_NN-reproduce_results_ssym.ipynb: notebook with all the code necessary to replicate the results on Ssym. 
  
  -ACDC_NN-reproduce_results_p53_myo.ipynb: notebook with all the code necessary to replicate the results on p53 and myoglobin. 

