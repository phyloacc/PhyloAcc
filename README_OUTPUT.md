# Format of output files
*prefix*\_rate_postZ\_*.txt: posterior mode of conserved and accelerated rate and laten state on each branch for each element. Columns in the file are:
1. element ID which is order of the element in the input bed file starting from zero
2. posterior mode of accelerated mutation rate
3. posterior mode of conserved mutation rate
4. from 4th column and on, we have four columns for each branch: *_0 indicates whether it's missing; *_1, *_2 and *_3 are the posterior probability in neutral, conserved and accelerated state respectively. The order of the branch is the same as that in *prefix*_elem_Z.txt.

*prefix*_elem_lik.txt: marginal logliklihood for all models (integrating out parameters). The columns are:
  * *No.*: The order of the element in the input bed file starting from zero
  * *ID*: The element name as in the iput bed file.
  * *loglik_NUll*: marginal logliklihood under null model.
  * *loglik_RES*: marginal logliklihood under accelerated model.
  * *loglik_all*: marginal logliklihood under full model.
  * *log_ratio*: Bayes factor between null and accelerated model. 
  * *loglik_Max1*: Maximum likelihood at $\hat r$ and $\hat Z$ under null model ($M_0$).
  * *loglik_Max2*: Maximum likelihood at $\hat r$ and $\hat Z$ under full model($M_1$).
  * *loglik_Max3*: Maximum likelihood at $\hat r$ and $\hat Z$ under accelerated model($M_2$).
  

*prefix*_1_elem_Z.txt: configurations of latent state Z with maximum loglikhood under full model, Z=-1(missing),0(neutral),1(conserved),2(accelerated); elements are ordered as the input bed file.
