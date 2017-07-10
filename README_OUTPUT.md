# Format of output files
*prefix*\_rate_postZ\_*.txt: posterior mode of conserved and accelerated rate and laten state on each branch for each element. Columns in the file are:
1. element ID which is order of the element in the input bed file starting from zero
2. posterior mode of accelerated mutation rate
3. posterior mode of conserved mutation rate
4. from 4th column and on, we have four columns for each branch: *_0 indicates whether it's missing; *_1, *_2 and *_3 are the posterior probability in neutral, conserved and accelerated state respectively. Column order of the branch is the same as that in *prefix*_elem_Z.txt. *

*prefix*_elem_lik.txt: marginal logliklihood for all models (integrating out parameters). The columns are:
  * *No.*: The order of the element in the input bed file starting from zero
  * *ID*: The element name as in the iput bed file
  * *loglik_NUll*: marginal logliklihood under null model
  * *loglik_RES*: marginal logliklihood under accelerated model
  * *loglik_all*: marginal logliklihood under full model
  * *log_ratio*: Bayes factor between null and accelerated model
  * *loglik_Max1, loglik_Max2, loglik_Max3*: Maximum likelihood at $\hat r$ and $\hat Z$ under null ($M_0$), full ($M_2$) and accelerated ($M_1$) model respectively.
  if updating hyperparameters, the algorithm will only compute the log-likelihood under full model. When the hyperparameters are updated, the log-likelihood for each element will be recomputed and concatenated to this file.
  

*prefix*_1_elem_Z.txt: maximum loglikhood configurations of latent state Z under full model, with Z=-1(if the element is 'missing' in the lineages of outgroup species),0(neutral),1(conserved),2(accelerated); each row is an element, ordered same as the input bed file.

*prefix*_hyper.txt: hyperparameters at each iteration, only meaningful if adopting full Bayesian approach by sampling hyperparameters. Columns are:
* *iter*: the number of iteration
* *nprior_a, nprior_b*: the shape and scale hyperparameter of the gamma prior of accelerated mutation rate
* *cprior_a, cprior_b*: the shape and scale hyperparameter of the gamma prior of conserved mutation rate
* *indel, indel2*: empty columns, place holders
* *grate*: rate of gain conservation, i.e. the prior of transition probability of $Z=0$ to $Z=1$
* *lrate*: rate of loss conservation, i.e. the prior of transition probability of $Z=1$ to $Z=2$

