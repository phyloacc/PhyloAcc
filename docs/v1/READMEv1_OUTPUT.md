# Format of output files
*prefix*\_rate_postZ\_M\*.txt: posterior median of conserved rate, accelerated rate, probability of gain and loss conservation ($\alpha = P(Z=0\rightarrow Z=1)$ and $\beta = P(Z=1\rightarrow Z=2)$), and posterior probability of being in each latent state on each branch for each element. Columns in the file are:
1. element No. which is the order of the element in the input bed file starting from zero
2. posterior median of accelerated substitution rate
3. posterior median of conserved substitution rate
4. posterior median of $\alpha$
5. posterior median of $\beta$
6. posterior median of $\beta_2 = P(Z=0\rightarrow Z=2)$, which is 0 in current implementation
7. from the 7th column and on, we have four columns for each branch:\*\_0 indicates whether it's "missing"; \*\_1, \*\_2 and \*\_3 are the posterior probability in the background, conserved and accelerated state respectively. The algorithm will prune "missing" branches within outgroup and set the latent states of them to -1 so that the three posterior probabilities are all zero. Column names indicate the branch right above the node and the order of the branch is the same as that in *prefix*\_elem_Z.txt. 

If sampling hyperparameters, the outputs under different hyperparameters will be concatenated to this file. If an element is filtered because of too many alignment gaps (criteria see [README_PARAMETER.md](README_PARAMETER.md)), it will not appear in this file.

*prefix*\_elem_lik.txt: marginal logliklihood for all models (integrating out parameters and latent states). The columns are:
  * *No.*: The order of the element in the input bed file starting from zero
  * *ID*: The element name as in the input bed file
  * *loglik_Null*: marginal logliklihood under the null model
  * *loglik_Acc*: marginal logliklihood under accelerated model
  * *loglik_Full*: marginal logliklihood under the full model
  * *logBF1*: log Bayes factor between null and accelerated model
  * *logBF2*: log Bayes factor between accelerated and full model
  * *loglik_Max_M0, loglik_Max_M1, loglik_Max_M2*: Maximum joint likelihood of Y (observed sequences), r (substitution rates) given Z (latent states) ($\max_{r, Z} P(Y, r|Z)$) under null ($M_0$), accelerated ($M_1$) and full ($M_2$) model respectively.

If updating hyperparameters, the algorithm will only compute the log-likelihood under the full model. When the hyperparameters are updated, the log-likelihoods for each element will be recomputed and concatenated to this file. If an element is filtered because of too many alignment gaps (criteria see [README_PARAMETER.md](README_PARAMETER.md)), all the columns will be zero. If the MCMC algorithm is trapped at some local modes or some other numerical errors occur for some elements, it will return NA.
  

*prefix*\_M\*_elem_Z.txt: maximum loglikhood configurations of latent state Z under null, accelerated and full model, with Z=-1(if the element is 'missing' in the branches of outgroup species),0(background),1(conserved),2(accelerated); each row is an element, ordered same as the input bed file. Output this file if not sample hyperparameters. If an element is filtered because of too many alignment gaps (criteria see [README_PARAMETER.md](README_PARAMETER.md)), all the columns will be zero. 

*prefix*\_hyper.txt: hyperparameters at each iteration, only meaningful if adopting full Bayesian approach by sampling hyperparameters. Columns are:
* *iter*: current iteration
* *nprior_a, nprior_b*: the shape and scale hyperparameter of the gamma prior of accelerated substitution rate
* *cprior_a, cprior_b*: the shape and scale hyperparameter of the gamma prior of conserved substitution rate
* *prior_lrate_a,prior_lrate_b*: hyperparameters of beta prior for probability of loss conservation
* *prior_grate_a,prior_grate_b*: hyperparameters of beta prior for probability of gain conservation

*prefix*\_mcmc_trace_M\_[0-2]_\*.txt: Output this file if verbose=T, which contains the trace of MCMC samples in each iteration for an element. Each row is one iteration and columns are: log-likelihood($P(Y|Z, r)$), accelerated substitution rate, conserved substitution rate, probability of gain and loss conservation and latent state Z of each branch. If sampling hyperparameters, MCMC samples will be concantenated together under different hyperparameters. In the file name, [0-2]: under null, accelerated and full model respectively. *: element No. .

The output files from the extended version including gBGC are slightly different. *prefix*\_rate_postZ\_\*.txt* also contains posterior mode of gBGC coefficient (*gBC* column) and posterior of having gBGC effect on each branch (*\*_B* columns).*prefix*_1_elem_Z.txt has the latent gBGC state on each branch in the first columns.
