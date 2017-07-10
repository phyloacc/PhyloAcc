# Format of output files
*prefix*\_rate_postZ\_*.txt: posterior mode of conserved and accelerated rate and laten state on each branch for each element. Columns in the file are:
1. element ID which is order of the element in the input bed file starting from zero
2. posterior mode of accelerated mutation rate
3. posterior mode of conserved mutation rate
4. from 4th column and on, we have four columns for each branch: *_0 indicates whether it's missing; *_1, *_2 and *_3 are the posterior probability in neutral, conserved and accelerated state respectively.

*prefix*_elem_lik.txt: marginal logliklihood for full/null models (integrating out parameters), log_ratio is the acceleration score. loglik_Max* is the maximum likelihood.
