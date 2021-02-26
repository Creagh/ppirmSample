# Poisson Process Infinite Relational Model with Transaction Amounts

## Summary

This repository contains a sample of work for an ongoing project. It includes an implementation of a Bayesian nonparametric model for clustering transactional networks. These networks can possess edges with unique timestamps and transaction amounts. The code framework is implemented in R.

## Examples

[`code/model/poisson_process_irm.R`](code/model/poisson_process_irm.R) contains the functions for forward-sampling from the model, include sampling from a Poisson process and Chinese Restaurant process.

[`code/samplers/sample_all_params.R`](code/samplers/sample_all_params.R) contains the driver code for generating Markov Chain Monte Carlo samples of the model parameters and hyperparameters.

## References

Blundell, C., Beck, J., & Heller, K. A. (2012). [Modelling reciprocating relationships with Hawkes processes.](http://www.gatsby.ucl.ac.uk/~ucgtcbl/papers/BluHelBec2012a.pdf) Advances in Neural Information Processing Systems, 25, 2600-2608.

Briercliffe, C. (2016) [Poisson Process Infinite Relational Model: a Bayesian nonparmetric model for transactional data.](https://open.library.ubc.ca/cIRcle/collections/ubctheses/24/items/1.0308711) Electronic Theses and Dissertations (ETDs), University of British Columbia.

