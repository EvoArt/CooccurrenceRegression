# CooccurrenceRegression

[![Build Status](https://github.com/evoart/CooccurrenceRegression.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/evoart/CooccurrenceRegression.jl/actions/workflows/CI.yml?query=branch%3Amaster)


[Bayesian estimation of co-occurrence affinity via dyadic regression](https://www.biorxiv.org/content/10.1101/2024.01.16.575941v1).

This package exports one main function `cooccurrence_regression` which performs Bayesian inference via Markov chain Monte Carlo (MCMC) to infer how co-occurrence relationships change in response to one or more explanatory variables, as detailed [here](https://www.biorxiv.org/content/10.1101/2024.01.16.575941v1). The basic inputs for a dataset with N sites and M species are one or many MxM arrays of explanatory variables and an NxM presence/absence matrix. **The presence/absence matrix must not contain any species that are present in either all or none of the sites**. Inference is performed using the [No U-turns (NUTS) sampler](https://jmlr.org/papers/volume15/hoffman14a/hoffman14a.pdf) implemented in [`Turing.jl`](https://turinglang.org/). The function returns a `Chains` struct. Currently, the package re-exports [`Turing.jl`](https://turinglang.org/), meaning that the user can easily process the [`Chains`](https://turinglang.org/MCMCChains.jl/stable/getting-started/) struct returned by the `cooccurrence_regression` function and use all posterior stats and diagnostics provided by [`MCMCChains.jl`](https://turinglang.org/MCMCChains.jl/stable/chains/).

```{julia}
cooccurrence_regression(
    X,                      MxM array of dyadic explanatory variables e.g. a distance matrix 
    Y,                      NxM Presence absence matrix with rows of sites and columns of species
    n_iter=1000;            Number of MCMC samples
    drop_lambda=false,      Whether to drop the random intercept parameters from the returned chain
    alpha_sd =4.0,          Standard deviation for the Gaussian prior for the global intercept
    beta_sd = 1.0,          Standard deviation for the Gaussian prior for the regression coefficients
    lambda_sd =1.0,         Standard deviation for the Gaussian prior for the random intercepts
    n_adapts=100,           Number of adaptation samples for the NUTS inference algorithm
    delta=0.65              Target acceptance rate for the NUTS inference algorithm
)

cooccurrence_regression(
    Xs,                     Vector of MxM arrays of dyadic explanatory variables e.g. distance matrices 
    "
    "
)
```

The returned `Chains` contains posterior samples for the intercept $\alpha$, regression coefficients for each explanatory variable in the order they are supplied $\beta$[i] and each species-level random intercept $\lambda$[i] in the order they appear in the presence absence matrix.

### Future plans

The current implementation depends on [`Turing.jl`](https://turinglang.org/) and [`Mooncake.jl`](https://chalk-lab.github.io/Mooncake.jl/stable/), which are heavy dependencies. I aim to streamline the package in the near future.

Depending on the level of interest, I will include additional model structures and options e.g. hierarchical models and splines.