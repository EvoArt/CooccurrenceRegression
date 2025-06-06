# CooccurrenceRegression

[![Build Status](https://github.com/evoart/CooccurrenceRegression.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/evoart/CooccurrenceRegression.jl/actions/workflows/CI.yml?query=branch%3Amaster)


[Bayesian estimation of co-occurrence affinity via dyadic regression](https://www.biorxiv.org/content/10.1101/2024.01.16.575941v1).


Currently, the package re-exports [`Turing.jl`](https://turinglang.org/), meaning that the user can easilly process the [`Chains`](https://turinglang.org/MCMCChains.jl/stable/getting-started/) struct returned by the `cooccurrence_regression` function
### Future plans

The current implementation depends on [`Turing.jl`](https://turinglang.org/) and [`Mooncake.jl`](https://chalk-lab.github.io/Mooncake.jl/stable/), so install times are not as quick as they could be. I aim to streamline the package in the near future.

Depending on the level of interest, I will include additional model structures and options e.g. hierarchical models and splines.