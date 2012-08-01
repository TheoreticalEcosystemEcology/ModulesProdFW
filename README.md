# How is ecosystem functioning affected by modules in food webs

Lab members involved: Tim Poisot & Dom Gravel   
Other people involved: Daniel Stouffer

## Synopsis

We assembled a large dataset of 137 binary food webs. For each we counted the frequency of 3-species motifs (community modules). Each food web is simulated using a classical and well parameterized model of food web dynamics, with and without allometric scaling, and the resulting functioning is recorded. We use an ANOVA approach to find out which modules influence functioning, and how.

## Main results

## Contents of the folders

## Technical stuff

The main document is in `LaTeX`, with figures using `pgfplots` and `tikz`, which should be installed by default if you use `TeXLive`.

Data analysis is done with the `xxx.R` script, which requires no additional packages.

The code for the simulations uses the `scipy`, `numpy`, and `networkx` python libraries, and runs under the 2.7 version of python.

The whole project is under version control, please begin each work session by a `git pull` and end each work session with a `git push`.
