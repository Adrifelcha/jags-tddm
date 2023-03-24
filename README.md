# JAGS-TDDM

Implements the Thurstonian extension of the circular drift diffusion model in JAGS.

`git clone`, and `./makedtddm.sh` should install the module.

In JAGS, `load tddm` as a module.

`X[1:2,] ~ dtddm(drift, bound, ter0, theta0)` where `X[1,]` is choice in radians and `X[2,]` is RT in seconds.

## Vagrantfile

If you use Vagrant, the Vagrantfile spawns a VM that installs and tests this module.
May require a GitHub account with associated keyfile.

## Pre-release

This is an untested pre-release.  Do not use unless you know what you are doing and how to evaluate the accuracy of results.