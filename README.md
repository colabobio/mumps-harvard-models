## Code for SEIR  mumps models

This folder contains R notebooks and scripts to calculate the Maximum Likelihood Estimation (MLE) of the parameters in a set of epidemiological models for the mumps outbreaks at Harvard University in 2016. The package [POMP](https://kingaa.github.io/pomp/) was used for the estimation.

### How to use

The .Rmd notebooks can be run from within RStudio, using the desired properties files. Otherwise, the .R scripts for intervention and baseline models can be run from the command line as follows:

```
./run_inter.sh fast.properties
```

```
./run_base fast.properties
```

### Versions used

* R: 4.0.5
* pomp: 3.3