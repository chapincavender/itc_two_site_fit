# ITC two site fit
Fit ITC data to a model with two interdependent and non-equivalent binding sites.

# Dependencies
This code requires python3 with the numpy and scipy packages.

# Repository structure

## fit_itc_model.py

The script used to fit ITC data, described in detail below.

## 2021_preq1_riboswitch_data

The data and fit parameters for the preQ<sub>1</sub> riboswitches described in this paper. A link to the paper will be provided here upon publication.

# Usage

Usage is `fit_itc_model.py [options] itc_data [itc_data_2 ...]` where `itc_data` is a whitespace-delimited text file containing two columns. The first column contains the volume of each injection in μL, and the second column contains the measured enthalpies in cal mol<sup>-1</sup>.

Options are:
- `-a` Use the approximate treatment of dilution used in the MicroCal PEAQ software. Implies `-i`.
- `-b int` Number of bootstrapping iterations to estimate confidence intervals. Default is `0`.
- `-d float` Width of the prior for the L2 regularization term for log K<sub>D</sub> in units of k<sub>B</sub> T. Default is `1.0`.
- `-e float` Width of the prior for the L2 regularization term for ΔH in units of k<sub>B</sub> T. Default is `1.0`.
- `-g str` Name of file conataining an initial guess for the fit parameters. Default for independent sites is to use `1.0` for the nuisance parameter and `0.0` for all other parameters. Default for interdependent sites is to use the result of the fit with independent sites.
- `-i` Treat binding sites as independent and equivalent.
- `-l float` Initial concentration of ligand in the syringe in μM.
- `-n int` Number of binding sites. Default is `1`.
- `-p float` Penalty for L2 regularization terms in kcal<sup>2</sup> mol<sup>-2</sup>. Default is `0` (i.e. no regularization).
- `-r float` Initial concentration of receptor in the ITC cell in μM.
- `-s int` Number of injections to skip. Skipped injections are included in calculations of volume displaced and heat evolved but are not used to fit parameters. Default is `0`.
- `-t float` Temperature in K. Only used to set prior widths and report ΔG values.
- `-v float` Volume of ITC cell in μL.
- `--print_cost` Print the value of the cost function for the initial guess. Useful for cross validation.
- `--save_bootstrap=str` Name of file to which to writebootstrap samples. Useful for debugging. Default is `None`.

You must give the initial concentrations of ligand and receptor and the volume of the ITC cell using the options `-l`, `-r`, and `-v`.

To perform a global fit of multiple experiments simultaneously, provide additional ITC data files as additional arguments to the script. The options `-l`, `-r`, `-s`, and `-v` can take a comma-separated list to provide different options to each experiment.

# Examples

`fit_itc_model.py -l 100 -r 50 -v 1000 my_itc_data`

Fit ITC data in the file `my_itc_data` using an initial ligand concentration of 100 μM, an initial receptor concentration of 50 μM, and an ITC cell volume of 1000 μL. The parameter estimates and confidence intervals will be reported using the least-squares estimator and the Jacobian of the cost function at the least-squares solution.

`fit_itc_model.py -b 10000 -p 10 -l 100 -r 50 -v 1000 my_itc_data`

Fit ITC data for the same experiment using an L2 regularization penalty of 10 kcal<sup>2</sup> mol<sup>-2</sup>. The parameter estimates and confidence intervals will be reported using the median and (2.5, 97.5) percentiles of the bootstrapped distribution.

`fit_itc_model.py -l 100,100,80 -r 50,40,30 -v 1000 -s 1 my_itc_data my_itc_data_2 my_itc_data_3`

Perform a global fit for three experiments. The ligand concentration is 100 μM in the first experiment, 100 μM in the second experiment, and 80 μM in the third experiment. The receptor concentration is 50 μM in the first experiment, 40 μM in the second experiment, and 30 μM in the third experiment. The ITC cell volume is 1000 μL for all three experiments. The first injection will be skipped in all three experiments.


