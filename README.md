# PhytoplanktonStoichiometry
Code for manuscript: Aarons and Follows (2026), Common macromolecular allocation strategies drive phytoplankton elemental stoichiometry across taxonomy and nutrient regimes. Manuscript submitted for publication.

## Running Code:
Download the associated phytoplankton stoichiometry literature compilation (https://doi.org/10.5281/zenodo.17289136) and save in same folder as code.

ParameterFitting.py estimates model parameter values using a Bayesian Markov chain Monte Carlo (MCMC) approach, based on observations from the literature compilation. Parameter estimation is performed using several alternative data subsets, with the resulting values available in the MCMC_Output folder.

Main.py runs the study analysis and generates all plots and tables. The cellular allocation model is applied to the data compilation, using the model parameters fit above. The full equations for the modified CFM-Phyto model are available in the CellModel.py file.
