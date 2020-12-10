# Drifty

Drifty is a tool created for the purposes of quickly and easily performing simulation-driven kinematic corrections to correct for any systematic drift introduced by the reconstruction of data.

It is a combination of two python classes which have been designed to work as an interactive "package" when used within Jupyter Notebook or Jupyter Lab, which allow python scripts to be ran in your browser and serves as a shortcut to a GUI.  

## Approach
Kinematic corrections are carried out by determining the relationship between a kinematic variable, and the 'drift' which is observed in simulated data when comparing what was generated with what is reconstructed.

The approach is as follows:
   1. Plot a 2D histograms of the difference(reconstructed - generated) vs. reconstructed for a kinematic variable which we want to correct (eg. the x-component of the momentum of a given particle in a given detector).

   2. Take N-slices in the x-axis (reconstructed), and plot the projections on the y-axis (difference).

   3. Fit each of these slices with a Gaussian distribution to ascertain the mean value of the difference in that projection.

   4. Plot the fitted means against the center position of the slices, and fit the resulting distribution.

   5. The resulting fit gives us a relationship with which we can correct for any systematic effect introduced by the reconstruction.


## Usage
The CalcCorrections class has has been designed to be largely automated.  When an instance of the class is created, the full process is executed, and results are displayed.  An example of the output is included below. 

![Example Output](ExampleResult.png)

### Possible adjustments
After the initial result is made, the user then has several options to fine-tune the result interactively by using the following methods of the class.

Action | Syntax
------------ | -------------
Limit fit range of slice | `set_slice_fitrange(slc_idx, x_low, x_high)`
Rebin a slice | `rebin_slice(slc_idx, factor)`
Hard-code initial fit parameters of a slice (if fit failed) | `set_slice_fitparams(slc_idx, params[dict])`
Restore a slice to its original state (pre-adjustments) | `restore_slice(slc_idx)`
Completely exclude a slice from final result | `veto_slice(slc_idx)`
Undo exclusion of a slice from final result | `unveto_slice(slc_idx)`
Choose a fit model for final result | `change_fit_model("<model>")`
Choose a user-defined function to model the final result | `user_def_model(function)`
 
 
Where:
   * `slc_idx` is the index of the target slice, displayed in red at the top-left of the projections.
   * `change_fit_model()` accepts `"linear"`, `"pol2"`, `"pol3"`, `"pol4"`.

## Getting Started
The tool has been designed such that even with minimal Python experience, a user will be able to quickly set up and perform these kinematic corrections.   Most of the code exists "under the hood" in two files ([SimKinCorr.py](./SimKinCorr.py), [utils.py](./utils.py)) which need to be imported into your notebook

[CalcCorr_Tutorial.ipynb](./CalcCorr_Tutorial.ipynb) contains a working example of how to use the class.  This can be used to explore the functionality, or used as a template for performing an analysis.  For the JLab user, the data required for this tutorial can be downloaded from: `/work/clas12/pnaidoo/CalcCorr/tutorial_data.root`

[GettingSetUp.md](./GettingSetUp.md) exists for the purposes of those unfamiliar with Python or Jupyter. It contains all the steps required to set up a virtual environment with the packages required by the class and offers an explanation on how to run CalcCorr_Tutorial.ipynb in Jupyter-Notebook.
