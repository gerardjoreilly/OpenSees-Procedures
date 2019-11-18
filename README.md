# Useful OpenSees Procedures

This repository contains a few handy procedures written in TCL for analysis in OpenSees.

If you find any bugs, let me know.

This readme file is not that complete so some things are not documented. I try my best to keep some description of things here though

## Units

This is a list of SI units that can be loaded at the beginning of a model to maintain consistency with units throughout. Similar files exist on the OpenSees Wiki page but this is done to keep everything in kN and m.

## runNRHA3D

This will allow you to run the non-linear analysis in a relatively efficient manner. If the model does not converge initially, it will try many different tolerance tests, algorithms and reduced time steps before finally giving up and marking the case as a non-converging case. This avoids wasting hours of your computer working overtime to converge a model that, when you actually look at the final results, make you realise that that last effort was indeed futile. It also allows you to specify a relatively large time step and conduct analyses much quicker.

In addition, it allows you to specify a collapse drift limit that will stop the analysis should it be exceeded. This works by checking the relative deformations between a specific set of nodes and in the event the collapse drift limit is exceeded, the analysis is stopped and the case is marked as a collapse case. This is quite convenient when doing IDA, for example.

Lastly, the procedure will print a log file of the main results of the analysis. It will let you know what the final state of the model was after the ground motion. 0 means it finished without any problems, 1 means that the collapse drift limit was exceeded and the analysis stopped and -1 means that the analysis could not be converged.

## runNRHA1D

This is the same as above but for SDOF structures. Very useful for doing extensive studies on SDOF oscillators.

## singlePush

This will perform a pushover in a single command. If the analysis doesn't converge, it will try a few options to achieve convergence. If the building loses its lateral capacity from loss of strength of P-Delta collapse, the analysis will be stopped.

## cyclicPush

Same as singlePush but performs full cycles at the target displacement. Can also perform up to 6 cycles at the same target displacement.

## LP_SteelSection

This creates a lumped plasticity element of a standard steel section. List of common sections are included to be called from, including European, UK and US sections.

## LumpedHingeElement

Similar to the LP_SteelSection command but for a more generic situation (i.e. RC sections) where the dimensions of the section and the flexural capacities are known.

## modalAnalysis

This is a very useful file to do a modal analysis of a structural model. It will print the periods on screen and provide a file with the values and also eigenvectors. The outputs of this function are needed for the Matlab-based plotter.

## getSaT

Even though its name says SaT, this function returns a few ground motion intensity measure characteristics. Things like Sa(T1), PGA, PGV can be obtained. Useful to use when you want to interactively get the intensity of a ground motion during analysis (i.e. think IDA)

## IDA_HTF

The TCL-oriented implementation of the Matlab-based script of the same name by the original developers of IDA. This can be used directly in OpenSees for different kinds of intensity measure types. All the user needs to do is provide a working built model, point to a directory of ground motions, define some collapse criteria (i.e. the demand that will signal a collapse in the structure) and the function will automatically carry out an IDA. 

# Licensing
Copyright (C) 2019  Gerard J. O'Reilly

All of these programs are copyrighted under the GNU General Public License as published by the Free Software Foundation, version 2. In short, you can employ them freely (assuming you cite the original source and the relevant publication) but if you want to build upon, extend or re-distribute them, then the derivative software products will also have to be covered under the GPL (i.e. be free software). Other licensing options are negotiable.

# Disclaimer
These programs are distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
