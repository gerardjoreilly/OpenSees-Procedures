# Useful OpenSees Procedures

This repository contains a few handy procedures written in TCL for analysis in OpenSees.

If you find any bugs, let me know.

## Units

This is a list of SI units that can be loaded at the beginning of a model to maintain consistency with units throughout. Similar files exist on the OpenSees Wiki page but this is done to keep everything in kN and m.

## runNRHA3D

This will allow you to run the non-linear analysis in a relatively efficient manner. If the model does not converge initially, it will try many different tolerance tests, algorithms and reduced time steps before finally giving up and marking the case as a non-converging case. This avoids wasting hours of your computer working overtime to converge a model that, when you actually look at the final results, make you realise that that last effort was indeed futile. It also allows you to specify a relatively large time step and conduct analyses much quicker.

In addition, it allows you to specify a collapse drift limit that will stop the analysis should it be exceeded. This works by checking the relative deformations between a specific set of nodes and in the event the collapse drift limit is exceeded, the analysis is stopped and the case is marked as a collapse case. This is quite convenient when doing IDA, for example.

Lastly, the procedure will print a log file of the main results of the analysis. It will let you know what the final state of the model was after the ground motion. 0 means it finished without any problems, 1 means that the collapse drift limit was exceeded and the analysis stopped and -1 means that the analysis could not be converged.

## singlePush

This will perform a pushover in a single command. If the analysis doesn't converge, it will try a few options to achieve convergence. If the building loses its lateral capacity from loss of strength of P-Delta collapse, the analysis will be stopped.

## cyclicPush

Same as singlePush but performs full cycles at the target displacement. Can also perform up to 6 cycles at the same target displacement.



# Licensing
Copyright (C) 2017  Gerard J. O'Reilly

All of these programs are copyrighted under the GNU General Public License as published by the Free Software Foundation, version 2. In short, you can employ them freely (assuming you cite the original source and the relevant publication) but if you want to build upon, extend or re-distribute them, then the derivative software products will also have to be covered under the GPL (i.e. be free software). Other licensing options are negotiable.

# Disclaimer
These programs are distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
