# Fermionization

This folder contains the Matlab scripts and resulting data used to create the figures in the publication

*Fermionization of a few-body Bose system immersed into a Bose-Einstein condensate*  
Tim Keller, Thomás Fogarty, Thomas Busch  
SciPost Phys. **15**, 095 (2023).  
[doi: 10.21468/SciPostPhys.15.3.095](https://doi.org/10.21468/SciPostPhys.15.3.095)

Simulations were performed with *MATLAB R2019a* and can be executed e.g. via

	cd(fullfile(pwd,'matlab/'))
	fermionization_main()

Depending on the available hardware the outer loops should be parallelized as much as possible, e.g. by running separate instances for each combination of N = [2, 3] and opt = ["superfluid", "pinned"] as well as further replacing the loop over gMIX by

    gMIX = gMIXarr(${SLURM_ARRAY_TASK_ID})  
	
in each instance when running the code as an array job on a cluster using slurm. For the given position grid, simulations require at least 

	N = 2: 8GB of RAM
	N = 3: 24GB of RAM

and 4 to 8 CPU cores to benefit from Matlab's built-in multi-threading for its base functions.  


External dependencies are the additional files *'fftdef.m'* and *'v2struct.m'* which are included in the *matlab/* folder. 
The colormap used in Fig. 2 is [*'Smooth Cool Warm'*](https://www.kennethmoreland.com/color-advice/) by Kenneth Moreland.

The figures are created via tikz in the .tex file for the publication itself. 
Tikz settings need to be included in the preamble via

	\input{figures/tikz_settings.tex}

The figures are then created at the desired locations via

	\begin{figure}
	\centering
	\input{figures/fermionization_Fig1.tikz}
	\caption{Caption}
	\label{fig:Fig1}
	\end{figure}
	
For questions please contact tim.keller@aoni.waseda.jp	