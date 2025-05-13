# Correcting time-variation in of energy gain of AGATA

Fixing the time-varying energy gain of AGATA requires following steps:
- note down which crystals are exhibiting varying gain
- create TH2F time(X) vs energy(Y) matrices
- find correction parameters 
- modify postPSAFilter in ```gen_conf.py``` to include the ```TimeEvoCC.conf``` files   

## Creating time-energy matrices

Executable ```matTimeEvo_AGATA``` is made for this exact purpose, use the ```--help``` for full list of input options. Matrices are constructed from trees produced by TreeBuilder consumer of narval/femul. Matrices are constructed **only for the core hits**.

This executable is looking for trees in the hardcoded path ```run_XXXX/Out/Analysis/Tree_*.root"``` and writes matrix of each crystal into a separate output file to ```timeEvo/temat_XXXX_CRY.root```, unless specified otherwise. Therefore, if you following standard AGATA folder naming and arrangement conventions, **you should be within your ```Replay/``` folder to run this executable**.

Time range of the matrices is plotted in minutes and is deduced from the first and last timestamp in the data. Binning is by default set to 30 seconds/bin, but this can changed using ```--Tbinning``` switch. Energy range is defaulted to 32 000 bins spread from 0 to 8 000keV, but can be changed using the ```--Ebinning``` switch.

**Note 1** Size of resulting matrices can be large enough to crash your computer. If this happens, you normally get the not-very-verbosive "Killed" message in terminal. If that happens, either change binning of the matrices or, preferably, select less crystals for which the matrices are being made for. 

**Note 2** If you want to use the *"poor man's grid search"* of the TimeEvo solver (see next section), it is highly recommended to create large matrices in this step. Size of the input matrix has a big effect on level of improvement one can achieve and grid search is also re-binning the input matrix.

## Solving the time evolution

The ```solveTimeEvo_AGATA.cpp``` is a driver code for the ```libccm``` library made to fix the apparent change of energy in time. See the ```--help``` switch of the ```solveTimeEvo``` executable to see all possible options available. In order to run the code, it is advised to be in the ```Replay``` directory with runs named as ```run_XXXX```, as per AGATA standard. After completion, the output diagnostics and final `TimeEvoCC.conf` configuration files are stored in the ```timeEvo/``` directory. The final parameters can be simply copied to the main directory, as the are already in the correct folder structure `run_XXXX/Conf/CRY/TimeEvoCC.conf`.

The CCM offers a few parameters that can influence level of improvement you can achieve, for details see top level readme. Code can be run in 2 modes:
- With ```--super_settings``` flag; only a single hard-coded set of parameters is used. Very fast but probably no the best result you can achieve.
- Without the flag above, code will run all possible combinations of hard-coded parameters. For this option to work, you need to specify the ```--fit_peak [1] [2] [3]``` to obtain a figure of merit of the applied corrections. Specifically, the FWFM of the specified peak is used. A peak that is NOT in the ROI should be used to avoid over-fitting. This executable takes significantly longer to run, but at least you get a print out as it moves along, which is nice. Parameters reaching best FoM are used to generate the file with final corrections, no true minimization is performed, as I still haven't found a minimizer that can work with discrete values/options. 

Code is running only 1 crystal at the time, because it needs to create a lot of copies of the input matrix that can be quite memory-heavy. To multithread the code for more crystals I recommend running the code with fork in bash, multithreading it in C++ is not possible due to ROOT using thread-unsafe minimizers.  

