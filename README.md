# The Cross-correlation Correction Method (CCM)

Matus Balogh, Slovak Academy of Sciences
matus.balogh@cern.ch

Prerequisities: ROOT (see https://root.cern.ch/building-root), tested on ROOT6.18
The code is rather hard to read, but if you miss some functionality, let me know. Description of the algorithm is given in the NIM A paper, copy is in this repo. 

# Implementation


### Calculating the displacement

As mentioned in the NIM paper, identifying the maximum of dot product for a set of displacement of given ROI and test spectrum yields only integer value, which provides only limited precision. Therefore, in this code a Gaussian fit is used to extract the displacement. The fit is performed on limited number of points around the dot product maximum (*D_M*). The number of accepted points must be whitin the interval *(D_M-D_A/2, D_M+D_A/2)*, where *D_A* is the average dot product calculated from the whole set. The Gaussian function is chosen, as peak-like shape is expected to be within the ROI, and convolution/cross-correlation of two Gaussian functions yiealds another Gaussian. 

Result of the fit provide beside the displacement also the sigma and chi2, which can be used to capture and ignore incorrect fits. Wrong result can be also captured by checking if the resulting displacements falls into pre-selected interval. If any of these conditions is not met, the displacement value of given ROI and test spectrum is flaged and handeled by procedure for extrating the correction parameters, see the *Correction parameters - fitting function* subsection. 


Function listed below are optional, but can be used to set the rules to flag bad fits. They should be used *after* running the CCM.
- *SetRules_sigma()*
- *SetRules_chi2()*
- *AutoRules()*
- *AutoRules_allROIs()*

After setting rules, you need to apply them with 

- ApplyRules()

I advise to use them if you experience bad results and even. Bad results may be however also be caused by the bad settings, usuall suspect is time-binning or presence of another peak/spectral shape similar to the one encompassed by the ROI within the displacement range.

### Correction parameters, the fitting function

Test spectrum correction parameters are extracted by fitting the displacements as a function of energy of their respective ROIs. Currently only hardcoded linear, quadratic and cubic function can be used. Fit function is selected using the:

- *SelectFitFunction(int selectFunction, bool allowLowOrderFit)*

Following options are available for the *selectFunction* based on number of free parameters:

- 0 - autoselect, spline if less than 5 parameters, otherwise cubic
- 1 - linear function with offset = 0
- 2 - linear function
- 3 - quadratic function
- 4 - cubic function

If the extracted displacement was flaged as wrong (see subsection *Calculating the displacement*) the number of fitting points might be less than the number of function's free parameters. This can be resolved by *allowLowOrderFit* allowing use of lower order fit to take place. If use of low order function is not allowed, but enough fit points are available the fit will proceed skipping the flaged points. However, if number of fit points is less than number of free parameters of the function, the flaged fit points are used anyway.

### Getting the results

If you need to use specific function to get the correction parameters you can get the shift table by calling the *BuildShiftTable()* after running the CCM. Its structure will be printed once called. 

Similarily, the fit table can be called by *BuildFitTable()* to get correction parameters directly. 

Very usefull is to fix the matrix used as an input by the *FixMatrix()* function. In some versions of the ROOT might not display the matrix content properly (or at all), but the matrix projections will work.

Events stored in the TTree* can be corrected individually and stored in a new "fixedTree" TTree. Currently, fixedTree will only contain the energy and time branches.  


# HOW TO USE

Example use is in main_example.cxx, available public functions are listed in CCM.h. Replace the "USER_SOURCE" with your own source code in Makefile. Run with ./CCM. 

**All inputs inputs into CCM are in terms of X or Y bin values of the time-energy matrix, not in actual values of energy and time!**

## Functions


### Constructor

**CCM(TH2D* matrix, vector<uint>& reh, vector<uint>& rel, vector<uint>& rp, vector<uint>& rm, uint reference_time_low, uint reference_time_high)** 
- *reh* - ROI energy high, vector-array of high energies (bin value) for each ROI
- *rel* - ROI energy low, vector-array of low energies (bin value) for each ROI
- the *rel - reh* region defines the dimension of the test/reference vectors (vectors in mathematical sense)
- *rp, rm* - ROI plus/minus, vector-array defining the maximum/minimum relative displacement for each ROI 
- reference_time_low, reference_time_high - time window (bin value) for the reference vector 

### Work function

**void StartCCM(unsigned int number_of_threads=4)**
- number_of_threads - specify number of worker threads


### Functions before starting the CCM


**void SetDrawDP(bool val)**
- option to crate PNG files of dot product vs relative displacement graphs
- png images are great for quick visual check of the relative displacement graphs
- it will significantly increase computation time

**void SelectFitFunction(int selectFunction, bool allowLowOrderFit)**
- detailed description in [this subsection](#correction-parameters,-the-fitting-function)

**void ChangeReferenceTime(double reference_time_low, double reference_time_high)**
- change the reference spectrum

### Functions after starting the CCM

**void SetRules_dp(int ROI, double dp)**
- set acceptance threshold for the dot product value extracted by Gaussian fit for given ROI
- if extracted dot product is less than the threshold, fit is flaged as bad
- [see the Calculating the displacement subsection](#calculating-the-displacement)

**void SetRules_chi2(int ROI, double chi2)**
- set acceptance threshold for the chi2 extracted by Gaussian fit for selected ROI
- if extracted chi2 is less than the threshold, fit is flaged as bad
- [see the Calculating the displacement subsection](#calculating-the-displacement)

**void SetRules_sigma(int ROI, double sigma_low, double sigma_high)**
- set acceptance limits for the sigma of the Gaussian dot product fit for given ROI
- if extracted sigma is outside of these limit, fit is flaged as bad
- [see the Calculating the displacement subsection](#calculating-the-displacement)

**void AutoRules(int ROI, double sigma_width_acceptance = 3, double dp_width_acceptance = 3)**
- creates histograms of the dot product and sigma values extracted by Gaussian fit and set appropriate limits based on user input for specified ROI
- [see the Calculating the displacement subsection](#calculating-the-displacement)

**void AutoRules(int ROI, double sigma_width_acceptance = 3, double dp_width_acceptance = 3)**
- creates histograms of the dot product and sigma values extracted by Gaussian fit and set appropriate limits based on user input for all ROIs
- not really a good idea to use, as settings for different ROIs should vary!
- [see the Calculating the displacement subsection](#calculating-the-displacement)

**void ApplyRules()**
- go through specified rules and flag bad fits

### Functions for getting the results

**void BuildShiftTable()**
- builds table of relative displacements
- table header:
```
number_of_ROIs
ROI1_energy ROI2_energy, ...
```
- table body:
```
timeBin_testSpectrum1 ROI1_isOK(bool) ROI1_displacement ROI1_isOK(bool) ROI1_displacement ...
timeBin_testSpectrum2 ROI1_isOK(bool) ROI1_displacement ROI1_isOK(bool) ROI1_displacement ...
...
```

**void BuildFitTable()**
- builds table with correction coefficients
- table header:

```
default_fit_function
```

- table body:

```
timeBin_testSpectrum1  usedFitFcn fit_par_0  fit_par_1  fit_par_2 ...
timeBin_testSpectrum2  usedFitFcn fit_par_0  fit_par_1  fit_par_2 ...
...
```

**void FixMatrix()**
- crates corrected copy of the time-energy matrix in *CCM_files/result.root*

**void FixTree(char* TFile_name, TTree* tree, char* e_branch, char* t_branch)**
- crates corrected tree
- *TFile_name* - name the file for the new tree
- *tree* - pointer to the tree carrying the data
- *e_branch* - name of the energy branch in the *tree*
- *t_branch* - name of the time branch in the *tree*

 



	
