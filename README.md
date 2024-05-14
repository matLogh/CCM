# The Cross-correlation Correction Method (CCM)

Any mistakes in the codes are purely mine, please report them to <matus.balogh@cern.ch>.

Basic description of the algorithm is given [in the NIM A paper](https://www.sciencedirect.com/science/article/pii/S0168900221003521), but the code here provides a few additional features, such as spline-based corrections. Advantage of the CCM algorithm is that it can be applied for time-varying spectra (such as beta-decays, isomeric decays etc.) and the region used to evaluate energy shifts does not need to be a peak but rather an arbitrary spectral feature unique it shape in it's vicinity. 

# Build
Prerequisities: ROOT (see https://root.cern.ch/building-root), tested on ROOT6.26/15

```bash
git clone https://github.com/matLogh/CCM.git
mkdir build
cd build
cmake ../
make
```

Executables ```simple_example``` and ```optimizer``` will be created that are using the data set from the ```data/``` directory. Source code ```simple_example.cpp``` showcases the most basic usage of the code. The ```optimizer.cpp``` demonstrate an automated brute-force approach to find the ideal parameters of the CCM in order to obtain the best result - in this case defined as a lowest FWHM for given peak. 


# How to use the code
## TL;DR - basics
Algorithm requires following inputs:
- time (X axis) vs energy (Y axis) matrix (TEMAT),
- reference time window,
- at least one Region Of Interest (ROI),
- at least one energy-correcting function.

**TEMAT** is an input matrix with the data. It is required that the matrix has a time on X axis and energy on Y. Both energy- and time-binning are important parameters that significantly affects the resulting correction. It is advised to tune these for the best results.

**Reference time window** denotes a single stretch of time, in which the energy is stable/do not vary in time. It is advised to make reference time window as wide as possible.

Similar to energy calibration, **ROI** defines the energy region which is used to obtain energy shift in order to convert "old" energy to "aligned" energy using the **correction function**. It should encapsulate a peak or another spectral feature that is used to evaluate offsets. Energies are aligned to the ones defined by the reference time window. Multiple ROIs can be defined. Each ROI is defined by
- energy window that completely encapsules a peak or spectral feature you are using for alignment. Window doesn't need to be symmetric, I used usually 50-100 bin wide window.
- energy displacement window that should be at least slightly larger then maximum energy offset of the peak/spectral feature.
- pointer to TEMAT to perform bin conversions


**Correction function** is an equivalent to the function used for energy calibration. Usually, simple polynomial ```0 + [1]*x``` is sufficient for HPGe detectors. It must given as TF1 object that **MUST be build using TFormula constructor**, e.g.:

```cpp
TF1 fcn("gain_fcn", "[0]*x", 0, 4000); // range can be arbitrary, it is set by CCM
TF1 f("my_fcn","[0] + sqrt(x)*[1] + x*[2]", 0 ,1); 
```
Multiple correction function can be specified, correction is performed with the first function that has number of degrees of freedom equal or less then number valid ROIs. 

Once these objects are set, it is necessary to call ```CalculateEnergyShifts(int nthreads)``` function that determine the energy offsets of ROIs for each time-bin of the TEMAT. Correction, or pseudo-energy-calibration functions for each time-slice is calculated by calling ```PerformFits( bool valid_only, bool use_spline)``` After that is done, it is up to you to decide what output do you prefer, see next sections.
***
## Results
To evaluate effectiveness of the CCM, the fastest way is [matrix correction](#matrix-corrections). It is useful especially for evaluation of the best parameters for CCM. 

CCM can also apply corrections for an arbitrary new [TTree](#tree-corrections).

Further, [text-based tables](#table-corrections) can be produced containing either just the shifts of each ROI, or a complete "recipe" containing list of correction functions and their parameters for each time period. 

Lastly, more complete set of results can be also saved into [CCM ROOT file](#root-file).
### Matrix corrections

The most basic output of the  ```CCM fix(...)``` object can be produced by correcting the input matrix
```cpp
TH2D *fixed_matrix = fix.FixMatrix();
```

Another way to fix the matrix is to utilize [spline interpolation](#spline-interpolation). To use it, one needs to create an identical matrix (data-wise) to original TEMAT, but with finer binning and pass it to the 
```cpp
TH2D *FixMatrix(const TH2D *input_mat, const bool valid_only = true)
``` 
function. [Spline interpolation](#spline-interpolation) is used to evaluate energy offsets for additional time bins. 

### Tree corrections
Event by event corrections can be applied by user-supplied TTree. CCM creates a clone of the tree (including unused branches) and apply correction on the content of the energy branch. TTree is expected to contain a single-value branches for time(stamp) and energy. For enhanced corrections a [spline interpolation](#spline-interpolation) can be used by specifying ```time_subdivision``` parameter, which essentially divides the range of each time-bin of the TEMAT into given number of sub-ranges with unique correction function.  

```cpp
void FixTree(const std::string &tfilename, const std::string &treename, const std::string &e_branchname,
                         const std::string &ts_branchname, const bool valid_only = true,
                         const int time_subdivision = 1);
```
### Table corrections
#### Shift table

Shift table is created by calling 
```cpp
void SaveShiftTable(const std::string &table_filename);
```
in following format (written also in the table)
```
# number of ROIs (regions of interest) in the matrix
# desired_energy_ROI0 desired_energy_ROI1 ...
#########
# time_bin time_start time_end ROI0_ok ROI0_shift ROI1_ok ROI1_shift ...
```

#### Fit table
Fit table is created by calling
```cpp
void SaveFitTable(const std::string &data_filename, const std::string &detector_name)
```
in following format (written also in the table)
```
# detector_name
# time range covered by this file
# number of correction functions used
# fcn_name number_of_parameters function_equation(TFormula style)
##############################
# time_start time_stop fcn_name number_of_parameters par0 par1 ...
```

#### ROOT file

A ROOT file containing the most detailed information about of the calculations can be created by calling
```cpp
void SaveToRootFile(const std::string &outroot_file = "ccm_output.root");
```
It contains 
- original TEMAT,
- corrected TEMAT (if created), 
- TTree for each ROI,
- list of TF1 correction functions used for energy corrections,
- TTree with parameters for energy corrections for each time bin of original TEMAT.

The ROI trees grants access to full information on the cross-correlation method used to evaluate energy shift for all time bins of original TEMAT. This includes:
- *fit_valid*; flag for validity of the calculated shift,
- *bin_shift*; shift in "bin-energy values",
- *energy_shift*; shift in energy units,
- *dot_product*; value of maximum dot product, 
- *dp_vector*; n-dimensional vector of dot products calculated between reference vector and test vector (n-dimension = number of shifts of test vector),
- *gfit_chi2*; chi2 of the [Gaussian fit](#validity-of-the-calculated-shift) of the *dp_vector* around the region of *dot_product*,
- *gfit_sigma*; sigma of the [Gaussian fit](#validity-of-the-calculated-shift) of the *dp_vector* around the region of *dot_product*,  
- *gfit_mu*; mu of the [Gaussian fit](#validity-of-the-calculated-shift) of the *dp_vector* around the region of *dot_product*,
- *time*; time


***
# More detailed description
## Calculating displacement

As mentioned in the NIM paper, the maximum dot product for a set of displacement of a given ROI-vector and reference-vector yields only integer value, which provides only limited precision. A floating point value can be obtained by performing a fit of the dot product distribution. However, precise fit would require apriori knowledge of the analytical function resulting from convolution of the "spectral features" (defined in the ROI bounds) with itself, but this is not possible to achieve in general. Therefore, a simple *pol2* fit of 9 points around the maximum dot product and floating point value of the shift is calculated from the maximum of the fitted function. 

Although it is not possible to know analytical function describing the convolution, it is reasonable to expect that resulting dot product distribution is peak-shaped. Therefore, a Gaussian fit of limited number of values around the maximum is performed in addition to the polynomial but it is not used as a default value of the shift. This can be changed by calling 
```cpp
void UseGaussianResult();
```


## Validity of the calculated shift
By default all shifts are marked as valid, unless the energy projection of given time-bin is 0; such bins are expected as one wants to make sure that all the data are inside the matrix. 

User can however mark some shifts as invalid manually by calling 
```cpp
void SetInvalidResult(const int ROI_no, const int time_index);
```
where the input is the ROI number (indexing from 0) and the index of time-projected spectrum (index from 0, basically bin number minus 1). The decision can be made based on the data stored in the ```ResCont``` structure obtained by calling 
```cpp
const ResCont *GetResultContainer(const int ROI_no, const int time_index);
```
Function returns ```nullptr``` in case of wrong inputs. This container holds 
- *fit_valid*; flag for validity of the calculated shift,
- *bin_shift*; shift in "bin-energy values",
- *energy_shift*; shift in energy units,
- *dot_product*; value of maximum dot product, 
- *dp_vector*; n-dimensional vector of dot products calculated between reference vector and test vector (n-dimension = number of shifts of test vector)

Moreover, it contain also details on Gaussian fit of the dp_vector points around the dot product maximum. As mentioned in the previous section, it is not possible to know analytical function describing the convolution, however it is reasonable to expect that resulting dot product distribution is peak-shaped. Therefore, a Gaussian fit of limited number of values around the maximum is reasonable. Idea behind this is that values like chi2 and sigma should be within a reasonable range - you can see this by quickly plotting them (individually for each ROI) from the [ROOT file](#root-file). The list of variables available in the *ResCont* are 

- *gfit_chi2*; chi2 of the [Gaussian fit](#validity-of-the-calculated-shift) of the *dp_vector* around the region of *dot_product*,
- *gfit_sigma*; sigma of the [Gaussian fit](#validity-of-the-calculated-shift) of the *dp_vector* around the region of *dot_product*,  
- *gfit_mu*; mu of the [Gaussian fit](#validity-of-the-calculated-shift) of the *dp_vector* around the region of *dot_product*,
- *time*; time

## Spline interpolation
Simply put, what this any any other code for time correction does, is that it subdivides the data set into smaller ones and "recalibrate" them individually. In case of CCM, the subdivision is given by the number of time-bins of the TEMAT. Number of time-bins/subdivisions is however limited by available statistics. 

The observed time instability is (at least in my experience) usually continuous (see exception in the example data). If that is the case, it can be used to our advantage: for a given ROI, we can interpolate between the calculated shifts to artificially obtain much higher subdivision in time. If enabled, the new "artificial" shifts are calculated using a *TSpline3*, which is constructed from the available data. The spline can be also used to replace [invalid shifts](validity-of-the-calculated-shift), but only if there exist valid points around it (in lower and higher time slices). Spline is never used to extrapolate the values outside of the first and last valid time range (given by time of the respective bin centers).
	
## Own reference vector

In case you need to correct multiple detectors with same properties (e.g. gamma detectors in an array) you manually set your reference vector used in the calculations by calling
```cpp
void SetReferenceVector(const std::vector<double> &own_reference_vector)
```


## Fitter
In the source files I added "my" pure C++ fitter using the Theuerkauf peak model used in the [HDTV program](https://github.com/janmayer/hdtv/tree/master). Large part of the fitter code is borrowed from there.