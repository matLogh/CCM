#include <iostream>
#include <tuple>
#include <vector>

#include "TGraphSmooth.h"
#include "TMultiGraph.h"
#include "TSpline.h"
#include "TTree.h"

#include "CCMInterpolator.h"
#include "variables.h"

namespace TEC
{

class CCM
{
  public:
    CCM() = delete;
    /// @brief Constructor for CCM class with input parameters for all ROIs (regions of interest), that are used for
    /// corrections.
    /// @param matrix Input TH2D matrix. X-axis expects to be time, Y-axis expects to be energy.
    /// @param _ROIs Vector of regions of interest, that defines regions used for correction (e.i. like peaks when one
    /// does energy calibration)
    CCM(const TH2D &matrix, const std::vector<Region_of_interest> &_ROIs, const double reference_time_low,
        const double reference_time_high);

    /// @brief Constructor for CCM class with input parameters for all ROIs (regions of interest), that are used for
    /// corrections. If you use this constructor, you HAVE TO provide reference vectors using
    /// SetReferenceProjection(...) SetReferenceVector(...) functions, see their documentation for more details.
    CCM(const TH2D &matrix, const std::vector<Region_of_interest> &_ROIs);

    ~CCM();
    /// @brief Calculates energy offsets/shifts
    /// @param fNthreads
    void CalculateEnergyShifts(const unsigned int fNthreads = 8);

    // void AutoRules(int ROI, double sigma_width_acceptance = 3, double dp_width_acceptance = 3);
    // void AutoRules_allROIs(double sigma_width_acceptance = 3, double dp_width_acceptance = 3);

    /// @brief Provide primary function used for correction of the energy, ideally use is polynomial. For correction of
    /// a given "time slice", a first function is used that has number of degrees of freedom less or equal to the number
    /// of valid ROIs
    /// @param fcn Function needs to be defined using TFormula!
    /// @param fit_options Fit options, see
    /// https://root.cern.ch/doc/master/classTGraph.html#aa978c8ee0162e661eae795f6f3a35589
    void SetCorrectionFunction(const TF1 &fcn, const std::string &fit_options);

    /// @brief Save table containing shifts for each ROI for each time slice of time.
    void SaveShiftTable(const std::string &table_filename = "shift_table.dat");

    /// @brief From ROI displacements calculates the correction functions for each time slice. By default only shifts
    /// marked as valid  are used - this can be changed by with calling ConfigureShiftInterpolator() function. In case a
    /// number of valid ROIs is lower then the number of primary function parameters, a first fallback function with
    /// npar=nROIs is used.
    void CalculateCorrectionFits(int time_subdivision = 1);

    /// @brief Fix a provided matrix based on the calculations. CCM object is owner of this matrix and manages it's
    /// deletion.
    /// @return
    TH2D *FixMatrix();

    /// @brief Returns fixed-copy of matrix passed as an input based on the corrections calculated on prior. Correction
    /// functions are calculated from interpolated offsets using TSpline3. The CCM object does not own the returned
    /// matrix, you need to delete it.
    /// @param input_mat
    /// @param valid_only
    /// @return
    TH2D *FixMatrix(const TH2D *input_mat);

    /// @brief Fix the tree using the correction functions - these has to be calculated prior! Provided tree is
    /// cloned into a new TFile and branches are corrected. Normally, energy corrections are done for each time
    /// slice, but time_subdivision will subdivide each time slice into smaller part, correction function for each
    /// part is calculated from the TSpline5 of the shift.
    /// @param tfilename filename containing the tree
    /// @param treename tree name
    /// @param e_branchname energy branch name
    /// @param ts_branchname timestamp branch name - timestamps are expected to be stored as Long64_t!!!
    /// @param time_subdivision
    void FixTree(const std::string &tfilename, const std::string &treename, const std::string &e_branchname,
                 const std::string &ts_branchname, const bool valid_only = true, const int time_subdivision = 1);

    void SaveFitTable(const std::string &data_filename = "fit_table.dat", const std::string &detector_name = "noname");
    void SaveToRootFile(const std::string &outroot_file = "ccm_output.root");

    const ResCont *GetResultContainer(const int ROI_no, const int time_index) const noexcept;

    const int GetNumberOfTimeIndices() const noexcept
    {
        return fXbins;
    };

    const int GetNumberOfROIs() const noexcept
    {
        return V.number_of_ROIs;
    };

    /// @brief To invalidate result for given ROI at given time index. All interpolators are reset.
    /// @param ROI_no
    /// @param time_index
    void SetInvalidResult(const int ROI_no, const int time_index);

    /// @brief Use Gaussian result instead of the polynomial one to get the shift. All interpolators are reset.
    void UseGaussianResult();

    /// @brief Use Gaussian result instead of the polynomial one to get the shift. All interpolators are reset.
    void UsePolynomialResult();

    /// @brief Set reference projection for the shifts. This is useful if you want to use the same shifts for multiple
    /// matrices. Be aware that this will work only if the defined ROI's are " behaving" the same in between the
    /// matrices.
    /// @param projection Projection of a matrix that is otherwise automatically generated using the reference time in
    /// constructor
    void SetReferenceProjection(const TH1 *projection);

    /// @brief Reference vector is calculated automatically from the matrix, but here you have an option to set it by
    /// hand for given ROI
    /// @param ROI_index
    /// @param own_reference_vector
    void SetReferenceVector(const unsigned int ROI_index, const std::vector<double> &own_reference_vector);

    /// @brief Get all the shifts for selected ROI as a function of time. This is useful if you want to modify
    /// calculated shifts using smoothing functions or set own interpolation. If valid_only is set to true, only valid
    /// ROIs are shown
    /// @param roi_index ROI number to be shown
    /// @param valid_only user is owner of the returned TGraph
    TGraph *GetROIShifts(const int roi_index, const bool valid_only = true);

    /// @brief Get shifts for all ROIs at given time. Useful for determination of ideal correction/calibration function.
    /// @param time_bin
    /// @return user is owner of the returned TGraph
    TGraph *GetShiftProfile(const int time_bin, const bool valid_only = true);

    /// @brief produces a TGraph of energy shifts for given ROI but using values from interpolator, not calculated
    /// shifts
    TGraph *GetInterpolationGraph(const int ROI_index, const int subdivide = 10, const bool valid_only = true);

    /// @brief Configure interpolators used to interpolate shifts across the time for given ROI.
    /// @param ROI_index
    /// @param type possible options LINEAR interpolation;
    /// 1) "POLYNOMIAL" interpolation, to be used for small number of points since introduces large
    /// oscillations;
    /// 2) "CSPLINE" cubic spline with natural boundary conditions;
    /// 3) "CSPLINE_PERIODIC" cubic spline with periodic boundary conditions;
    /// 4) "AKIMA", Akima spline with natural boundary conditions ( requires a minimum of 5 points);
    /// 5) "AKIMA_PERIODIC", Akima spline with periodic boundaries ( requires a minimum of 5 points);
    void ConfigureShiftInterpolator(const int ROI_index, const std::string type = "AKIMA",
                                    const bool valid_only = true);
    /// @brief Uses shift of the closest calculated (time) point, shift interpolation disabled
    void DisableInterpolation(const int ROI_index);
    /// @brief Uses interpolation to calculate the shift value in between two calculated (time) points
    void EnableInterpolation(const int ROI_index);

    /// @brief Smoothens the shift values, exact shift values are no longer used. For smoothing the LOWESS (Locally
    /// Weighted Scatterplot Smoothing) algorithm is used.
    /// @param ROI_index which ROI should be smoothed
    /// @param lowess_span Meaning: The span parameter determines the proportion of the data used to fit each local
    /// regression. It is a value between 0 and 1. Effect: A smaller span means that fewer data points are used for each
    /// local regression, leading to a more flexible fit that can capture more detail but may also be more sensitive to
    /// noise. A larger span results in a smoother fit that is less sensitive to noise but may miss finer details.
    /// @param iter Meaning: The iter parameter specifies the number of robustifying iterations to perform. Effect: In
    /// each iteration, the algorithm assigns weights to the data points based on their residuals from the previous
    /// iteration. This helps to reduce the influence of outliers. More iterations can improve robustness but also
    /// increase computation time.
    /// @param delta Meaning: The delta parameter is used to speed up the computations by specifying a distance within
    /// which the weights are considered constant. Effect: When delta is set to a positive value, the algorithm avoids
    /// recalculating weights for points that are close to each other, thus reducing the number of computations. This
    /// can significantly speed up the process for large datasets.
    void SmoothShifts_Lowess(const int ROI_index, const double lowess_span = 0.05, const int iter = 3,
                             const double delta = 0.0);

    /// @brief Smoothens the shift values, exact shift values are no longer used. For smoothing the Kernel Smoother is
    /// used
    /// @param ROI_index which ROI should be smoothed
    /// @param bandwidth Meaning: The bandwidth parameter, often denoted as h, specifies the width of the kernel
    /// function. It determines the range of data points that influence the smoothed value at each point. Effect:
    /// 1) A small bandwidth results in a narrow kernel, which means that only nearby points have a significant
    /// influence on the smoothed value.This can capture more detail and variability in the data but may also lead to a
    /// noisier estimate. 2) A large bandwidth results in a wider kernel, which means that more distant points also
    /// influence the smoothed value. This produces a smoother estimate that is less sensitive to noise but may miss
    /// finer details in the data.
    void SmoothShifts_KernelSmoother(const int ROI_index, const double bandwidth = 2.);

  public:
    const std::string EMPTY_FUNCTION_NAME{"EMPTY_FUNCTION"};

  private:
    int fXbins;
    int fYbins;

    const int fMINIMUM_SMOOTHING_POINTS{3};
    const std::string fDEFAULT_INTERPOLATOR{"AKIMA"};

    TH2D *fFixedTEMAT{nullptr};

    std::vector<std::pair<TF1 *, std::string>> fCorrectionFunctions;
    std::map<int, CCMInterpolator> fInterpolator;
    bool fInterpolatorReset{false};

    std::atomic<int> fThreadTask;
    int fNthreads;
    VarManager V;
    bool fFitDone{false};
    ResCont **ResVec;
    // FitCont *FitVec;
    std::map<double, FitCont> fCorrectionFits;

    void CheckReferenceVector(int ROI_index);
    void CheckReferenceVectors();

    void CopyMatrixContent(TH2D *matrix);
    void CreateReferenceVector(const int ROI_index, const double ROI_time_low, const double ROI_time_high);
    void Normalize(std::vector<double> &v);
    std::pair<TF1 *, std::string> *FindCorrectionFunction(const int nrois);

    /// @brief Get timestamp value from time-bin number.
    /// @param time_slice_index time-slice index, function automatically adds +1 to compensate for root binning (bin
    /// 0 is underflow, bin 1 is first real bin)
    /// @return
    double GetMatrixTime(const int time_slice_index)
    {
        return V.TEMAT->GetXaxis()->GetBinCenter(time_slice_index + 1);
    };

    void BuildInterpolator(const int ROI_index);
    void BuildInterpolators();

    const FitCont CalculateCorrectionFit(const double time);
};

} // namespace TEC
