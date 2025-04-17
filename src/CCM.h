#pragma once

#include <atomic>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall"
#pragma GCC diagnostic ignored "-Wpedantic"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wdouble-promotion"
#include <TGraph.h>

#pragma GCC diagnostic pop

#include "CCMInterpolator.h"
#include "RegionOfInterest.h"
#include "variables.h"

namespace TEC
{
static const int                             MINIMUM_SMOOTHING_POINTS{3};
static const ROOT::Math::Interpolation::Type DEFAULT_INTERPOLATOR{
    ROOT::Math::Interpolation::Type::kAKIMA};

enum SmootherType
{
    KERNEL,
    LOWESS,
    SUPER,
    NONE
};

class CCM
{
  public:
    CCM() = delete;
    /// @brief Constructor for CCM class with input parameters for all ROIs (regions of
    /// interest), that are used for corrections.
    /// @param matrix Input TH2 matrix. X-axis expects to be time, Y-axis expects to be
    /// energy.
    /// @param _ROIs Vector of regions of interest, that defines regions used for
    /// correction (e.i. like peaks when one does energy calibration)
    CCM(std::shared_ptr<TH2>                matrix,
        const std::vector<RegionOfInterest> _ROIs,
        const double                        reference_time_low,
        const double                        reference_time_high);

    /// @brief Constructor for CCM class with input parameters for all ROIs (regions of
    /// interest), that are used for corrections. If you use this constructor, you HAVE TO
    /// provide reference vectors using SetReferenceProjection(...)
    /// SetReferenceVector(...) functions, see their documentation for more details.
    CCM(std::shared_ptr<TH2> matrix, const std::vector<RegionOfInterest> _ROIs);

    CCM(const CCM &other);
    CCM(CCM &&other);

    ~CCM();
    /// @brief Calculates energy offsets/shifts
    /// @param fNthreads
    void CalculateEnergyShifts(const unsigned int fNthreads = 8);

    // void AutoRules(int ROI, double sigma_width_acceptance = 3, double
    // dp_width_acceptance = 3); void AutoRules_allROIs(double sigma_width_acceptance = 3,
    // double dp_width_acceptance = 3);

    /// @brief Provide primary function used for correction of the energy, ideally use is
    /// polynomial. For correction of a given "time slice", a first function is used that
    /// has number of degrees of freedom less or equal to the number of valid ROIs
    /// @param fcn Function needs to be defined using TFormula!
    /// @param fit_options Fit options, see
    /// https://root.cern.ch/doc/master/classTGraph.html#aa978c8ee0162e661eae795f6f3a35589
    void SetCorrectionFunction(const TF1 &fcn, const std::string &fit_options);

    /// @brief Save table containing shifts for each ROI for each time slice of time.
    void SaveShiftTable(const std::string &table_filename = "ccm_shift_table.dat");

    /// @brief From ROI displacements calculates the correction functions for each time
    /// slice. By default only shifts marked as valid  are used - this can be changed by
    /// with calling ConfigureShiftInterpolator() function. In case a number of valid ROIs
    /// is lower then the number of primary function parameters, a first fallback function
    /// with npar=nROIs is used.
    void CalculateCorrectionFits(int time_subdivision = 1);

    /// @brief Fix a provided matrix based on the calculations. CCM object is owner of
    /// this matrix and manages it's deletion.
    /// @return
    std::shared_ptr<TH2> FixMatrix();

    /// @brief Returns fixed-copy of matrix passed as an input based on the corrections
    /// calculated on prior. Correction functions are calculated from interpolated offsets
    /// using TSpline3. The CCM object does not own the returned matrix, you need to
    /// delete it.
    /// @param input_mat
    /// @param valid_only
    /// @return
    std::shared_ptr<TH2> FixMatrix(const TH2 *input_mat);

    /// @brief Fix the tree using the correction functions - these has to be calculated
    /// prior! Provided tree is cloned into a new TFile and branches are corrected.
    /// Normally, energy corrections are done for each time slice, but time_subdivision
    /// will subdivide each time slice into smaller part, correction function for each
    /// part is calculated from the TSpline5 of the shift.
    /// @param tfilename filename containing the tree
    /// @param treename tree name
    /// @param e_branchname energy branch name
    /// @param ts_branchname timestamp branch name - timestamps are expected to be stored
    /// as Long64_t!!!
    /// @param time_subdivision
    void FixTree(const std::string &tfilename,
                 const std::string &treename,
                 const std::string &e_branchname,
                 const std::string &ts_branchname,
                 const bool         valid_only       = true,
                 const int          time_subdivision = 1);

    void SaveFitTable(const std::string &data_filename = "fit_table.dat",
                      const std::string &detector_name = "noname");
    void SaveToRootFile(const std::string &outroot_file = "ccm_output.root");

    const ResCont *GetResultContainer(const size_t ROI_no,
                                      const size_t time_index) const noexcept;

    size_t GetNumberOfTimeIndices() const noexcept { return fXbins; };

    size_t GetNumberOfROIs() const noexcept { return V.ROIs.size(); };

    const std::shared_ptr<TH2> GetInputMatrix() const noexcept { return V.TEMAT; };

    /// @brief Option to manually set valid/invalid status to the result for given ROI and
    /// time (expressed as a bin number of the input matrix)
    /// @param ROI_no
    /// @param time_index
    /// @param valid
    void SetResultStatus(const size_t ROI_no, const size_t time_index, const bool valid);

    /// @brief For each time and ROI cross correlation object calculates an array of dot
    /// products. In theory, maximum DP=best match between reference vector and test
    /// vector and so index of the maximum DP represents the shift one needs to apply for
    /// the correction. However, this is a discrete value limited by the energy binning of
    /// the input matrix. One way to overcome this limitation is to perform Gaussian of
    /// the several points around the maximum DP to get a better shift value. This assumes
    /// that the distribution of DPs around maxDP is gaussian, which it usually is and you
    /// can verify this by checking root file with the results that contain dot product
    /// arrays for each time and ROI. This function ensures that gaussian mu is used as
    /// shift value used for correction.
    ///
    /// The polynomial fit is already performed, so there is little overhead in using
    /// this function.
    ///
    /// Alternative option is calling UsePolynomialResult() function or UseMaxDPResult()
    ///
    /// Note: sigma of the gaussian fit and chi2 are stored in the result container, you
    /// can use them to check the quality of the fit and maybe mark the fits as invalid
    /// using SetResultStatus() function.
    ///
    /// Result of this fit used by default.
    void UseGaussianResult();

    /// @brief For each time and ROI cross correlation object calculates an array of dot
    /// products. In theory, maximum DP=best match between reference vector and test
    /// vector and so index of the maximum DP represents the shift one needs to apply for
    /// the correction. However, this is a discrete value limited by the energy binning of
    /// the input matrix. One way to overcome this limitation is to perform polynomial2 of
    /// the several points around the maximum DP to get a better shift value. Naturally,
    /// this does NOT assume that the distribution of DPs around maxDP is gaussian. This
    /// function ensures that maximum of the polynomial fit is used as shift value used
    /// for correction.
    ///
    /// The polynomial fit is already performed, so there is little overhead in using
    /// this function.
    ///
    /// Alternative option is calling UseGaussianResult() function or
    /// UseMaxDPResult()
    void UsePolynomialResult();

    /// @brief For each time and ROI cross correlation object calculates an array of dot
    /// products. In theory, maximum DP=best match between reference vector and test
    /// vector and so index of the maximum DP represents the shift one needs to apply for
    /// the correction. However, this is a discrete value limited by the energy binning of
    /// the input matrix. To overcome this limitation consider using UseGaussianResult()
    /// or UsePolynomialResult() functions instead.
    void UseMaxDPResult();

    /// @brief Set reference projection for the shifts. This is useful if you want to use
    /// the same shifts for multiple matrices. Be aware that this will work only if the
    /// defined ROI's are " behaving" the same in between the matrices.
    /// @param projection Projection of a matrix that is otherwise automatically generated
    /// using the reference time in constructor
    void SetReferenceProjection(const TH1 *projection);

    /// @brief Reference vector is calculated automatically from the matrix, but here you
    /// have an option to set it by hand for given ROI
    /// @param ROI_index
    /// @param own_reference_vector
    void SetReferenceVector(const unsigned int        ROI_index,
                            const std::vector<float> &own_reference_vector);

    /// @brief Get all the shifts for selected ROI as a function of time. This is useful
    /// if you want to modify calculated shifts using smoothing functions or set own
    /// interpolation. If valid_only is set to true, only valid ROIs are shown
    /// @param roi_index ROI number to be shown
    /// @param valid_only user is owner of the returned TGraph
    std::unique_ptr<TGraph> GetROIShifts(const size_t roi_index,
                                         const bool   valid_only = true);

    /// @brief Get shifts for all ROIs at given time. Useful for determination of ideal
    /// correction/calibration function.
    /// @param time_bin
    /// @return user is owner of the returned TGraph
    std::unique_ptr<TGraph> GetShiftProfile(const int  time_bin,
                                            const bool valid_only = true);

    /// @brief produces a TGraph of energy shifts for given ROI but using values from
    /// interpolator, not calculated shifts
    std::unique_ptr<TGraph> GetInterpolationGraph(const size_t ROI_index,
                                                  const int    subdivide  = 10,
                                                  const bool   valid_only = true);

    /// @brief Configure interpolators used to interpolate shifts across the time for
    /// given ROI.
    /// @param ROI_index
    /// @param type possible options LINEAR interpolation;
    /// 1) "POLYNOMIAL" interpolation, to be used for small number of points since
    /// introduces large oscillations; 2) "CSPLINE" cubic spline with natural boundary
    /// conditions; 3) "CSPLINE_PERIODIC" cubic spline with periodic boundary conditions;
    /// 4) "AKIMA", Akima spline with natural boundary conditions ( requires a minimum of
    /// 5 points); 5) "AKIMA_PERIODIC", Akima spline with periodic boundaries ( requires a
    /// minimum of 5 points);
    /// @param valid_only if set to true, only valid ROIs are used for interpolation
    void ConfigureShiftInterpolator(const size_t      ROI_index,
                                    const std::string type,
                                    const bool        valid_only = true);

    /// @brief Configure interpolators used to interpolate shifts across the time for
    /// given ROI.
    /// @param ROI_index
    /// @param type possible options LINEAR interpolation;
    /// 1) "POLYNOMIAL" interpolation, to be used for small number of points since
    /// introduces large oscillations; 2) "kCSPLINE" cubic spline with natural boundary
    /// conditions; 3) "kCSPLINE_PERIODIC" cubic spline with periodic boundary conditions;
    /// 4) "kAKIMA", Akima spline with natural boundary conditions ( requires a minimum of
    /// 5 points); 5) "kAKIMA_PERIODIC", Akima spline with periodic boundaries ( requires
    /// a minimum of 5 points);
    /// @param valid_only if set to true, only valid ROIs are used for interpolation
    void ConfigureShiftInterpolator(const size_t                          ROI_index,
                                    const ROOT::Math::Interpolation::Type type,
                                    const bool valid_only = true);

    /// @brief Configure interpolators used to interpolate shifts across the time for
    /// all ROIs.
    /// @param type possible options LINEAR interpolation;
    /// 1) "POLYNOMIAL" interpolation, to be used for small number of points since
    /// introduces large oscillations; 2) "CSPLINE" cubic spline with natural boundary
    /// conditions; 3) "CSPLINE_PERIODIC" cubic spline with periodic boundary conditions;
    /// 4) "AKIMA", Akima spline with natural boundary conditions ( requires a minimum of
    /// 5 points); 5) "AKIMA_PERIODIC", Akima spline with periodic boundaries ( requires a
    /// minimum of 5 points);
    /// @param valid_only if set to true, only valid ROIs are used for interpolation
    void ConfigureShiftInterpolator(const std::string type       = "AKIMA",
                                    const bool        valid_only = true);

    /// @brief Configure interpolators used to interpolate shifts across the time for
    /// all ROIs.
    /// @param type possible options LINEAR interpolation;
    /// 1) "POLYNOMIAL" interpolation, to be used for small number of points since
    /// introduces large oscillations; 2) "kCSPLINE" cubic spline with natural boundary
    /// conditions; 3) "kCSPLINE_PERIODIC" cubic spline with periodic boundary conditions;
    /// 4) "kAKIMA", Akima spline with natural boundary conditions ( requires a minimum of
    /// 5 points); 5) "kAKIMA_PERIODIC", Akima spline with periodic boundaries ( requires
    /// a minimum of 5 points);
    /// @param valid_only if set to true, only valid ROIs are used for interpolation
    void ConfigureShiftInterpolator(const ROOT::Math::Interpolation::Type type,
                                    const bool valid_only = true);

    /// @brief Uses shift of the closest calculated (time) point, shift interpolation
    /// disabled
    void DisableInterpolation(const size_t ROI_index);

    /// @brief Uses shift of the closest calculated (time) point, shift interpolation
    /// disabled
    void DisableInterpolation();

    /// @brief Uses interpolation to calculate the shift value in between two calculated
    /// (time) points for a specific ROI. Should be used if the offsets are mostly smooth
    /// and you want to achieve higher precision.
    /// @param ROI_index which ROI should be interpolated
    void EnableInterpolation(const size_t ROI_index);

    /// @brief Uses interpolation to calculate the shift value in between two calculated
    /// (time) points for all previously provided ROIs. Should be used if the offsets are
    /// mostly smooth and you want to achieve higher precision.
    void EnableInterpolation();

    /// @brief Smoothens the calculated shift values that are used to calculate the
    /// correction(/calibration) function. Exact shift values calculated using CCM are no
    /// longer used. Could result in better performance if observed jumps/calibration
    /// changes are not very sharp.
    /// @param smoother Meaning: Type of smoother to be used. Options are
    /// TEC::SmootherType::LOWESS, TEC::SmootherType::KERNEL and TEC::SmootherType::SUPER
    /// . For details please check ROOT documentation.
    /// @param smoother_parameter Meaning: The parameter that defines the smoothing
    /// strength. For LOWESS, it is the span parameter between 0 (minimum smoothing) and 1
    /// (maximum smoothing). For KERNEL, it is the bandwidth parameter that goes from 0
    /// (minimum smoothing) to any value (maximum smoothing, probably not ideal to go
    /// higher than 50, but go ahead and try). For SUPER, it is bass parameter that
    /// controls the smoothness of the fitted curve. Values of up to 10 indicate
    /// increasing smoothness.
    void SmoothShifts(const SmootherType smoother,
                      const double       smoother_parameter,
                      const size_t       ROI_index);

    /// @brief Smoothens the calculated shift values (for all ROIs) that are used to
    /// calculate the correction(/calibration) function. Exact shift values calculated
    /// using CCM are no longer used. Could result in better performance if observed
    /// jumps/calibration changes are not very sharp.
    /// @param smoother Meaning: Type of smoother to be used. Options are
    /// TEC::SmootherType::LOWESS, TEC::SmootherType::KERNEL and TEC::SmootherType::SUPER
    /// . For details please check ROOT documentation.
    void SmoothShifts(const SmootherType smoother, const double smoother_parameter);

    std::map<double, FitCont> GetCorrectionFits() const noexcept
    {
        return fCorrectionFits;
    };

    const FitCont GetCorrectionFit(const double time);

  public:
    static const std::string EMPTY_FUNCTION_NAME;

  private:
    size_t fXbins;
    size_t fYbins;

    std::shared_ptr<TH2> fFixedTEMAT{nullptr};

    std::vector<std::pair<TF1 *, std::string>> fCorrectionFunctions;

    bool fForceRebuildInterpolators{false};

    std::atomic<int> fThreadTask;
    uint             fNthreads;
    VarManager       V;
    bool             fFitDone{false};
    ResCont        **ResVec;

    // protect root fit procedures from concurrent access
    static std::mutex fMtx_ROOTfit;

    std::map<double, FitCont> fCorrectionFits;

  private:
    void CheckReferenceVector(const size_t ROI_index);
    void CheckReferenceVectors();

    void                           CopyMatrixContent();
    void                           CreateReferenceVector(const uint   ROI_index,
                                                         const double ROI_time_low,
                                                         const double ROI_time_high);
    void                           Normalize(std::vector<float> &v);
    std::pair<TF1 *, std::string> *FindCorrectionFunction(const int nrois);

    /// @brief Get timestamp value from time-bin number.
    /// @param time_slice_index time-slice index, function automatically adds +1 to
    /// compensate for root binning (bin 0 is underflow, bin 1 is first real bin)
    /// @return
    double GetMatrixTime(const size_t time_slice_index)
    {
        return V.TEMAT->GetXaxis()->GetBinCenter((int)time_slice_index + 1);
    };

    void BuildInterpolator(const size_t ROI_index);
    void BuildInterpolators();
};

} // namespace TEC
