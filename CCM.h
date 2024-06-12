#include <iostream>
#include <tuple>
#include <vector>

#include "TMultiGraph.h"
#include "TSpline.h"
#include "TTree.h"

// #include "CCMInterpolator.h"
#include "variables.h"

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
    /// @brief Calculates offsets/displacements
    /// @param fNthreads
    void CalculateEnergyShifts(const unsigned int fNthreads = 8);

    void AutoRules(int ROI, double sigma_width_acceptance = 3, double dp_width_acceptance = 3);
    void AutoRules_allROIs(double sigma_width_acceptance = 3, double dp_width_acceptance = 3);

    /// @brief Provide primary function used for correction of the energy, ideally use is polynomial
    /// @param fcn
    /// @param fit_options Fit options, see
    /// https://root.cern.ch/doc/master/classTGraph.html#aa978c8ee0162e661eae795f6f3a35589
    void SetCorrectionFunction(const TF1 &fcn, const std::string &fit_options);

    /// @brief Provide secondary functions used for correction of the energy. Multiple functions can be provided by call
    /// this function multiple times. Functions will be called
    /// @param fcn should have LESS parameters than the primary function!
    /// @param fit_options Fit options, see
    /// https://root.cern.ch/doc/master/classTGraph.html#aa978c8ee0162e661eae795f6f3a35589
    void SetFallbackCorrectionFunction(const TF1 &fcn, const std::string &fit_options);

    void SaveShiftTable(const std::string &table_filename = "shift_table.dat");

    /// @brief From ROI displacements calculates the correction functions for each time slice. If valid_only is set to
    /// true, only ROIs marked as valid are used. In case of valid_only=true AND a ROI is marked as invalid AND it is
    /// possible to interpolate invalid ROI offset from neighbors using TSpline5 (only in case they exist and are
    /// valid). In case a number of valid ROIs is lower then the number of primary function parameters, a first fallback
    /// function with npar=nROIs is used.
    /// @param valid_only
    /// @param use_spline
    void PerformFits(const bool valid_only = true, const bool use_spline = false);

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
    TH2D *FixMatrix(const TH2D *input_mat, const bool valid_only = true);

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

    void SetInvalidResult(const int ROI_no, const int time_index);

    /// @brief Use Gaussian result instead of the polynomial one to get the shift
    void UseGaussianResult();

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
    /// @param valid_only
    TGraph GetROIShifts(const int roi_index, const bool valid_only = true);

    /// @brief Get shifts for all ROIs at given time. Useful for determination of ideal correction/calibration function.
    /// @param time_bin
    /// @return
    TGraph GetShiftProfile(const int time_bin, const bool valid_only = true);

  private:
    int fXbins;
    int fYbins;

    TH2D *fFixedTEMAT{nullptr};

    std::vector<std::pair<TF1 *, std::string>> fCorrectionFunctions;
    // std::vector<CCMInterpolator> fInterpolator;

    std::atomic<int> fThreadTask;
    int fNthreads;
    VarManager V;
    bool fFitDone{false};
    ResCont **ResVec;
    FitCont *FitVec;

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

    void BuildInterpolators(const bool valid_only = true);
};