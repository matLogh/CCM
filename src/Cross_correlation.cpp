/*
list of throw exceptions:

1 = vectors fed to the constructor are not the same size
2 = normalization of sample vector error = its field of zeroes
3 = vectors in dot product are not the same size
*/

// #include "CheckCCM.h"
#include "Cross_correlation.h"
#include "variables.h"

void TEC::CrossCorrel::Process(unsigned int      thread_id,
                               std::atomic<int> *thread_task,
                               std::mutex       &mtx_task,
                               std::mutex       &mtx_fit)
{
    static const int JOB_BATCH_SIZE = 100;

    fThread_id = thread_id;

    int start_task;
    int end_task;
    while (true)
    {

        mtx_task.lock();
        // if no tasks left
        if (*thread_task >= (int)V->total_tasks)
        {
            mtx_task.unlock();
            break;
        }

        start_task   = *thread_task;
        end_task     = *thread_task + JOB_BATCH_SIZE <= (int)V->total_tasks
                           ? *thread_task + JOB_BATCH_SIZE
                           : (int)V->total_tasks;
        *thread_task = end_task;

        // current_task = ++(*thread_task);
        // if (current_task >= V->total_tasks)
        // {
        //     mtx_task.unlock();
        //     break;
        // }

        // std::cout << "Thread " << thread_id << " is processing tasks from " <<
        // start_task
        //           << " to " << end_task << std::endl;

        current_task = start_task;

        // std::cout << std::setprecision(3)
        //           << "Progress: " << (float)current_task / (float)V->total_tasks * 100.
        //           << "%                \r" << std::flush;
        mtx_task.unlock();

        for (current_task = start_task; current_task < end_task; current_task++)
        {

            ROIAnalysis(current_task % (int)V->ROIs.size(),
                        current_task / (int)V->ROIs.size(), mtx_fit);
        }
    }
}

/****************************************************
value -999 passed to _dp_vec indicates error and this
must be disregarded later
****************************************************/
void TEC::CrossCorrel::ROIAnalysis(const int   ROI_index,
                                   const int   time,
                                   std::mutex &mtx_fit)
{
    // create data vector with size that includes whole area around the floating vector
    std::vector<float> data_vec = GetDataVec(ROI_index, time);

    // vector of dot products
    _dp_vec.clear();
    _dp_vec.resize((size_t)V->ROIs[(uint)ROI_index].displacement_steps);
    // temporary vector holding current data
    std::vector<float> temp_data(V->ROIs[(uint)ROI_index].vector_dimension);
    // cycle through through the whole region, calculate cross correlation
    for (std::size_t shift = 0; shift < V->ROIs[(uint)ROI_index].displacement_steps; shift++)
    {
        std::memcpy(temp_data.data(), &data_vec[shift],
                    V->ROIs[(uint)ROI_index].vector_dimension * sizeof(float));
        // temp_data = std::vector<float>(
        //     &data_vec[shift],
        //     &data_vec[shift + V->ROIs[(uint)ROI_index].vector_dimension]);

        if (this->Normalize(temp_data) == 0) // 0 == its OK
            _dp_vec[shift] =
                this->DotProduct(temp_data, V->sample_vector[(uint)ROI_index]);
        else { _dp_vec[shift] = -999; }
    }

    this->SaveToContainer(time, ROI_index, mtx_fit);
}

std::vector<float> TEC::CrossCorrel::GetDataVec(const int ROI_index, const int time)
{
    std::vector<float> data_vec(
        &V->TEMATarr[ROI_index][time][0],
        &V->TEMATarr[ROI_index][time][V->ROIs[static_cast<size_t>(ROI_index)].displacement_range]);
    return data_vec;
}

/****************************************************
retval:
   0 = OK
  -1 = input vector is full zeroes
****************************************************/
int TEC::CrossCorrel::Normalize(std::vector<float> &v)
{
    _norm = 0;
    for (uint i = 0; i < v.size(); i++) { _norm += (v[i] * v[i]); }
    if (_norm == 0) { return -1; }

    _norm = 1.f / sqrtf(_norm);
    for (uint i = 0; i < v.size(); i++) { v[i] = v[i] * _norm; }
    return 0;
}

float TEC::CrossCorrel::DotProduct(const std::vector<float> &v1,
                                    const std::vector<float> &v2)
{
    _dp = 0;
    if (v1.size() != v2.size())
    {
        std::cerr << "WRONG VECTOR SIZE" << std::endl;
        throw 3;
        exit(3);
    }
    for (uint i = 0; i < v1.size(); i++) { _dp += v1[i] * v2[i]; }
    return _dp;
}

void TEC::CrossCorrel::SaveToContainer(const int   time,
                                       const int   ROI_index,
                                       std::mutex &mtx_fit)
{
    // check if dp has valid values
    if (static_cast<std::size_t>(std::count(_dp_vec.begin(), _dp_vec.end(), -999)) == _dp_vec.size())
    {
        ResVec[ROI_index][time].isValid = false;
        // std::cout << "\nTime " << time << " contains " << std::count(_dp_vec.begin(),
        // _dp_vec.end(), -999)
        //           << " invalid values" << std::endl;
        return;
    }
    std::replace_if(
        _dp_vec.begin(), _dp_vec.end(), [](double i) { return i == -999; }, 0.0);

    // copy dp to container
    ResVec[ROI_index][time].dp_vec.resize(_dp_vec.size());
    memcpy(&ResVec[ROI_index][time].dp_vec[0], &_dp_vec[0], _dp_vec.size() * sizeof(float));

    // get shift using pol2 fit of 9 points
    this->GetShift_Poly2(ResVec[ROI_index][time].poly_shift, ResVec[ROI_index][time].dp,
                         mtx_fit);
    // fit with gaussian to estimate property of the peak
    mtx_fit.lock();
    this->GetShift_Gaussian(ResVec[ROI_index][time].gfit_chi2,
                            ResVec[ROI_index][time].gfit_sigma,
                            ResVec[ROI_index][time].gfit_mu);
    mtx_fit.unlock();
    // set the default shift value - use gaussian as default
    ResVec[ROI_index][time].bin_shift =
        ResVec[ROI_index][time].gfit_mu + V->ROIs[static_cast<size_t>(ROI_index)].base_shift_value;
    ResVec[ROI_index][time].energy_shift =
        V->TEMAT->GetYaxis()->GetBinWidth(1) * ResVec[ROI_index][time].bin_shift;
}

void TEC::CrossCorrel::GetShift_Poly2(double &shift, double &dp, std::mutex &mtx_fit)
{
    const static int ndim = 9;
    static_assert(ndim % 2 == 1, "ndim must be odd");

    int index =
        static_cast<int>(std::distance(_dp_vec.begin(), std::max_element(_dp_vec.begin(), _dp_vec.end())));
    if (_dp_vec.size() < ndim)
    {
        shift = static_cast<double>(index);
        dp    = _dp_vec[static_cast<std::size_t>(index)];
        return;
    }

    double arr_x[ndim];
    double arr_y[ndim];
    int    correction = 0;
    if (index + ndim / 2 > static_cast<int>(_dp_vec.size()) - 1)
        correction = (static_cast<int>(_dp_vec.size()) - 1) - (index + ndim / 2);
    if (index - ndim / 2 < 0) correction = ndim / 2 - index;
    // for (int i = -ndim / 2 + correction; i <= ndim / 2 + correction; i++)
    for (int i = 0; i < ndim; i++)
    {
        arr_x[i] = index + i + correction - ndim / 2;
        arr_y[i] = _dp_vec[static_cast<std::size_t>(arr_x[i])];
    }
    TGraph gr(ndim, arr_x, arr_y);

    mtx_fit.lock();
    {
        TF1 fcn(Form("shift_fitfcn_%i", fThread_id), "pol2", arr_x[0], arr_x[ndim - 1]);
        fcn.SetNpx(ndim * 1000);
        // for (int i = 0; i < ndim; i++)
        // std::cout << arr_x[i] << " " << arr_y[i] << " ";
        // std::cout << std::endl;
        gr.Fit(&fcn, "RQNC");
        // std::cout << "linear fit end" << std::endl;
        // if (fcn.GetMaximum(arr_x[0], arr_x[ndim - 1]) < _dp_vec[index])
        // {
        //     shift = static_cast<double>(index);
        //     dp = _dp_vec[index];
        //     mtx_fit.unlock();
        //     return;
        // }

        shift = fcn.GetMaximumX();
        dp    = fcn.GetMaximum();
    }
    mtx_fit.unlock();
}

void TEC::CrossCorrel::GetShift_Gaussian(double &rchi2, double &sigma, double &mu)
{
    int    nbins = static_cast<int>(_dp_vec.size());
    TH1D   h(Form("h_%d", current_task), "", nbins, -0.5, nbins - 0.5);
    float mean_dp = 0;
    for (uint i = 0; i < (uint)nbins; i++)
    {
        mean_dp += _dp_vec[i];
        h.SetBinContent((int)i + 1, _dp_vec[i]);
    }

    mean_dp = mean_dp / (float)nbins;

    uint left_bin   = 0;
    uint right_bin  = static_cast<uint>(nbins);
    uint center_bin = static_cast<uint>(h.GetMaximumBin());

    // find fist bin containing sub-mean-value LEFT from the peak
    for (uint i = center_bin; i >= 1; i--)
    {
        left_bin = i;
        if (_dp_vec[i - 1] - mean_dp < 0) { break; }
    }

    // find fist bin containing sub-mean-value RIGHT from the peak
    for (uint i = static_cast<uint>(center_bin); i <= (uint)nbins; i++)
    {
        right_bin = i;

        if (_dp_vec[i - 1] - mean_dp < 0) { break; }
    }

    // return values of result:
    // param     0 - reduced chi2
    // param     1 - sigma
    // param     2 - mu
    double *result = new double[3];
    // check if the fitting range is bigger that 5 bins, if its 5 bins or less return -1
    // as sigma
    if (right_bin - left_bin <= 5)
    {
        result[0] = -1;
        result[1] = -1;
        result[2] = -1;
    }

    // get axis values for given bins
    double low    = h.GetXaxis()->GetBinCenter(static_cast<int>(left_bin));
    double high   = h.GetXaxis()->GetBinCenter(static_cast<int>(right_bin));
    double middle = h.GetXaxis()->GetBinCenter(static_cast<int>(center_bin));

    // set fit function with parameters
    // param:   0 - amplitude
    // param:   1 - mu/center
    // param:   2 - sigma
    // param:   3 - bcg offset
    // param:   4 - bcg gain/linear
    TF1 fcn("gauss_with_background", TEC::gaussianWithBackGround, low, high, 5);
    fcn.SetNpx(1000);
    fcn.SetParLimits(0, 0., 1.);
    fcn.SetParameter(0, 0.2);
    fcn.SetParLimits(1, low, high);
    fcn.SetParameter(1, middle);
    fcn.SetParLimits(2, 0., 1.E6);
    fcn.SetParameter(2, 10.);
    fcn.SetParLimits(3, 0., 1.);
    fcn.SetParameter(3, mean_dp);
    // fcn.SetParLimits(4, 0, 10);
    fcn.SetParameter(4, 0);

    fcn.SetParName(0, "amplitude");
    fcn.SetParName(1, "mean");
    fcn.SetParName(2, "sigma");
    fcn.SetParName(3, "offset");
    fcn.SetParName(4, "gain/linear");

    h.Fit("gauss_with_background", "RQN");

    double chi2 = fcn.GetChisquare();
    rchi2       = chi2 / (double)(nbins);
    sigma       = fcn.GetParameter(2);
    mu          = fcn.GetParameter(1);
}
