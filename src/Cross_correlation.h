#pragma once

#include <algorithm> // copy
#include <atomic>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <string>
#include <vector>

// Root
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall"
#pragma GCC diagnostic ignored "-Wpedantic"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wdouble-promotion"
#pragma GCC diagnostic ignored "-Wfloat-equal"
#include <TF1.h>
#include <TGraph.h>
#include <TH1D.h>
#pragma GCC diagnostic pop

#include "RegionOfInterest.h"
#include "variables.h"

namespace TEC
{

#define SQRT_2_PI 2.5066282746310002

inline Double_t gaussianWithBackGround(const Double_t *x, const Double_t *par)
{
    return par[0] / (par[2] * SQRT_2_PI) *
               exp(-(((x[0] - par[1]) * (x[0] - par[1])) / (2 * par[2] * par[2]))) +
           par[3] + par[4] * x[0];
}

class CrossCorrel
{

  public:
    CrossCorrel() = delete;
    CrossCorrel(VarManager *v, ResCont **ResV) : V(v), ResVec(ResV){};
    void Process(unsigned int      thread_id,
                 std::atomic<int> *thread_task,
                 std::mutex       &mtx_task,
                 std::mutex       &mtx_fit);

    int    Normalize(std::vector<float> &v);
    double DotProduct(const std::vector<float> &v1, const std::vector<float> &v2);

  private:
    int current_task;

    double norm;
    double dp;
    int    retval;

    unsigned int fThread_id;

    // vector of dot products calculated using cross correlation
    std::vector<float> dp_vec;
    VarManager        *V;
    ResCont          **ResVec;

    void ROIAnalysis(const int ROI, const int time, std::mutex &mtx_fit);

    void SaveToContainer(const int time, const int ROI, std::mutex &mtx_fit);

    std::vector<float> GetDataVec(const int ROI, const int time);
    void               GetShift_Gaussian(double &rchi2, double &sigma, double &mu);
    void               GetShift_Poly2(double &shift, double &dp, std::mutex &mtx_fit);
};

} // namespace TEC