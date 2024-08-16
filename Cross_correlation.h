#pragma once

#include <algorithm> // copy
#include <atomic>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <string>
#include <vector>

// Root
#include "TApplication.h"
#include "TF1.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TSpline.h"

#include "variables.h"

namespace TEC
{

#define SQRT_2_PI 2.5066282746310002

inline Double_t gaussianWithBackGround(const Double_t *x, const Double_t *par)
{
    return par[0] / (par[2] * SQRT_2_PI) * exp(-(((x[0] - par[1]) * (x[0] - par[1])) / (2 * par[2] * par[2]))) +
           par[3] + par[4] * x[0];
}

class CrossCorrel
{

  public:
    CrossCorrel() = delete;
    CrossCorrel(VarManager *v, ResCont **ResV) : V(v), ResVec(ResV){};
    void Process(unsigned int thread_id, std::atomic<int> *thread_task, std::mutex &mtx_task, std::mutex &mtx_fit);

  private:
    int current_task;

    double norm;
    double dp;
    int retval;

    std::vector<double> dp_vec;
    VarManager *V;
    ResCont **ResVec;

    void ROIAnalysis(const int ROI, const int time, std::mutex &mtx_fit);
    int Normalize(std::vector<double> &v);
    double DotProduct(std::vector<double> &v1, std::vector<double> &v2);

    void SaveToContainer(const int time, const int ROI, std::mutex &mtx_fit);

    std::vector<double> GetDataVec(const int ROI, const int time);
    void FitControlGaussian(double &rchi2, double &sigma, double &mu);
    void GetShift(double &shift, double &dp, std::mutex &mtx_fit);
};

} // namespace TEC