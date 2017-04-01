/*
list of throw exceptions:

1 = vectors fed to the constructor are not the same size
2 = normalization of sample vector error = its field of zeroes
3 = vectors in dot product are not the same size
*/
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>
#include <thread>
#include <atomic>
#include <stdexcept>      // std::out_of_range
#include <fstream>//open file
#include <stdio.h> //strcat
#include <algorithm>    // copy 
#include <fstream>
#include <mutex>

// Root
#include "TAxis.h"
#include "TImage.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TStyle.h"
#include "TPad.h"
#include "TROOT.h"
#include "TColor.h"
#include "TGFrame.h"
#include "TSystem.h"
#include "TVirtualPad.h"


#include "variables.h"
#include "Cross_correlation.h"
#include "CheckCCM.h"

using namespace std;


void CrossCorrel::SetVariables(VarManager *v){
	V = v;
}

void CrossCorrel::Process(unsigned int thread_id, atomic<int>* thread_task, 
                  VarManager *v, std::mutex &mtx, ResCont** ResV)
{
  ResVec = ResV;
  V = v;
  while(true)
  {
    mtx.lock();
    current_task = ++(*thread_task);
    mtx.unlock();

    if(current_task >= V->total_tasks) 
    {
      break;
    }
    printf("THREAD %d PROCESSING-TASK %d \n", thread_id, current_task);
    ROIAnalysis(current_task % V->number_of_ROIs, current_task/V->number_of_ROIs, mtx);
    printf("THREAD %d  FINISHED-TASK  %d \n", thread_id, current_task);
  }
  printf("\n");
}

/****************************************************
value -999 passed to dp_vec indicates error and this 
must be disregarded later
****************************************************/
void CrossCorrel::ROIAnalysis(unsigned int ROI, int time, std::mutex &mtx)
{
  //create data vector with size that includes whole area around the floating vector
  vector<double> data_vec = GetDataVec(ROI, time);
  //vector of dot products
  dp_vec.resize(V->displ_range[ROI] - V->vector_dimension[ROI]);
  //temporery vector holding current data
  std::vector<double> temp_data;
  //cycle through through the whole region, calculate cross correlation
  for(int shift=0; shift < V->displ_range[ROI] - V->vector_dimension[ROI]; shift++)
  {
    temp_data = std::vector<double>(&data_vec[shift], &data_vec[shift + V->vector_dimension[ROI]]);
    if(Normalize(temp_data) == 0)
      dp_vec[shift] = DotProduct(temp_data, V->sample_vector[ROI]);
    else
      dp_vec[shift] = -999; 
  }
  mtx.lock();
  double *arr = FitSigma(dp_vec, V->setting.create_fit_PNG, ROI, time);
  mtx.unlock();
  SaveToContainer(time, ROI, dp_vec, arr);
  SaveDPtoFile(time, ROI, dp_vec);

  delete arr;
}

std::vector<double> CrossCorrel::GetDataVec(int ROI, int time)
{
  vector<double> data_vec(&V->TEMATarr[ROI][time][0], &V->TEMATarr[ROI][time][V->displ_range[ROI]]);
  return data_vec;
}

/****************************************************
retval:
   0 = OK
  -1 = input vector is full zeroes
****************************************************/
int CrossCorrel::Normalize(std::vector<double>& v)
{
  retval = 0;
  norm = 0;
  for(int i=0; i < v.size(); i++){
    norm = norm + (v[i]*v[i]);
  }
  if(norm <=0) {
    norm = 1;
    retval = -1;
  }

  norm = sqrt(norm);
  for(int i=0; i < v.size(); i++){
    v[i] = v[i]/(double)norm;
  }
  return retval;
}

double CrossCorrel::DotProduct(std::vector<double>& v1, std::vector<double>& v2)
{ 
  dp = 0;
  if(v1.size() != v2.size()) 
  {
    cout << "WRONG VECTOR SIZE" << endl;
    throw 3;
    return -999;
  }
  for(int i=0; i < v1.size(); i++)
  {
    dp += v1[i]*v2[i];
  }
  return dp;
}

void CrossCorrel::SaveToContainer(int time, int ROI, std::vector<double>& dp_vec, double* arr)
{
  //var->finnal_shift_values[ROI][time] = VectorMaximumIndex(v);
  //  var->finnal_dp_values[ROI][time] = v[var->finnal_shift_values[ROI][time]];
  //  var->finnal_shift_values[ROI][time] = var->finnal_shift_values[ROI][time] - var->base_shift_value[ROI];

  ResVec[ROI][time].shift      = IndexOfMaxValue(dp_vec);
  ResVec[ROI][time].dp         = dp_vec[ResVec[ROI][time].shift];
  ResVec[ROI][time].shift      = ResVec[ROI][time].shift - V->base_shift_value[ROI];
  ResVec[ROI][time].fit_chi2   = arr[0];
  ResVec[ROI][time].fit_sigma  = arr[1];
  ResVec[ROI][time].fit_mi     = arr[2];
  printf("ROI: %d, time: %d, DP %f shift %d , base_shift %d \n", ROI, time, ResVec[ROI][time].dp, ResVec[ROI][time].shift, V->base_shift_value[ROI]);
}

void CrossCorrel::SaveDPtoFile(int time, int ROI, std::vector<double>& v)
{ 
    //write whole dot product array to the file
    ofstream output;
    char str[100];
    sprintf(str, "CCM_files/dot_products/DP_coef_ROI%i_T%i.dat",ROI, time);
    output.open(str);
   //output << "timeBin "<<time << '\n';
    for(int i=0; i<v.size(); i++) 
    {
      output << i - V->base_shift_value[ROI] <<"\t" << v[i] <<'\n';
    }
    output << '\n';
    output.close();

    ofstream f;
    char st[100];
    sprintf(st, "CCM_files/dot_products/MAX_ROI%i_T%i.dat", ROI, time);
    f.open(st);
    f << ResVec[ROI][time].shift;
    f.close();
}

int CrossCorrel::IndexOfMaxValue(std::vector<double>& v){
  return std::distance(v.begin(), std::max_element(v.begin(), v.end())); 
}

double* CrossCorrel::FitSigma(std::vector<double> dp_vec, bool createPNG, int ROI, int time)
{
  int bins = dp_vec.size();

  char str[50];
  sprintf(str, "h_%d", current_task);

  TH1D* h = new TH1D(str, str, bins, 0, 1);
  double mean_dp = 0;
  for(int i=0; i<bins; i++)
  {
    mean_dp += dp_vec[i];
    h->SetBinContent(i+1, dp_vec[i]);
  }
  mean_dp = mean_dp/(double)bins;
  int left_bin = 0;
  int right_bin= bins;
  int center_bin = h->GetMaximumBin();

  //find fist bin containing sub-mean-value LEFT from the peak
  for(int i=center_bin; i >= 1; i--)
  {
    if(dp_vec[i] - mean_dp < 0)
    {
      left_bin = i;
      break;
    }
  }

  //find fist bin containing sub-mean-value RIGHT from the peak
  for(int i=center_bin; i <= bins; i++)
  {
    if(dp_vec[i] - mean_dp < 0)
    {
      right_bin = i;
      break;
    }
  } 

  //return values of result:
  //param     0 - reduced chi2
  //param     1 - sigma
  //param     2 - mi
  double* result = new double[3];
  //check if the fitting range is bigger that 5 bins, if its 5 bins or less return -1 as sigma
  if(right_bin - left_bin <= 5)
  {
    result[0] = -1;
    result[1] = -1;
    result[2] = -1;
  }

  //get axis values for given bins
  double low = h->GetXaxis()->GetBinCenter(left_bin);
  double high = h->GetXaxis()->GetBinCenter(right_bin);
  double middle = h->GetXaxis()->GetBinCenter(center_bin);

  //set fit function with parameters
  //param:   0 - amplitude
  //param:   1 - mi/center
  //param:   2 - sigma
  //param:   3 - bcg offset
  //param:   4 - bcg gain/linear
  TF1 *g = new TF1("g", gaussianWithBackGround, low, high, 5);
  g->SetParLimits(0, 0, 1); 
  g->SetParameter(0, 0.2);
  g->SetParLimits(1, low, high);
  g->SetParameter(1, middle);
  g->SetParLimits(2, 0, 1E6);
  g->SetParameter(2, 10);
  g->SetParLimits(3, 0, 1);
  g->SetParameter(3, mean_dp);
  g->SetParLimits(4, 0, 10);
  g->SetParameter(4, 0);
  g->SetParName(0, "amplitude");
  g->SetParName(1, "mean");
  g->SetParName(2, "sigma");
  g->SetParName(3, "offset");
  g->SetParName(4, "gain/linear");

  h->Fit("g","RQ");

  if(createPNG)
  {
    char str[100];
    sprintf(str, "CCM_files/fit_images/DP_coef_ROI%i_T%i.png",ROI, time);
    SaveAsPNG(h, str, g);
  }

  double chi2 = g->GetChisquare();
  result[0] = chi2/(double)(right_bin - left_bin);
  result[1] = g->GetParameter(2);
  result[2] = g->GetParameter(1);

  delete h;
  delete g;
  return result;
}

void CrossCorrel::SaveAsPNG(TH1D* h, char* name)
{
  TCanvas *c1 = new TCanvas("c1","",200,10,1920,1080);    
  c1->SetFillColor(0);
  c1->SetGrid();
  h->Draw();

  TImage *img = TImage::Create();
  img->FromPad(c1);
  img->WriteImage(name);

  delete img;
  delete c1;
}

void CrossCorrel::SaveAsPNG(TH1D* h, char* name, TF1* funct)
{
  TCanvas *c1 = new TCanvas("c1","",200,10,1920,1080);    
  c1->SetFillColor(0);
  c1->SetGrid();
  h->Draw();
  funct->Draw("same");

  TImage *img = TImage::Create();
  img->FromPad(c1);

  img->WriteImage(name);

  delete img;
  delete c1;
}