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

#include "Process.h"
#include "variables.h"
#include "Cross_correlation.h"
#include "CheckCCM.h"

using namespace std;

CCM::CCM(){
  printf("Naked class succesfully Initialized");
}

CCM::CCM(TH2D* matrix, vector<double>& rmh, vector<double>& rml, vector<double>& rp, vector<double>& rm, double sample_time_low, double sample_time_high)
{
  
  if((rmh.size() != rml.size()) || (rp.size() != rm.size()) || (rp.size() != rmh.size())) 
  {
    printf("vectors fed to the constructor are not the same size \n");
    throw 1;
  }
  sample_time_start = sample_time_low;
  sample_time_end = sample_time_high;

  //inicialization of vectors
  V.number_of_ROIs = rmh.size();

  V.ewin_high.reserve(V.number_of_ROIs);
  V.ewin_low.reserve(V.number_of_ROIs);
  V.displ_high.reserve(V.number_of_ROIs);
  V.displ_low.reserve(V.number_of_ROIs);
  V.vector_dimension.reserve(V.number_of_ROIs);
  V.sample_vector.reserve(V.number_of_ROIs);

  //convert vectors to bins
  V.time_width = (sample_time_high - sample_time_low);
  V.time_width = 1;
  for(int ROI=0; ROI<V.number_of_ROIs; ROI++)
  {
    V.ewin_high.push_back(rmh[ROI]);
    V.ewin_low.push_back(rml[ROI]);
    V.displ_high.push_back(rp[ROI]);
    V.displ_low.push_back(rm[ROI]);
  }

  //fill the displacement_range
  for(int ROI=0; ROI<V.number_of_ROIs; ROI++)
  {
    V.displ_range.push_back(V.displ_high[ROI] - V.displ_low[ROI]);  
  }

  HandleMatrix(matrix);
  

  //fill the vector_dimension
  for(int ROI=0; ROI<V.number_of_ROIs;ROI++)
  {
    V.vector_dimension.push_back(V.ewin_high[ROI] - V.ewin_low[ROI]);  
  }

  //create sample_vector
  for(int ROI=0; ROI<V.number_of_ROIs; ROI++)
  {
    CreateSampleVector(V.TEMAT, ROI, sample_time_low, sample_time_high);
  }
  
  //total_tasks = number_of_ROIs*Xbins/time_width;
  V.total_tasks = Xbins*V.number_of_ROIs/V.time_width;
  cout << "TOTAL OF " << V.total_tasks << " TASKS" << endl;

  //reserve space for result container structure
  ResVec = new ResCont*[V.number_of_ROIs];
  for(int ROI=0; ROI<V.number_of_ROIs; ROI++)
  {
    ResVec[ROI] = new ResCont[V.time_bins];
  }
  //calculate base shift vector
  for(int i=0; i< V.number_of_ROIs; i++){
    V.base_shift_value.push_back(V.ewin_low[i] - V.displ_low[i]);
  }

  //inicialize RULES substruct
  V.rule.sigma_low = new double[V.number_of_ROIs];
  V.rule.sigma_high = new double[V.number_of_ROIs];
  V.rule.dp = new double[V.number_of_ROIs];
  V.rule.chi2 = new double[V.number_of_ROIs];

  V.rule.isSigma = new bool[V.number_of_ROIs];
  V.rule.isDP = new bool[V.number_of_ROIs];
  V.rule.isChi2 = new bool[V.number_of_ROIs];
  for(int i=0; i<V.number_of_ROIs; i++)
  {
    V.rule.isSigma[i] = false;
    V.rule.isDP[i] = false;
    V.rule.isChi2[i] = false;
  }

  cout << "Class succesfully Initialized" << endl;
}

VarManager* CCM::GetVar()
{
  return &V;
}

void CCM::HandleMatrix(TH2D* matrix)
{
  V.TEMAT = matrix;
  Xbins = V.TEMAT->GetXaxis()->GetNbins();
  V.time_bins = Xbins;

  V.TEMATarr = new double**[V.number_of_ROIs];
  

  for(int ROI=0; ROI<V.number_of_ROIs; ROI++)
  {
    V.TEMATarr[ROI] = new double*[Xbins];
    for(int j=0; j<Xbins; j++)
    { 
      V.TEMATarr[ROI][j] = new double[V.displ_range[ROI]];
    }
  }

  for(int ROI=0; ROI<V.number_of_ROIs; ROI++)
  {
    for(int x=0; x<Xbins; x++)
    {
      for(int y=0; y<V.displ_range[ROI]; y++)
      {
          V.TEMATarr[ROI][x][y] = V.TEMAT->GetBinContent(x+1, y + V.displ_low[ROI]);
      }
    }
  }
}

void CCM::CreateSampleVector(TH2D* matrix, int ROI_number, double ROI_time_low, double ROI_time_high)
{
  std::vector<double> vec(V.vector_dimension[ROI_number]);
  
  int vector_iterator;
  for(int t = V.ewin_low[ROI_number], vector_iterator = 0; t < V.ewin_high[ROI_number]; t++, vector_iterator++)
  {
    vec[vector_iterator] = 0;
    for(int i=ROI_time_low; i<ROI_time_high; i++){
      vec[vector_iterator] += matrix->GetBinContent(i,t); 
    }
  }
  Normalize(vec);
  V.sample_vector.push_back(vec);
}

void CCM::Normalize(std::vector<double>& v)
{
  double norm = 0;
  for(int i=0; i < v.size(); i++)
  {
    norm = norm + (v[i]*v[i]);
  }

  if(norm <=0) 
  {
    printf("\nNormalization ERROR, sample vector cannot cosist of zeros!\n");
    throw 2;
  }

  norm = sqrt(norm);
  for(int i=0; i < v.size(); i++)
  {
    v[i] = v[i]/(double)norm;
  }
}

void CCM::StartCCM(unsigned int threads = 1)
{ 
  V.setting.use_container_flag = true;

  number_of_threads = threads;

  std::vector<std::thread> t;
  t.reserve(number_of_threads);

  printf("******----CCM starting----****** \n");
  
  thread_task = -1; //so that the first job start at zero

  CrossCorrel starting_class[number_of_threads];

  std::mutex mtx;
  //start threads
  if(number_of_threads > V.total_tasks)
  {
    for(int i=0; i<V.total_tasks; i++)
    {
      //pass thread task, VarManger variable, mutex and ResVec as reference
      t.push_back(std::thread(&CrossCorrel::Process, starting_class[i], i+1, &thread_task, &V, std::ref(mtx), ResVec));
    }
  }
  else
  {
    for(int i=0; i<number_of_threads; i++)
    {
      //pass thread task, VarManger variable, mutex and ResVec as reference
      t.push_back(std::thread(&CrossCorrel::Process, starting_class[i], i+1, &thread_task, &V, std::ref(mtx), ResVec));
    }
  }
  for (auto& th : t) th.join();
}

void CCM::ChangeSampleTime(double sample_time_low, double sample_time_high){
  sample_time_start = sample_time_low;
  sample_time_end = sample_time_high;

  V.time_width = (sample_time_high - sample_time_low);

  DeleteSampleVectors();

  for(int i=0; i< V.number_of_ROIs; i++){
    CreateSampleVector(V.TEMAT, i, sample_time_low, sample_time_high);
    cout << "sample VECTOR SIZE " << V.sample_vector.size() << endl;
  }
}

void CCM::DeleteSampleVectors(){
    V.sample_vector.clear();
}


void CCM::SetRules_dp(int ROI, double dp)
{
  SetRules(ROI, &V, 
          false, 0, 0, 
          true, dp,
          false, 0);
}


void CCM::SetRules_sigma(int ROI, double sigma_low, double sigma_high)
{
  SetRules(ROI, &V, 
          true, sigma_low, sigma_high, 
          false, 0,
          false, 0);
}


void CCM::SetRules_chi2(int ROI, double chi2)
{
  SetRules(ROI, &V, 
          false, 0, 0, 
          false, 0,
          true, chi2);
}

void CCM::SetDrawDP(bool val)
{
  V.setting.create_fit_PNG = val;
}

void CCM::BuildShiftTable()
{
  cout << endl;
  cout << "*********** BUILDING SHIFT TABLE ***********" << endl;
  cout << "Shift Table format:" << endl;
  cout << "number_of_ROIs" << endl;
  cout << "ROI1_center ROI2_center, ..." << endl;
  cout << "time_bin1 shift1OK(bool) shiftROI1 shift2OK(bool) shiftROI2 ..." << endl;
  cout << "time_bin2 shift1OK(bool) shiftROI1 shift2OK(bool) shiftROI2 ..." << endl;
  cout << "time_bin3 shift1OK(bool) shiftROI1 shift2OK(bool) shiftROI2 ..." << endl;
  cout << "etc" << endl;
  cout << endl;
  
  char name[] = "shift_table.txt";
  ofstream output;
  output.open(name);
  output << V.number_of_ROIs << '\n';
  for(int ROI=0; ROI<V.number_of_ROIs; ROI++)
  {
     output << V.ewin_low[ROI] + (V.ewin_high[ROI] - V.ewin_low[ROI])/(double)2 <<'\t';
  }
  output << '\n';

  for(int time=0; time < V.time_bins; time++)
  {
    output << time ;      
    for(int ROI=0; ROI<V.number_of_ROIs; ROI++){
      output << '\t' <<ResVec[ROI][time].isOk << '\t' << ResVec[ROI][time].shift;
    }
    output << '\n';
  }
}
void CCM::AutoRules(int ROI, double sigma_width_acceptance, double dp_width_acceptance)
{
  SetAutoRules(ROI, sigma_width_acceptance, dp_width_acceptance, &V, ResVec);
}

void CCM::AutoRules_allROIs(double sigma_width_acceptance, double dp_width_acceptance)
{
  for(int ROI =0; ROI<V.number_of_ROIs; ROI++)
  {
    SetAutoRules(ROI, sigma_width_acceptance, dp_width_acceptance, &V, ResVec);
  }
}
void CCM::ApplyRules()
{
  ApplyRls(&V, ResVec);
}

void CCM::BuildFitTable()
{ 
  cout << endl;
  cout << "************ BUILDING FIT TABLE ************" << endl;
  cout << "Fit Table format:" << endl;
  cout << "number_of_ROIs" <<endl;
  cout << "ROI1_center, ROI2_center, ..." << endl;
  cout << "time_bin1  par[0]  par[1]  par[2] ..." << endl;
  cout << "time_bin2  par[0]  par[1]  par[2] ..." << endl;
  cout << "time_bin3  par[0]  par[1]  par[2] ..." << endl;
  cout << "etc" << endl;
  cout << endl;
  if(V.setting.use_container_flag)
    BuildFitTable_fromContainer(&V, ResVec);
}