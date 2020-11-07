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
#include <bits/stdc++.h> 
#include <iostream> 
#include <sys/stat.h> 
#include <sys/types.h> 

// Root
#include "TTree.h"
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

#include "CCM.h"
#include "variables.h"
#include "Cross_correlation.h"
#include "CheckCCM.h"

using namespace std;

CCM::CCM(){
}

CCM::CCM(TH2D* matrix, vector<double>& reh, vector<double>& rel, vector<double>& rp, vector<double>& rm, double reference_time_low, double reference_time_high)
{
  
  if((reh.size() != rel.size()) || (rp.size() != rm.size()) || (rp.size() != reh.size())) 
  {
    printf("vectors fed to the constructor are not the same size \n");
    throw 1;
  }
  sample_time_start = reference_time_low;
  sample_time_end = reference_time_high;

  //inicialization of vectors
  V.number_of_ROIs = reh.size();

  V.ewin_high.reserve(V.number_of_ROIs);
  V.ewin_low.reserve(V.number_of_ROIs);
  V.displ_high.reserve(V.number_of_ROIs);
  V.displ_low.reserve(V.number_of_ROIs);
  V.vector_dimension.reserve(V.number_of_ROIs);
  V.sample_vector.reserve(V.number_of_ROIs);

  //convert vectors to bins
  V.time_width = (reference_time_high - reference_time_low);
  V.time_width = 1;
  for(int ROI=0; ROI<V.number_of_ROIs; ROI++)
  {
    if(rp[ROI]<0 || rm[ROI]>0)
    {
      printf("Wrong sign for the displacement regions!! \n");
      throw 1;
    }
    V.ewin_high.push_back(reh[ROI]);
    V.ewin_low.push_back(rel[ROI]);
    V.displ_high.push_back(rp[ROI]+reh[ROI]);
    V.displ_low.push_back(rel[ROI]+rm[ROI]);
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
    CreateSampleVector(V.TEMAT, ROI, reference_time_low, reference_time_high);
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
  //reserve space for fit container - fit results and other variables
  FitVec = new FitCont[V.time_bins];
  for(int time=0; time<V.time_bins; time++)
  {
    FitVec[time].coef.resize(4);
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
    for(uint j=0; j<Xbins; j++)
    { 
      V.TEMATarr[ROI][j] = new double[V.displ_range[ROI]];
    }
  }

  for(int ROI=0; ROI<V.number_of_ROIs; ROI++)
  {
    for(uint x=0; x<Xbins; x++)
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
  //loop over energy
  for(int t = V.ewin_low[ROI_number], vector_iterator = 0; t < V.ewin_high[ROI_number]; t++, vector_iterator++)
  {
    vec[vector_iterator] = 0;
    //loop over time, since sample vector should be formed of more than 1 time-bin width of TEMAT
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
  for(uint i=0; i < v.size(); i++)
  {
    norm = norm + (v[i]*v[i]);
  }

  if(norm <=0) 
  {
    printf("\nNormalization ERROR, sample vector cannot cosist of zeros!\n");
    throw 2;
  }

  norm = sqrt(norm);
  for(uint i=0; i < v.size(); i++)
  {
    v[i] = v[i]/(double)norm;
  }
}

void CCM::StartCCM(unsigned int threads = 4)
{ 

  //create directories
  mkdir("CCM_files", 0777);
  mkdir("CCM_files/dpFit_images", 0777);
  mkdir("CCM_files/dpFit_root", 0777);
  mkdir("CCM_files/dot_products", 0777);

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

void CCM::ChangeReferenceTime(double reference_time_low, double reference_time_high)
{
  sample_time_start = reference_time_low;
  sample_time_end = reference_time_high;

  V.time_width = (reference_time_high - reference_time_low);

  DeleteSampleVectors();

  for(int i=0; i< V.number_of_ROIs; i++){
    CreateSampleVector(V.TEMAT, i, reference_time_low, reference_time_high);
    cout << "sample VECTOR SIZE " << V.sample_vector.size() << endl;
  }
}

void CCM::DeleteSampleVectors()
{
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
  if(V.number_of_ROIs == 0)
  {
    cout << "ERROR: 0 ROIs added" << endl;
    return;
  } 
  cout << endl;
  cout << "************ BUILDING FIT TABLE ************" << endl;
  cout << "Fit Table format:" << endl;
  cout << "default_fit_function" <<endl;
  cout << "time_bin1  used_fit_fcn par[0]  par[1]  par[2] ..." << endl;
  cout << "time_bin2  used_fit_fcn par[0]  par[1]  par[2] ..." << endl;
  cout << "time_bin3  used_fit_fcn par[0]  par[1]  par[2] ..." << endl;
  cout << "etc" << endl;
  cout << endl; 
  if(V.setting.use_container_flag)
  { 
    //if(V.number_of_ROIs == 1 || V.setting.selectFitFunct == 1)
      //BuildFitTable_fromContainer_zeroCross(&V, ResVec);
    //else
     // BuildFitTable_fromContainer(&V, ResVec);
      FillFitContainer(&V, ResVec, FitVec);
      SaveFitTable(&V, FitVec);
  }
}

void CCM::SelectFitFunction(int selectFunction, bool allowLowOrderFit)
{ 
  if(selectFunction <0 || selectFunction > 4)
  { 
    cout << endl;
    cout << "****************** ERROR ******************" << endl;
    cout << "WRONG FUNCTION SLECTED!" << endl;
    cout << "Please select function in range 0-4" << endl;
    cout << "0 for autoselect (always spline)" << endl;
    cout << "1-4 for fit/spline or lower function order IF" << endl;
    cout << "IF allowLowOrderFit is true (default is false)" << endl << endl;
  }
  else
  {
    V.setting.selectFitFunct = selectFunction;
    V.setting.allowLowOrderFit = allowLowOrderFit;
  }
}

void CCM::FixMatrix()
{
  FixMatrix_fromFile(&V);
}

void CCM::FixTree(char* TFile_name, TTree* tree, char* e_branch, char* t_branch)
{
  FixTree_fromContainer(&V, tree, e_branch, t_branch, FitVec, TFile_name);
}