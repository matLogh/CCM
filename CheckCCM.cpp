


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
#include <numeric>      // std::accumulate

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

#include "CheckCCM.h"
#include "variables.h"

using namespace std;




Double_t gaussian(Double_t *x, Double_t *par) 
{
  return par[0]*exp(-(((x[0]-par[1])*(x[0]-par[1]))/(2*par[2]*par[2])));
}

Double_t background(Double_t *x, Double_t *par) 
{
  return par[0] + par[1]*x[0];
}

Double_t gaussianWithBackGround(Double_t *x, Double_t *par) 
{
  return gaussian(x,par) + background(x,&par[3]);
}

Double_t linear(Double_t *x, Double_t *par) 
{
  return par[0] + par[1]*x[0];
}

Double_t quad(Double_t *x, Double_t *par) 
{
  return linear(x,par) + par[2]*x[0]*x[0];
}

Double_t cubic(Double_t *x, Double_t *par) 
{
  return quad(x,par) + par[3]*x[0]*x[0]*x[0];
}

Double_t zero_cross(Double_t *x, Double_t *par) 
{
  return par[0]*x[0];
}
/*
double* FitSigma(std::vector<double> dp_vec, bool createPNG, int ROI, int time)
{
  int bins = dp_vec.size();

  TH1D* h = new TH1D("h", "h", bins, 0, 1);
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
  g->SetParameter(3, mean_dp-0.2);
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
*/
void SaveAsPNG(TH1D* h, char* name)
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

void SaveAsPNG(TH1D* h, char* name, TF1* funct)
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

void ApplyRls(VarManager *V, ResCont** rc)
{
  for(int ROI=0;  ROI<V->number_of_ROIs; ROI++)
  {
    for(int time=0; time<V->time_bins;time++)
    {
      //dot product check
      if(V->rule.isDP[ROI])
      {
        if(!(rc[ROI][time].dp >= V->rule.dp[ROI]))
        {
          rc[ROI][time].isOk = false;
          continue;
        }
      }
      //sigma interval check
      if(V->rule.isSigma[ROI])
      {
        if(!( (rc[ROI][time].fit_sigma >= V->rule.sigma_low[ROI]) && (rc[ROI][time].fit_sigma <= V->rule.sigma_high[ROI])))
        {
          rc[ROI][time].isOk = false;
          continue;
        }
      }
      //chi2 check
      if(V->rule.isChi2[ROI])
      {
        if(!(rc[ROI][time].fit_chi2 <= V->rule.chi2[ROI]))
        {
          rc[ROI][time].isOk = false;
          continue;
        }
      }
      rc[ROI][time].isOk = true;
    }
  }
    
  
}

void SetRules(int ROI, 
              VarManager *V,
              bool isSigma,
              double sigma_low,
              double sigma_high,
              bool isDP,
              double dp,
              bool isChi2,
              double chi2)
{
  if(isSigma)
  {
    V->rule.isSigma[ROI] = true;
    V->rule.sigma_low[ROI] = sigma_low;
    V->rule.sigma_high[ROI] = sigma_high;
  }
  if(isDP)
  {
    V->rule.isDP[ROI] = true;
    V->rule.dp[ROI] = dp;
  }
  if(isChi2)
  {
    V->rule.isChi2[ROI] = true;
    V->rule.chi2[ROI] = chi2;
  }
}

void SetAutoRules(int ROI, double sigma_width_acceptance, double dp_width_acceptance, VarManager *V, ResCont** ResVec)
{ 

  TH1D* h_s = new TH1D("h_s", "sigma distribution", V->setting.SigmaDistro_bins, 
                                                    V->setting.SigmaDistro_from, 
                                                    V->setting.SigmaDistro_to);
  TH1D* h_dp = new TH1D("h_dp", "dot product distribution", V->setting.dpDistro_bins, 
                                                            V->setting.dpDistro_from, 
                                                            V->setting.dpDistro_to);
  for(int time=0; time<V->time_bins; time++)
  { 
    cout << ResVec[ROI][time].shift << '\t' << ResVec[ROI][time].dp << '\t' << ResVec[ROI][time].fit_sigma << endl;
    h_dp->Fill(ResVec[ROI][time].dp);
    h_s->Fill(ResVec[ROI][time].fit_sigma);
  }

  ofstream out_s, out_dp;

  char str[50];
  sprintf(str, "CCM_files/sigma_distro_ROI_%d.dat", ROI);
  out_s.open(str);
  sprintf(str, "CCM_files/dp_distro_ROI_%d.dat", ROI);
  out_dp.open(str);
  for(int i=0; i<V->setting.SigmaDistro_bins; i++){
    out_s <<i << '\t' << h_s->GetBinContent(i) << '\n';
  }
  for(int i=0; i<V->setting.dpDistro_bins; i++){
    out_dp <<i << '\t' << h_s->GetBinContent(i) << '\n';
  }
  out_dp.close();
  out_s.close();

  //get fit for sigma and store result
  double* sigma_result;//[0] - chi2, [1] - sigma, [2] - mean
  sprintf(str, "sigma_distro");
  sigma_result = SimpleFit(h_s, str);
  sprintf(str, "CCM_files/sigma_distro_ROI_%d.png", ROI);
  SaveAsPNG(h_s, str);
  //get fit for DP and store result
  double* dp_result;//[0] - chi2, [1] - sigma, [2] - mean
  sprintf(str, "dp_distro");
  dp_result = SimpleFit(h_dp, str);
  sprintf(str, "CCM_files/dp_distro_ROI_%d.png", ROI);
  SaveAsPNG(h_dp, str);

  cout << "********** AUTOMATIC RULE SETTING **********" << endl << endl;
  cout << "Sigma fit: " << endl;
  cout << '\t' << "chi2" << '\t' << sigma_result[0] << endl;
  cout << '\t' << "sigma" << '\t' << sigma_result[1] << endl;
  cout << '\t' << "mi" << '\t' << sigma_result[2] << endl;

  cout << endl;

  cout << "Dot product fit: " << endl;
  cout << '\t' << "chi2" << '\t' << dp_result[0] << endl;
  cout << '\t' << "sigma" << '\t' << dp_result[1] << endl;
  cout << '\t' << "mi" << '\t' << dp_result[2] << endl;

  cout << endl;

  SetRules(ROI, V, true, 
              sigma_result[2] - (sigma_result[1] * sigma_width_acceptance * 0.5),
              sigma_result[2] + (sigma_result[1] * sigma_width_acceptance * 0.5),
              true,
              dp_result[2] - (dp_result[1] * dp_width_acceptance * 0.5),
              false,
              0);

  delete[] dp_result;
  delete[] sigma_result;

  delete h_s;
  delete h_dp;

  
}

double* SimpleFit(TH1D* h, char* name)
{
  TCanvas *c = new TCanvas("c","",200,10,1920,1080);  

  TH1* t = h->ShowBackground(20, "nocompton");
  TH1D *no_bcg = (TH1D*)h->Clone("h1");

  //fill no_bcg with vals substracted form bcground
  for(int i=1; i<h->GetNbinsX(); i++){
    double value = h->GetBinContent(i) - t->GetBinContent(i);
    if(value > 0) no_bcg->SetBinContent(i, value);
    else no_bcg->SetBinContent(i, 0);
  }

  int low_bin = 1;
  int high_bin = h->GetNbinsX();
  int middle_bin = no_bcg->GetMaximumBin();

  double lowered_amplitude = no_bcg->GetBinContent(no_bcg->GetMaximumBin())/(double)3;

  //find fist bin containing lowered amplitude LEFT from the peak
  for(int i=middle_bin; i >= 1; i--){
    if(no_bcg->GetBinContent(i) < lowered_amplitude){
      low_bin = middle_bin - (middle_bin - i)*3;
      break;
    }
  }

  //find fist bin containing lowered amplitude RIGHT from the peak
  for(int i=middle_bin; i <= h->GetNbinsX(); i++){
    if(no_bcg->GetBinContent(i) < lowered_amplitude){
      high_bin = (i - middle_bin)*3 + middle_bin;
      break;
    }
  } 

  double low = h->GetXaxis()->GetBinCenter(low_bin);
  double high = h->GetXaxis()->GetBinCenter(high_bin);
  double middle = h->GetXaxis()->GetBinCenter(middle_bin);

  TF1 *g = new TF1("g", gaussianWithBackGround, low, high, 5);
  g->SetParameter(0,0.2);
  g->SetParLimits(0,0,1E6);
  g->SetParLimits(1, low, high);
  g->SetParameter(1, middle);
  g->SetParLimits(2,0,1E6);
  g->SetParameter(2,10);
  g->SetParameter(3,0.5);
  g->SetParameter(4,0);
  g->SetParName(0, "amplitude");
  g->SetParName(1, "mean");
  g->SetParName(2, "sigma");
  g->SetParName(3, "offset");
  g->SetParName(4, "linear");

  no_bcg->Fit("g","RQ");

  //return values of result:
  //param     0 - reduced chi2
  //param     1 - sigma
  //param     2 - mi
  double* result = new double[3];

  double chi2 = g->GetChisquare();
  result[0] = chi2/(double)(low_bin - high_bin);
  result[1] = g->GetParameter(2);
  result[2] = g->GetParameter(1);

  char str[50];
  if(name)
  {
    sprintf(str, "CCM_files/fit_%s.png",name);
  } 
  else 
  {
    sprintf(str, "CCM_files/fit_.png");
  }
  SaveAsPNG(h, str, g);

  delete no_bcg;
  delete t;
  delete g;
  delete c;

  return result;
}

void BuildFitTable_fromContainer(VarManager *V, ResCont** ResVec)
{
  char s[50];
  sprintf(s, "fit_coef_table.dat");

  ofstream output;
  output.open(s);

  std::vector<double> ROI_center;
  output << V->number_of_ROIs << '\n';
  for(int j=0; j<V->number_of_ROIs; j++)
  {
    ROI_center.push_back(V->ewin_low[j] + (V->ewin_high[j]-V->ewin_low[j])/(double)2);
    output << ROI_center.at(j);
    output << '\t';
  }
  output << '\n';

  TF1 *g;
  if(V->number_of_ROIs == 1 || V->setting.selectFitFunct == 1)
    g = new TF1("g", zero_cross, ROI_center[0], ROI_center[V->number_of_ROIs], 1);
  else if (V->number_of_ROIs == 2 || V->setting.selectFitFunct == 2)
    g = new TF1("g", linear, ROI_center[0], ROI_center[V->number_of_ROIs], 2);
  else if (V->number_of_ROIs == 3 || V->setting.selectFitFunct == 3)
    g = new TF1("g", quad, ROI_center[0], ROI_center[V->number_of_ROIs], 3);
  else if (V->number_of_ROIs >= 4 || V->setting.selectFitFunct == 4)
    g = new TF1("g", cubic, ROI_center[0], ROI_center[V->number_of_ROIs], 4);

  int bad_counter;
  int counter;
  double x[V->number_of_ROIs];
  double y[V->number_of_ROIs];
  for(int time=0; time < V->time_bins; time++)
  {
    bad_counter = 0;
    counter = 0;
    for(int ROI=0; ROI<V->number_of_ROIs; ROI++)
    {
      if(!ResVec[ROI][time].isOk) bad_counter++; // get number of bad points
      else                                        //if OK, add points to the array
      { 
        x[counter] = ROI_center[ROI];
        y[counter] = ResVec[ROI][time].shift;
        counter++;
      }
    } 
    

    output << time; 
    if(bad_counter == 0)
    {
      TGraph* gr = new TGraph(V->number_of_ROIs,x,y);
      gr->Fit("g","RQ");
      for(int ROI=0; ROI<V->number_of_ROIs; ROI++)
        output << '\t' << g->GetParameter(ROI);
      output << '\n';
    }
    else if(bad_counter < V->number_of_ROIs)
    {
      if(V->setting.allowLowOrderFit)
      {
        TGraph* gr = new TGraph(V->number_of_ROIs-bad_counter,x,y);
        TF1 *f;
        if(V->number_of_ROIs-bad_counter       == 1 || V->setting.selectFitFunct == 1)
          f = new TF1("f", zero_cross, ROI_center[0], ROI_center[V->number_of_ROIs], 1);
        else if (V->number_of_ROIs-bad_counter == 2 || V->setting.selectFitFunct == 2)
          f = new TF1("f", linear, ROI_center[0],     ROI_center[V->number_of_ROIs], 2);
        else if (V->number_of_ROIs-bad_counter == 3 || V->setting.selectFitFunct == 3)
          f = new TF1("f", quad, ROI_center[0],       ROI_center[V->number_of_ROIs], 3);
        else if (V->number_of_ROIs-bad_counter >= 4 || V->setting.selectFitFunct == 4)
          f = new TF1("f", cubic, ROI_center[0],      ROI_center[V->number_of_ROIs], 4);

        gr->Fit("f","RQ");
        for(int ROI=0; ROI<V->number_of_ROIs-bad_counter; ROI++)
          output << '\t' << f->GetParameter(ROI);
        for(int ROI=V->number_of_ROIs-bad_counter; ROI<V->number_of_ROIs; ROI++)
          output << '\t' << 0;
        output << '\n';
        delete f;

      }
      else 
      {
        for(int ROI=0; ROI<V->number_of_ROIs; ROI++)
          output << '\t' << 0;
        output << '\n';
      }
    }
    else 
    {
      for(int ROI=0; ROI<V->number_of_ROIs; ROI++)
        output << '\t' << 0;
      output << '\n';
    }
  }
}

void BuildFitTable_fromContainer_zeroCross(VarManager *V, ResCont** ResVec)
{
  char s[50];
  sprintf(s, "fit_coef_table.dat");

  ofstream output;
  output.open(s);

  std::vector<double> ROI_center;
  output << V->number_of_ROIs << '\n';

  ROI_center.push_back(0);
  for(int j=0; j<V->number_of_ROIs; j++)
  {
    ROI_center.push_back(V->ewin_low[j] + (V->ewin_high[j]-V->ewin_low[j])/(double)2);
    output << ROI_center.at(j);
    output << '\t';
  }
  output << '\n';

  TF1 *g;
  g = new TF1("g", zero_cross, ROI_center[0], ROI_center[V->number_of_ROIs], 1);
 
  int bad_counter;
  int counter;
  double x[V->number_of_ROIs+1];
  double y[V->number_of_ROIs+1];
  for(int time=0; time < V->time_bins; time++)
  {
    bad_counter = 0;
    x[0] = 0;
    y[0] = 0;
    counter = 1;
    for(int ROI=0; ROI<V->number_of_ROIs; ROI++)
    {
      if(!ResVec[ROI][time].isOk) bad_counter++; // get number of bad points
      else                                        //if OK, add points to the array
      { 
        x[counter] = ROI_center[ROI];
        y[counter] = ResVec[ROI][time].shift;
        counter++;
      }
    } 
    
    output << time; 
    if(bad_counter <V->number_of_ROIs)
    {
      TGraph* gr = new TGraph(V->number_of_ROIs+1-bad_counter,x,y);  
      gr->Fit("g","RQ");
      for(int ROI=0; ROI<V->number_of_ROIs; ROI++)
        output << '\t' << g->GetParameter(ROI);
      output << '\n';
    }
    else 
    {
      for(int ROI=0; ROI<V->number_of_ROIs; ROI++)
        output << '\t' << 0;
      output << '\n';
    }
  }
}