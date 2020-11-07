
#ifndef CHECK_FUNCT
#define CHECK_FUNCT

#include <string>
#include <thread>
#include <atomic>
// Root
#include "TTree.h"
#include "TAxis.h"
#include "TH1D.h"
#include "TH2D.h"
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

using namespace std;
	
Double_t gaussian(Double_t *x, Double_t *par);
Double_t background(Double_t *x, Double_t *par);
Double_t gaussianWithBackGround(Double_t *x, Double_t *par);

double* FitSigma(std::vector<double> dp_vec, bool createPNG, int ROI, int time);
void SaveAsPNG(TH1D* h, char* name);
void SaveAsPNG(TH1D* h, char* name, TF1* funct);
void ApplyRls(VarManager *V, ResCont** rc);
void SetRules(int ROI, 
              VarManager *V,
              bool isSigma,
              double sigma_low,
              double sigma_high,
              bool isDP,
              double dp,
              bool isChi2,
              double chi2);
void SetAutoRules(int ROI, double sigma_width_acceptance, double dp_width_acceptance, VarManager *V, ResCont** ResVec);
double* SimpleFit(TH1D* h, char* name = NULL);
void BuildFitTable_fromContainer(VarManager *V, ResCont** ResVec);
void FixMatrix_fromFile(VarManager *V);
void FillFitContainer(VarManager *V, ResCont** ResVec, FitCont* FitVec);
void SaveFitTable(VarManager *V,FitCont* FitVec);
void FixTree_fromContainer(VarManager *V, TTree* tree, char* e_branch, char* t_branch, FitCont* FitVec, char* tree_name);
TH2D* GetEmptyClone(TH2D* mat, char* root_name);

#endif