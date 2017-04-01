#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <chrono> //measure time

// Root
#include "TAxis.h"
#include "TFile.h"
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
#include "TVirtualPad.h"

#include "Process.h"
#include "Cross_correlation.h"
#include "variables.h"
#include "CheckCCM.h"

using namespace std;
using namespace std::chrono;

int main(int argc, char **argv)  
{

   TFile *file = new TFile("simulated_gerlach.root"); 
   //file->ls();
   TH2D* TEMAT = (TH2D*)file->Get("dTE");
   cout << "File OPENED" <<endl;

   vector<double> ewin_high;
   vector<double> ewin_low;
   vector<double> displ_low;
   vector<double> displ_high;

   //ROI 0
   ewin_low.push_back(5800);
   ewin_high.push_back(5900);
   displ_low.push_back(ewin_low[0] - 50);
   displ_high.push_back(ewin_high[0] + 50);

   //ROI 1
  /* ewin_low.push_back(2620);
   ewin_high.push_back(2700);
   displ_low.push_back(ewin_low[1] - 50);
   displ_high.push_back(ewin_high[1] + 50);
*/

   int model_vector_bgn = 4400;
   int model_vector_end = 4420;

   int time_mod_from = 4400;
   int time_mod_to = 4420;
   
   CCM fix(TEMAT, ewin_high, ewin_low, displ_high, displ_low, model_vector_bgn, model_vector_end);

   high_resolution_clock::time_point t1 = high_resolution_clock::now();
   fix.SetDrawDP(false);
   fix.StartCCM(4);
   fix.AutoRules_allROIs();
   fix.ApplyRules();
   fix.BuildShiftTable();
   fix.BuildFitTable();

   high_resolution_clock::time_point t2 = high_resolution_clock::now();
   auto duration = duration_cast<microseconds>( t2 - t1 ).count();
   cout << "TOTAL DURATION OF " << duration/(double)1E6 << " seconds" << endl;

   file->Close();
   return 0;
}