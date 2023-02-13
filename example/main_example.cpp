#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <string>
#include <chrono> //measure time

// Root
#include "TTree.h"
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
#include "TChain.h"

#include "CCM.h"
#include "Cross_correlation.h"
#include "variables.h"
#include "CheckCCM.h"

using namespace std;
using namespace std::chrono;
//g++ -std=c++0x `root-config --libs --cflags`  CheckCCM.cpp CCM.cpp Cross_correlation.cpp main_piestany_fakeDecay.cpp -o CCM


#ifndef DATA_PATH
   #define DATA_PATH ""
#endif

int main(int argc, char **argv)  
{
   high_resolution_clock::time_point t1 = high_resolution_clock::now();

   //Get tree from rootfiles
   TChain* tree = new TChain("data_tree");
   // char str[100];
   for(int i=3; i<14; i++)
   {
      //this is a bit chaotic due to size limitation of the github repository
      // sprintf(str,"data/DecayGammaSpectroscopy_timeUnstable_%i.root",i);
      std::string path = DATA_PATH;
      path += "/data/DecayGammaSpectroscopy_timeUnstable_" + std::to_string(i) + ".root";
      std::cout << path << std::endl;

      tree->Add(path.c_str());
   }

   Double_t e;
   Long64_t t;
   tree->SetBranchAddress("e", &e);
   tree->SetBranchAddress("t", &t);

   //Fill time-energy matrix (TEMAT) with tree events
   //Binning of the TEMAT is crucial!
   int bins = 1000;
   double time_start = 22E12;
   double time_end =  50E12;
   cout << "creating matrix " <<endl;
   Long64_t nentries = tree->GetEntries();
   TH2D *TEMAT = new TH2D("dte", "distorted T vs E matrix",bins, time_start, time_end, 3200,0,1600);
   for(Long64_t entry = 0; entry < nentries; entry++)
   {
      tree->GetEntry(entry);
      TEMAT->Fill(t,e);
   }
   cout << "MATRIX DONE" <<endl;

   //Setup ROI
   vector<double> ewin_high;
   vector<double> ewin_low;
   vector<double> displ_low;
   vector<double> displ_high;

   //ROI energy limits
   ewin_low.push_back(2890);  //the 2890 is the binY value of the provided TEMAT
   ewin_high.push_back(2950); //the 2950 is the binY value of the provided TEMAT
   //Displacement range
   displ_low.push_back(-50);  // -50 BINS, not energy units
   displ_high.push_back(50);  // 50 BINS, not energy units

   //Set reference vector
   double scale_factor = (time_end - time_start)/(double)bins;
   int reference_spectrum_bgn = (3.919e13-time_start)/scale_factor;
   int reference_spectrum_end = (3.9468e13-time_start)/scale_factor;
   //reference_spectrum_bgn/reference_spectrum_end refers to binX values of provided TEMAT


   if (reference_spectrum_bgn >= reference_spectrum_end) 
   {
      cout << "Wrong definition of referemce vector" << endl;
      throw 1;
   }

   //create CCM object
   CCM fix(TEMAT, ewin_high, ewin_low, displ_high, displ_low, reference_spectrum_bgn, reference_spectrum_end);

   fix.SetDrawDP(false);
   fix.SelectFitFunction(1,1);
   //START CCM
   fix.StartCCM(8);
   //fix.AutoRules_allROIs(5,5); //optionally set rules for accepting fit results
   //fix.ApplyRules();           //if rules are set, apply them!

   fix.BuildShiftTable();        //create text file with final displacements of test spectra
   fix.BuildFitTable();          //create text file with final correction parameters
   fix.FixMatrix();              //create tresult.root file with original and corrected matrix
   fix.FixTree("fixtree.root", tree, "e","t"); //create new Tree with corrected events

   high_resolution_clock::time_point t2 = high_resolution_clock::now();
   auto duration = duration_cast<microseconds>( t2 - t1 ).count();
   cout << "TOTAL DURATION OF " << duration/(double)1E6 << " seconds" << endl;

   return 0;
}