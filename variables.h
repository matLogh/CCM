#include <iostream>
#include <vector>
#include <atomic>
#include "TH2D.h"

#ifndef HEADER_H
#define HEADER_H

extern std::atomic<int> thread_task;
extern std::atomic<bool> wait;

using namespace std;

struct Settings
{  
   //automatic settings
   bool use_container_flag = false;
   //user settings

   bool create_fit_PNG = false;
   bool create_distro_PNG = true;
   
   bool use_shift = true;
   bool use_mi = false;

   int SigmaDistro_bins = 1100;
   int SigmaDistro_from = -1;
   int SigmaDistro_to   = 100;

   int dpDistro_bins = 100;
   int dpDistro_from = 0;
   int dpDistro_to   = 1;

   int selectFitFunct = 0; //0-autochoose to make spline, 1-zero cross, 2-linnear, 3-quadratic, 4-cubic
   bool allowLowOrderFit = true;
};

struct Rules
{  
   bool* isSigma;
   bool* isDP;
   bool* isChi2;
   double* sigma_low;
   double* sigma_high;
   double* dp;
   double* chi2;
};


struct VarManager
{
   //INPUT VIABLES
   TH2D* TEMAT;
   double*** TEMATarr;              //TH2D loaded as a data cube [ROI][x][y]
                                    //keep in mind that coord x and y starts with 0 for every ROI
   double time_width;

   int number_of_ROIs;
   int total_tasks;
   int time_bins;
   //EVERYTHING IS IN BINS!!!
   //dimensions of these std::vector are for different ROIs
   std::vector<uint> ewin_low;          //low edge of energy window, basically defines the sample vector beginning 
   std::vector<uint> ewin_high;         //high edge of energy window, basically defines the sample vector end 
   std::vector<uint> displ_low;         //lowest displacement edge 
   std::vector<uint> displ_high;        //highest displacement edge
   std::vector<uint> vector_dimension;  //number of dimension of vector
   std::vector<uint> base_shift_value;  //ewin_low[i] - displ_low[i] -> for purpose of displaying correct displacements
   std::vector<uint> displ_range;       //displ_high[i] - displ_low[i]

   std::vector<vector<double> > sample_vector;

   Settings setting;
   Rules rule;


};
//definition of result container
struct ResCont
{
   int shift;
   double dp;
   double fit_sigma;
   double fit_mi;
   double fit_chi2;
   bool isOk = true;
}; 

struct FitCont
{
   vector<double> coef;
   int functionUsed; //1-4
   int defaultFunction;
   double time_width;
};




#endif