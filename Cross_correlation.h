#include <string>
#include <thread>
#include <atomic>
#include <mutex>
// Root
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

class CrossCorrel{

public:
	void Process(unsigned int thread_id, atomic<int>* thread_task, 
				VarManager *v, std::mutex &mtx, ResCont** ResV);
	void SetVariables(VarManager *v);
	void PassResVec(ResCont** rv);
	
private:
	int current_task;

	double norm;
	double dp;
	int retval;
	
	vector<double> dp_vec;
	VarManager *V;
	ResCont** ResVec;

	void ROIAnalysis(unsigned int ROI, int time, std::mutex &mtx);
	int Normalize(vector<double>& v);
	double DotProduct(vector<double>& v1, vector<double>& v2);

	void SaveDPtoFile(int time, int ROI, std::vector<double>& v);
	void SaveToContainer(int time, int ROI, std::vector<double>& v, double* arr);
	
	int IndexOfMaxValue(vector<double>& v);
	std::vector<double> GetDataVec(int ROI, int time);

	void SaveAsPNG(TH1D* h, char* name);
	void SaveAsPNG(TH1D* h, char* name, TF1* funct);
	double* FitSigma(std::vector<double> dp_vec, bool createPNG, int ROI, int time);
};
