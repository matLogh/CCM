#include <iostream>
#include <vector>

#include "variables.h"

using namespace std;

class CCM{
public:
	CCM();
	CCM(TH2D* matrix, vector<double>& rmh, vector<double>& rml, vector<double>& rp, vector<double>& rm, double model_time_low, double model_time_high);

	//void StartCCM();
	void StartCCM(unsigned int number_of_threads);
	VarManager* GetVar();
	void ChangeSampleTime(double model_time_low, double model_time_high);

	void AutoRules(int ROI, double sigma_width_acceptance = 3, double dp_width_acceptance = 3);
	void AutoRules_allROIs(double sigma_width_acceptance = 3, double dp_width_acceptance = 3);

	void SetRules_dp(int ROI, double dp);
	void SetRules_sigma(int ROI, double sigma_low, double sigma_high);
	void SetRules_chi2(int ROI, double chi2);

	void ApplyRules();

	void SetDrawDP(bool val);	//create dot product images

	void BuildShiftTable();
	void BuildFitTable();


private:


	unsigned int sample_time_start;
	unsigned int sample_time_end;
	unsigned int Xbins;
    unsigned int Ybins;

    atomic<int> thread_task;
	int number_of_threads;
	VarManager V;
	ResCont** ResVec;

	void DeleteSampleVectors();
	void HandleMatrix(TH2D* matrix);
	void CreateSampleVector(TH2D* matrix, int ROI_number, double ROI_time_low, double ROI_time_high);
	void Normalize(vector<double>& v);



};