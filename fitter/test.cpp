

#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "TApplication.h"
#include "TObject.h"
#include "TMath.h"
#include "TH1F.h"
#include "TF1.h"
#include "TColor.h"
#include "TCanvas.h"

#include "TheuerkaufPeak.hpp"

Double_t gauss(Double_t *x, Double_t *par)
{
	Double_t arg = 0;
	if (par[2] != 0) arg = (x[0] - par[1])/par[2];
	Double_t val = par[0]*TMath::Exp(-0.5*arg*arg);
	return val;
}


int main(int argc, char** argv)
{
    TApplication app("app", &argc, argv);

    TH1F *h1 = new TH1F("h1", "histo from a gaussian", 1200, 0, 200);
    TF1* fcn_gen = new TF1("fcn_gen", "gaus", 0, 200);
    fcn_gen->SetParameter(0,1);
    fcn_gen->SetParameter(1,30);
    fcn_gen->SetParameter(2,5);

    fcn_gen->SetLineColor(kGreen);
    // fcn_gen->Draw();
    

    auto fcn_bcg = std::make_unique<TF1>("fcn_bcg","pol2",0,200);
    fcn_bcg->SetParameter(0,30);
    fcn_bcg->SetParameter(1,3);
    fcn_bcg->SetParameter(2,0.4);
    h1->FillRandom("fcn_bcg",100000);

    auto fcn_bcg2 = std::make_unique<TF1>("fcn_bcg2","pol2",0,200);
    h1->Fit("fcn_bcg2","RQN");

    h1->FillRandom("fcn_gen",100000);
    fcn_gen->SetParameter(1,70);
    h1->FillRandom("fcn_gen",100000);
    fcn_gen->SetParameter(1,50);
    h1->FillRandom("fcn_gen",30000);

    h1->Sumw2();

    std::cout << "h11 " << h1->GetBinContent(50) << " " << h1->GetBinError(50) << " " << sqrt((double)h1->GetBinContent(50)) << std::endl; 

    TheuerkaufPeak pk(0,200,0,false,false,false);
    pk.SetParameter_Volume(sqrt(2.*M_PI),TheuerkaufPeak::ParamState::FREE,0, 1E9);
    pk.SetParameter_Position(150, TheuerkaufPeak::ParamState::FREE, 0,200);
    pk.SetParameter_Sigma(5., TheuerkaufPeak::ParamState::FREE,0,10);
    pk.SetParameter_TailLeft(5., TheuerkaufPeak::ParamState::FREE);

    TF1* tfcn = pk.GetFunction();
    h1->FillRandom(tfcn->GetName(), 10000);

    h1->Draw();








    // pk.SetParameter_TailLeft(2,TheuerkaufPeak::ParamState::FIXED,0, 1E9);
    // pk.SetParameter_StepHeight(-0.1,TheuerkaufPeak::ParamState::FIXED);
    // pk.SetParameter_StepWidth(1.,TheuerkaufPeak::ParamState::FIXED);

    // pk.Draw("SAME");

    


    // std::vector<TheuerkaufPeak> peak;
    // for(int i=0; i<20;i++)
    // {
    //     // peak.emplace_back(-10,10,i,false,false,false);
    //     auto& pkk = peak.emplace_back(0,100,0,false,false,false);

    //     pkk.SetParameter_Volume(sqrt(2.*M_PI),TheuerkaufPeak::ParamState::FREE,0, 1E9);
    //     pkk.SetParameter_Position(70, TheuerkaufPeak::ParamState::FREE, 0, 100);
    //     pkk.SetParameter_Sigma(1., TheuerkaufPeak::ParamState::FREE,0,10);
    //     // pk->SetParameter_TailLeft(i*0.1,TheuerkaufPeak::ParamState::FIXED,0, 1E9);
    //     // pk->SetParameter_TailRight(i*0.1,TheuerkaufPeak::ParamState::FIXED,0, 1E9);
    //     // pkk.SetParameter_StepHeight(-1.,TheuerkaufPeak::ParamState::FIXED);
    //     // pkk.SetParameter_StepWidth(i,TheuerkaufPeak::ParamState::FIXED);
    // }   

    // for(auto& pk : peak) pk.Draw("SAME");
    // pk.GetFunction()->GetYaxis()->SetRangeUser(-100,100);



    

    TheuerkaufFitter fitter(0, 200);
    // fitter.SetBacground(fcn_bcg2);
    fitter.SetBacgroundPoly(3);
    fitter.AddPeak(70,false,false,false);
    fitter.AddPeak(30,false,false,false);
    fitter.AddPeak(50,false,false,false);
    fitter.AddPeak(150,false,false,false);

    fitter.Fit(h1);
    fitter.Analyze(h1);

    fitter.GetPeak(0);


    app.Run();

    return 0;
}

void test()
{
    int argc;
    char** argv;
    main(argc, argv);
}