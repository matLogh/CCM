#include "CCMInterpolator.h"

#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TPad.h>
#include <iostream>
#include <vector>

using namespace TEC;

int main()
{
    TApplication app("App", nullptr, nullptr);
    app.SetReturnFromRun(true);

    std::vector<double> x{0, 1, 2, 3, 4, 5};
    std::vector<double> y{0, 1.1, 3.9, 9.2, 16.8, 23.1};

    CCMInterpolator interpolator(ROOT::Math::Interpolation::Type::kAKIMA, true);
    for (size_t i = 0; i < x.size(); ++i) { interpolator.AddPoint(x[i], y[i], true); }
    interpolator.EnableInterpolation();

    TMultiGraph *mg = new TMultiGraph();

    {
        std::vector<double> x_new;
        std::vector<double> y_new;

        for (double x = 0; x < 5; x += 0.1)
        {
            x_new.push_back(x);
            y_new.push_back(interpolator.Eval(x));
        }
        TGraph *gr = new TGraph(x_new.size(), x_new.data(), y_new.data());
        gr->SetName("Interpolated");
        gr->SetMarkerStyle(20);
        gr->SetMarkerColor(kBlue);
        gr->SetLineColor(kBlue);
        gr->SetLineWidth(2);
        gr->SetFillStyle(0);
        gr->SetDrawOption("AL");
        mg->Add(gr, "ALP");
        // gr->SetMarkerStyle(22);
    }

    {
        interpolator.SetType(ROOT::Math::Interpolation::Type::kCSPLINE_PERIODIC, true);
        std::vector<double> x_new;
        std::vector<double> y_new;

        for (double x = 0; x < 5; x += 0.1)
        {
            x_new.push_back(x);
            y_new.push_back(interpolator.Eval(x));
        }
        TGraph *gr = new TGraph(x_new.size(), x_new.data(), y_new.data());
        gr->SetName("Interpolated");
        gr->SetMarkerStyle(20);
        gr->SetMarkerColor(kOrange);
        gr->SetLineColor(kOrange);
        gr->SetLineWidth(2);
        gr->SetFillStyle(0);
        gr->SetDrawOption("AL");
        mg->Add(gr, "ALP");
        // gr->SetMarkerStyle(22);
    }

    TGraph *gr_data = new TGraph(x.size(), x.data(), y.data());
    gr_data->SetName("Data");
    gr_data->SetMarkerStyle(22);
    gr_data->SetMarkerColor(kRed);
    gr_data->SetLineColor(kRed);
    gr_data->SetLineWidth(2);
    gr_data->SetFillStyle(0);
    gr_data->SetDrawOption("AP");
    mg->Add(gr_data, "ALP");

    TCanvas *c = new TCanvas("c", "Interpolation Example", 800, 600);
    c->SetGrid();
    mg->SetTitle("Interpolation Example;X;Y");
    mg->Draw("A");

    app.Run();

    return 0;
}