#include <TFile.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TString.h>
#include <TSystem.h>

#include <iostream>

void combine_offset_gain_timeevo(int id_start = 0,
                                 int id_stop = 31,
                                 bool include_test = false)
{
    const char* input_dir = "../offset_and_gain";
    const char* gain_only_dir = "../gain_only";
    const char* test_dir = "../test";
    const char* output_dir = "output/combined_offset_and_gain_timeevo";

    const char* run_ranges[] = {
        "242083_243791",
        "243796_245219",
        "245222_245290",
        "246313_247251",
        "247253_248623"
    };

    const int run_first_all = 242083;
    const int run_last_all = 248623;
    const int n_x_bins = run_last_all - run_first_all + 1;

    const int y_bins_old = 32768;
    const int y_rebin = 2;
    const int y_bins_new = y_bins_old / y_rebin;

    gSystem->mkdir(output_dir, true);

    for (int id = id_start; id <= id_stop; ++id) {
        TString output_name = Form("%s/RunVsEnergy_all_ranges_det%d_correctedTimeEvo.root",
                                   output_dir, id);

        std::cout << "========================================" << std::endl;
        std::cout << "Detector Id" << id << std::endl;
        std::cout << "Output file: " << output_name << std::endl;

        if (!gSystem->AccessPathName(output_name)) {
            std::cout << "Overwriting existing output: " << output_name << std::endl;
        }

        TH2D* h_sum = new TH2D(Form("RunVsEnergy_det%d_all_ranges", id),
                               Form("Detector %d, all run ranges", id),
                               n_x_bins, run_first_all - 0.5, run_last_all + 0.5,
                               y_bins_new, 0.0, y_bins_old);

        TH2D* h_fixed_sum = new TH2D(Form("RunVsEnergy_det%d_fixed_all_ranges", id),
                                     Form("Detector %d fixed, all run ranges", id),
                                     n_x_bins, run_first_all - 0.5, run_last_all + 0.5,
                                     y_bins_new, 0.0, y_bins_old);

        TH2D* h_gain_only_sum = new TH2D(Form("RunVsEnergy_det%d_gain_only_all_ranges", id),
                                         Form("Detector %d gain-only fixed, all run ranges", id),
                                         n_x_bins, run_first_all - 0.5, run_last_all + 0.5,
                                         y_bins_new, 0.0, y_bins_old);

        TH2D* h_test_sum = nullptr;
        if (include_test) {
            h_test_sum = new TH2D(Form("RunVsEnergy_det%d_test_all_ranges", id),
                                  Form("Detector %d test fixed, all run ranges", id),
                                  n_x_bins, run_first_all - 0.5, run_last_all + 0.5,
                                  y_bins_new, 0.0, y_bins_old);
        }

        h_sum->GetYaxis()->SetRangeUser(3000.0, 4000.0);
        h_fixed_sum->GetYaxis()->SetRangeUser(3000.0, 4000.0);
        h_gain_only_sum->GetYaxis()->SetRangeUser(3000.0, 4000.0);
        if (h_test_sum) {
            h_test_sum->GetYaxis()->SetRangeUser(3000.0, 4000.0);
        }

        h_sum->GetXaxis()->SetNoExponent(true);
        h_fixed_sum->GetXaxis()->SetNoExponent(true);
        h_gain_only_sum->GetXaxis()->SetNoExponent(true);
        if (h_test_sum) {
            h_test_sum->GetXaxis()->SetNoExponent(true);
        }

        TGraph* g_roi0 = new TGraph();
        TGraph* g_roi1 = new TGraph();
        TGraph* g_roi0_gain_only = new TGraph();
        g_roi0->SetName(Form("shift_ROI_0_det%d_all_ranges", id));
        g_roi1->SetName(Form("shift_ROI_1_det%d_all_ranges", id));
        g_roi0_gain_only->SetName("ROI_0_gain_only");
        g_roi0->SetTitle("shifts for ROI 0, all run ranges");
        g_roi1->SetTitle("shifts for ROI 1, all run ranges");
        g_roi0_gain_only->SetTitle("gain-only shifts for ROI 0, all run ranges");

        int n_files_used = 0;
        int n_gain_only_files_used = 0;
        int n_test_files_used = 0;

        for (const char* runrange : run_ranges) {
            TString input_name = Form("%s/RunVsEnergy_%s_det%d_correctedTimeEvo.root",
                                      input_dir, runrange, id);
            TString gain_only_name = Form("%s/RunVsEnergy_%s_det%d_correctedTimeEvo.root",
                                          gain_only_dir, runrange, id);
            TString test_name = Form("%s/RunVsEnergy_%s_det%d_correctedTimeEvo.root",
                                     test_dir, runrange, id);

            std::cout << "Reading file: " << input_name << std::endl;

            TFile* f = TFile::Open(input_name, "READ");
            if (!f || f->IsZombie()) {
                std::cerr << "Skipping missing file: " << input_name << std::endl;
                continue;
            }

            TString h_name = Form("RunVsEnergy_det%d", id);
            TString h_fixed_name = Form("RunVsEnergy_det%d_fixed", id);
            TH2D* h = (TH2D*)f->Get(h_name);
            TH2D* h_fixed = (TH2D*)f->Get(h_fixed_name);
            TGraph* src_roi0 = (TGraph*)f->Get("shift_ROI_0");
            TGraph* src_roi1 = (TGraph*)f->Get("shift_ROI_1");

            std::cout << "  TH2D: " << h_name << std::endl;
            std::cout << "  TH2D: " << h_fixed_name << std::endl;
            std::cout << "  TGraph: shift_ROI_0" << std::endl;
            std::cout << "  TGraph: shift_ROI_1" << std::endl;

            if (!h || !h_fixed) {
                std::cerr << "  Missing TH2D, skipping file." << std::endl;
                f->Close();
                delete f;
                continue;
            }

            for (int ix = 1; ix <= h->GetNbinsX(); ++ix) {
                int run = static_cast<int>(h->GetXaxis()->GetBinCenter(ix) + 0.5);
                int ix_out = h_sum->GetXaxis()->FindBin(run);

                for (int iy = 1; iy <= h->GetNbinsY(); ++iy) {
                    double content = h->GetBinContent(ix, iy);
                    if (content == 0.0) {
                        continue;
                    }

                    int iy_out = ((iy - 1) / y_rebin) + 1;
                    h_sum->AddBinContent(h_sum->GetBin(ix_out, iy_out), content);
                }
            }

            for (int ix = 1; ix <= h_fixed->GetNbinsX(); ++ix) {
                int run = static_cast<int>(h_fixed->GetXaxis()->GetBinCenter(ix) + 0.5);
                int ix_out = h_fixed_sum->GetXaxis()->FindBin(run);

                for (int iy = 1; iy <= h_fixed->GetNbinsY(); ++iy) {
                    double content = h_fixed->GetBinContent(ix, iy);
                    if (content == 0.0) {
                        continue;
                    }

                    int iy_out = ((iy - 1) / y_rebin) + 1;
                    h_fixed_sum->AddBinContent(h_fixed_sum->GetBin(ix_out, iy_out), content);
                }
            }

            if (src_roi0) {
                for (int i = 0; i < src_roi0->GetN(); ++i) {
                    g_roi0->SetPoint(g_roi0->GetN(),
                                     src_roi0->GetX()[i],
                                     src_roi0->GetY()[i]);
                }
            } else {
                std::cerr << "  Missing shift_ROI_0" << std::endl;
            }

            if (src_roi1) {
                for (int i = 0; i < src_roi1->GetN(); ++i) {
                    g_roi1->SetPoint(g_roi1->GetN(),
                                     src_roi1->GetX()[i],
                                     src_roi1->GetY()[i]);
                }
            } else {
                std::cerr << "  Missing shift_ROI_1" << std::endl;
            }

            ++n_files_used;

            f->Close();
            delete f;

            std::cout << "Reading gain-only file: " << gain_only_name << std::endl;

            TFile* f_gain_only = TFile::Open(gain_only_name, "READ");
            if (!f_gain_only || f_gain_only->IsZombie()) {
                std::cerr << "Skipping missing gain-only file: " << gain_only_name << std::endl;
                if (f_gain_only) {
                    f_gain_only->Close();
                    delete f_gain_only;
                }
            } else {
                TH2D* h_gain_only_fixed = (TH2D*)f_gain_only->Get(h_fixed_name);
                TGraph* src_roi0_gain_only = (TGraph*)f_gain_only->Get("shift_ROI_0");

                std::cout << "  TH2D: " << h_fixed_name << std::endl;
                std::cout << "  TGraph: shift_ROI_0" << std::endl;

                if (!h_gain_only_fixed) {
                    std::cerr << "  Missing gain-only fixed TH2D, skipping file." << std::endl;
                } else {
                    for (int ix = 1; ix <= h_gain_only_fixed->GetNbinsX(); ++ix) {
                        int run = static_cast<int>(h_gain_only_fixed->GetXaxis()->GetBinCenter(ix) + 0.5);
                        int ix_out = h_gain_only_sum->GetXaxis()->FindBin(run);

                        for (int iy = 1; iy <= h_gain_only_fixed->GetNbinsY(); ++iy) {
                            double content = h_gain_only_fixed->GetBinContent(ix, iy);
                            if (content == 0.0) {
                                continue;
                            }

                            int iy_out = ((iy - 1) / y_rebin) + 1;
                            h_gain_only_sum->AddBinContent(h_gain_only_sum->GetBin(ix_out, iy_out), content);
                        }
                    }

                    if (src_roi0_gain_only) {
                        for (int i = 0; i < src_roi0_gain_only->GetN(); ++i) {
                            g_roi0_gain_only->SetPoint(g_roi0_gain_only->GetN(),
                                                       src_roi0_gain_only->GetX()[i],
                                                       src_roi0_gain_only->GetY()[i]);
                        }
                    } else {
                        std::cerr << "  Missing gain-only shift_ROI_0" << std::endl;
                    }

                    ++n_gain_only_files_used;
                }

                f_gain_only->Close();
                delete f_gain_only;
            }

            if (include_test) {
                std::cout << "Reading test file: " << test_name << std::endl;

                TFile* f_test = TFile::Open(test_name, "READ");
                if (!f_test || f_test->IsZombie()) {
                    std::cerr << "Skipping missing test file: " << test_name << std::endl;
                    if (f_test) {
                        f_test->Close();
                        delete f_test;
                    }
                    continue;
                }

                TH2D* h_test_fixed = (TH2D*)f_test->Get(h_fixed_name);

                std::cout << "  TH2D: " << h_fixed_name << std::endl;

                if (!h_test_fixed) {
                    std::cerr << "  Missing test fixed TH2D, skipping file." << std::endl;
                    f_test->Close();
                    delete f_test;
                    continue;
                }

                for (int ix = 1; ix <= h_test_fixed->GetNbinsX(); ++ix) {
                    int run = static_cast<int>(h_test_fixed->GetXaxis()->GetBinCenter(ix) + 0.5);
                    int ix_out = h_test_sum->GetXaxis()->FindBin(run);

                    for (int iy = 1; iy <= h_test_fixed->GetNbinsY(); ++iy) {
                        double content = h_test_fixed->GetBinContent(ix, iy);
                        if (content == 0.0) {
                            continue;
                        }

                        int iy_out = ((iy - 1) / y_rebin) + 1;
                        h_test_sum->AddBinContent(h_test_sum->GetBin(ix_out, iy_out), content);
                    }
                }

                ++n_test_files_used;

                f_test->Close();
                delete f_test;
            }
        }

        if (n_files_used == 0) {
            std::cerr << "No input files used for Id" << id << std::endl;
            delete h_sum;
            delete h_fixed_sum;
            delete h_gain_only_sum;
            delete h_test_sum;
            delete g_roi0;
            delete g_roi1;
            delete g_roi0_gain_only;
            continue;
        }

        TFile* out = TFile::Open(output_name, "RECREATE");
        if (!out || out->IsZombie()) {
            std::cerr << "Cannot create output ROOT file: " << output_name << std::endl;
            delete h_sum;
            delete h_fixed_sum;
            delete h_gain_only_sum;
            delete h_test_sum;
            delete g_roi0;
            delete g_roi1;
            delete g_roi0_gain_only;
            continue;
        }

        h_sum->Write();
        h_fixed_sum->Write();
        h_gain_only_sum->Write();
        if (h_test_sum) {
            h_test_sum->Write();
        }
        g_roi0->Write();
        g_roi1->Write();
        g_roi0_gain_only->Write();

        out->Close();
        delete out;

        std::cout << "Saved combined file: " << output_name << std::endl;
        std::cout << "Files used: " << n_files_used << std::endl;
        std::cout << "Gain-only files used: " << n_gain_only_files_used << std::endl;
        if (include_test) {
            std::cout << "Test files used: " << n_test_files_used << std::endl;
        }
        std::cout << "TH2D bins: x=" << n_x_bins
                  << ", y=" << y_bins_new
                  << " (rebinned by " << y_rebin << ")" << std::endl;
        std::cout << "Y-axis display range: 3000 to 4000" << std::endl;
        std::cout << "X-axis exponent disabled" << std::endl;
        std::cout << "shift_ROI_0 points: " << g_roi0->GetN() << std::endl;
        std::cout << "shift_ROI_1 points: " << g_roi1->GetN() << std::endl;
        std::cout << "ROI_0_gain_only points: " << g_roi0_gain_only->GetN() << std::endl;

        delete h_sum;
        delete h_fixed_sum;
        delete h_gain_only_sum;
        delete h_test_sum;
        delete g_roi0;
        delete g_roi1;
        delete g_roi0_gain_only;
    }
}
