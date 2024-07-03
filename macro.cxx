
// void macro()
{
    TFile *ff = TFile::Open("prisma_abberations.root");
    TH2D *mat = (TH2D *)ff->Get("CCM_matrix");

    TH1D *proj = mat->ProjectionY("proj", 770, 775);
    proj->Draw();
}