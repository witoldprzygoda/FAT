//#include "data.h"
#include "manager.h"

Manager::Manager(const char *name) : outfile{new TFile(name, "RECREATE")}
{
}

TH1F *Manager::createHistogram(const TH1F* ptr)
{
    TH1F *histogram = const_cast<TH1F*>(ptr);
    histograms.push_back(histogram);
    return histogram;
}

TH1F *Manager::createHistogram(const char *name, const char *title, Int_t nbinsx, Double_t xlow, Double_t xup, Bool_t sumw2)
{
    TH1F *histogram = new TH1F(name, title, nbinsx, xlow, xup);
    if (sumw2)
        histogram->Sumw2();
    histograms.push_back(histogram);
    return histogram;
}

TH2F *Manager::createHistogram(const char *name, const char *title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup, Bool_t sumw2)
{
    TH2F *histogram = new TH2F(name, title, nbinsx, xlow, xup, nbinsy, ylow, yup);
    if (sumw2)
        histogram->Sumw2();
    histograms2.push_back(histogram);
    return histogram;
}

TH3F *Manager::createHistogram(const char *name, const char *title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup, Int_t nbinsz, Double_t zlow, Double_t zup, Bool_t sumw2)
{
    TH3F *histogram = new TH3F(name, title, nbinsx, xlow, xup, nbinsy, ylow, yup, nbinsz, zlow, zup);
    if (sumw2)
        histogram->Sumw2();
    histograms3.push_back(histogram);
    return histogram;
}

TH1F* Manager::createHistogram(const char* name, const char* title, Int_t nbinsx, const Double_t* xbins, Bool_t sumw2)
{
    TH1F *histogram = new TH1F(name, title, nbinsx, xbins);
    if (sumw2)
        histogram->Sumw2();
    histograms.push_back(histogram);
    return histogram;
}


HNtuple *Manager::createNtuple(const char *name, const char *title)
{
    HNtuple *ntuple = new HNtuple(name, title);
    ntuple->setFile(outfile);
    ntuples.push_back(ntuple);
    return ntuple;
}

void Manager::writeData()
{
    writeHistograms();
    writeNtuples();
}

void Manager::writeHistograms()
{
    outfile->cd();
    for (auto &histogram : histograms)
    {
        histogram->Write();
    }
    for (auto &histogram : histograms2)
    {
        histogram->Write();
    }
    for (auto &histogram : histograms3)
    {
        histogram->Write();
    }
}

void Manager::writeNtuples()
{
    outfile->cd();
    for (auto &ntuple : ntuples)
    {
        ntuple->Write();
    }
}

Manager::~Manager()
{
    writeData();
    for (auto &histogram : histograms)
    {
        delete histogram;
    }
    for (auto &histogram : histograms2)
    {
        delete histogram;
    }
    for (auto &histogram : histograms3)
    {
        delete histogram;
    }
    // for (auto& ntuple : ntuples) {
    //     delete ntuple;
    // }
    outfile->Close();
    delete outfile;
}

double Manager::openingangle(const TLorentzVector &a, const TLorentzVector &b) 
{
    return TMath::ACos((a.Px() * b.Px() + a.Py() * b.Py() + a.Pz() * b.Pz()) / (a.Vect().Mag() * b.Vect().Mag()));
}

void Manager::normalize(TH1 *hist)
{
    for (Int_t j = 1; j < hist->GetNbinsX() + 1; ++j)
    {
        hist->SetBinContent(j, hist->GetBinContent(j) / hist->GetBinWidth(j));
        //         hist->SetBinError( j, TMath::Sqrt( hist->GetBinContent(j) ) );
        hist->SetBinError(j, hist->GetBinError(j) / hist->GetBinWidth(j));
    }
}

TH1* Manager::signal(const char *name, TH1 *hist, TH1 *back1, TH1 *back2)
{
    TH1 *ptr = (TH1 *)hist->Clone(name);
    for (Int_t j = 1; j < hist->GetNbinsX() + 1; ++j)
    {
        //ptr->SetBinContent(j, hist->GetBinContent(j) - back1->GetBinContent(j) - back2->GetBinContent(j));
        ptr->SetBinContent(j, hist->GetBinContent(j) - 2*TMath::Sqrt(back1->GetBinContent(j)*back2->GetBinContent(j)));
        ptr->SetBinError(j, TMath::Sqrt(hist->GetBinError(j) * hist->GetBinError(j) +
                                        back1->GetBinError(j) * back1->GetBinError(j) +
                                        back2->GetBinError(j) * back2->GetBinError(j)));
    }

    return ptr;
}
