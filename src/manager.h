#ifndef MANAGER_H
#define MANAGER_H

#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include "hntuple.h"


class Manager {
public:
    Manager(const char* name = "output.root"); // name of output file + init
    virtual ~Manager(); // finalize + write + deletes objects 

    TH1F* createHistogram(const TH1F* ptr);
    TH1F* createHistogram(const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup, Bool_t sumw2 = kTRUE);
    TH2F* createHistogram(const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup, Bool_t sumw2 = kTRUE);
    TH3F* createHistogram(const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup, Int_t nbinsz, Double_t zlow, Double_t zup, Bool_t sumw2 = kTRUE);
    // variable bin size:
    TH1F* createHistogram(const char* name, const char* title, Int_t nbinsx, const Double_t* xbins, Bool_t sumw2 = kTRUE);
    
    HNtuple* createNtuple(const char* name, const char* title = nullptr); // here only name and title c-tor

    static double openingangle(const TLorentzVector& a, const TLorentzVector& b);
    static void normalize(TH1 *hist);
    static TH1* signal(const char *name, TH1 *hist, TH1 *back1, TH1 *back2);

    virtual void initData() = 0;
    virtual void finalizeData() = 0;

private:
    void writeData();
    void writeHistograms();
    void writeNtuples();

    std::vector<TH1F*> histograms;
    std::vector<TH2F*> histograms2;
    std::vector<TH3F*> histograms3;
    std::vector<HNtuple*> ntuples;
    TFile* outfile {nullptr};
};


#endif // MANAGER_H

