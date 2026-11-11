// export_deck_data.C -- numbers and histograms of the deck, written to data/deck_numbers.json.
// Sources: plots/output/cb_study_sim_std_full.root (plots/cb_study_sim.C), the simulation file (legs for the
// beta_new vs q*p picture, parent statistics) and the simulation PID polygons v5.
// Run from this directory: root -l -b -q 'export_deck_data.C+'
#include <TFile.h>
#include <TH1.h>
#include <TH2D.h>
#include <TCutG.h>
#include <TTree.h>
#include <TMath.h>
#include <cmath>
#include <cstdio>
#include <map>
#include <set>
#include <array>
#include <string>
#include <vector>

static FILE* J = nullptr; static bool first = true;
static void key(const std::string& k) { fprintf(J, "%s\"%s\":", first ? "" : ",\n", k.c_str()); first = false; }
static void arr(const std::string& k, const std::vector<double>& v) { key(k); fputc('[', J); for (size_t i = 0; i < v.size(); ++i) fprintf(J, "%s%.10g", i ? "," : "", v[i]); fputc(']', J); }
static std::vector<double> cont(TH1* h) { std::vector<double> v; for (int b = 1; b <= h->GetNbinsX(); ++b) v.push_back(h->GetBinContent(b)); return v; }
static std::vector<double> errs(TH1* h) { std::vector<double> v; for (int b = 1; b <= h->GetNbinsX(); ++b) v.push_back(h->GetBinError(b)); return v; }

void export_deck_data(const char* study = "../../plots/output/cb_study_sim_std_full.root",
                      const char* simf = "/hdd3/przygoda/HADES/files/sim/smash_std_sim_epem.root",
                      const char* cuts = "/hdd3/przygoda/HADES/files/sim/pp45_pid_cuts_sim_betanew_v5.root",
                      const char* out = "data/deck_numbers.json")
{
    J = fopen(out, "w"); fputs("{\n", J);
    TFile fs(study);
    TH1* ref = (TH1*)fs.Get("Npm_all_all_oa0");
    std::vector<double> ed; for (int b = 1; b <= ref->GetNbinsX() + 1; ++b) ed.push_back(ref->GetXaxis()->GetBinLowEdge(b)); arr("edges", ed);
    for (const char* v : {"all", "best"}) for (double oa : {0., 2.}) {
        for (const char* n : {"Npm", "Npp", "Nmm", "Btrue", "Blep", "Bfake", "Sgen", "Sdalitz", "Sconv", "Same", "CBhat"}) {
            TString hn = Form("%s_all_%s_oa%g", n, v, oa); TH1* h = (TH1*)fs.Get(hn); if (!h) { printf("missing %s\n", hn.Data()); continue; } arr(hn.Data(), cont(h)); if (TString(n) == "CBhat") arr(std::string(hn.Data()) + "_err", errs(h)); }
        for (const char* p : {"pm", "pp", "mm"}) for (int c = 0; c < 15; ++c) { TString hn = Form("h_all_%s_oa%g_%s_cat%d", v, oa, p, c); TH1* h = (TH1*)fs.Get(hn); if (h) arr(hn.Data(), cont(h)); }
        for (const char* p : {"pm", "pp", "mm"}) for (int c = 2; c <= 5; ++c) { TString hn = Form("hc_all_%s_oa%g_%s_c%d", v, oa, p, c); TH1* h = (TH1*)fs.Get(hn); if (h) arr(hn.Data(), cont(h)); }
        for (int c = 2; c <= 5; ++c) { TString hn = Form("hcB_all_%s_oa%g_c%d", v, oa, c); TH1* h = (TH1*)fs.Get(hn); if (h) arr(hn.Data(), cont(h)); }
    }
    // PID polygons (e- mirrored in the plotting script)
    TFile fc(cuts);
    for (const char* q : {"ep", "em"}) for (const char* s : {"1sig", "3sig", "5sig"}) { TCutG* g = (TCutG*)fc.Get(Form("%s_cut_%s", q, s)); if (!g) { printf("no polygon %s %s\n", q, s); continue; }
        std::vector<double> x(g->GetX(), g->GetX() + g->GetN()), y(g->GetY(), g->GetY() + g->GetN()); arr(std::string("poly_") + q + "_" + s + "_x", x); arr(std::string("poly_") + q + "_" + s + "_y", y); }

    // legs of written e+e- hypotheses: beta_new vs q*p; fractions inside 1/3/5 sigma; parent statistics
    TFile f(simf); TTree* t = (TTree*)f.Get("EpEm_ID"); t->SetBranchStatus("*", 0);
    Float_t ev, vz, vx, vy, pa, tha, pha, bea, pb, thb, phb, beb, ida, idb, para, parb, proca, procb, ga, gb, gpa, gpb, svxa, svya, svza, svxb, svyb, svzb, trig;
    auto sb = [&](const char* n, Float_t* a) { t->SetBranchStatus(n, 1); t->SetBranchAddress(n, a); };
    sb("event", &ev); sb("eVertReco_x", &vx); sb("eVertReco_y", &vy); sb("eVertReco_z", &vz); sb("trigbit", &trig);
    sb("ep_p", &pa); sb("ep_theta", &tha); sb("ep_phi", &pha); sb("ep_beta_new", &bea); sb("ep_sim_id", &ida); sb("ep_sim_parentid", &para); sb("ep_sim_processid", &proca); sb("ep_sim_geninfo", &ga); sb("ep_sim_grandparentid", &gpa); sb("ep_sim_vertexx", &svxa); sb("ep_sim_vertexy", &svya); sb("ep_sim_vertexz", &svza);
    sb("em_p", &pb); sb("em_theta", &thb); sb("em_phi", &phb); sb("em_beta_new", &beb); sb("em_sim_id", &idb); sb("em_sim_parentid", &parb); sb("em_sim_processid", &procb); sb("em_sim_geninfo", &gb); sb("em_sim_grandparentid", &gpb); sb("em_sim_vertexx", &svxb); sb("em_sim_vertexy", &svyb); sb("em_sim_vertexz", &svzb);
    const int nx = 200, ny = 200; const double xlo = -1200, xhi = 1200, ylo = 0.75, yhi = 1.25;
    TH2D h2("h2", "", nx, xlo, xhi, ny, ylo, yhi); h2.SetDirectory(nullptr);
    TCutG* cut[2][3]; for (int q = 0; q < 2; ++q) { const char* sg[3] = {"1sig", "3sig", "5sig"}; for (int s = 0; s < 3; ++s) cut[q][s] = (TCutG*)fc.Get(Form("%s_cut_%s", q == 0 ? "ep" : "em", sg[s])); }
    double nLeg[2] = {0, 0}, inside[2][3] = {{0, 0, 0}, {0, 0, 0}}; std::map<std::array<int, 3>, double> tally;
    long nEvents = 0, nPairs = 0; Float_t ce = -1, cx = 0, cy = 0, cz = 0; bool have = false;
    std::set<std::array<float, 3>> P, M;
    auto addLeg = [&](int q, double p, double be) { h2.Fill(q == 0 ? p : -p, be); nLeg[q] += 1; for (int s = 0; s < 3; ++s) if (cut[q][s]->IsInside(p, be)) inside[q][s] += 1; };
    for (Long64_t i = 0; i < t->GetEntries(); ++i) { t->GetEntry(i);
        if (!have || ev != ce || vx != cx || vy != cy || vz != cz) { P.clear(); M.clear(); ce = ev; cx = vx; cy = vy; cz = vz; have = true; if (vz > -500) ++nEvents; }
        if (vz <= -500) continue; ++nPairs;
        if (P.insert({pa, tha, pha}).second) addLeg(0, pa, bea);
        if (M.insert({pb, thb, phb}).second) addLeg(1, pb, beb);
        if ((ida == 2 || ida == 3) && (idb == 2 || idb == 3) && ida != idb && svxa == svxb && svya == svyb && svza == svzb && para == parb && proca == procb) tally[{int(para), int(ga), int(gpa)}] += 1;
    }
    std::vector<double> z; for (int iy = 1; iy <= ny; ++iy) for (int ix = 1; ix <= nx; ++ix) z.push_back(h2.GetBinContent(ix, iy));
    arr("pid_h2", z); arr("pid_nx_ny", {double(nx), double(ny)}); arr("pid_x", {xlo, xhi}); arr("pid_y", {ylo, yhi});
    arr("legs", {nLeg[0], nLeg[1]}); arr("inside_ep", {inside[0][0], inside[0][1], inside[0][2]}); arr("inside_em", {inside[1][0], inside[1][1], inside[1][2]});
    arr("n_events_pairs", {double(nEvents), double(nPairs)});
    std::vector<double> tall; for (auto& kv : tally) { tall.push_back(kv.first[0]); tall.push_back(kv.first[1]); tall.push_back(kv.first[2]); tall.push_back(kv.second); } arr("genuine_tally", tall);
    fputs("\n}\n", J); fclose(J); printf("written %s; events %ld, pairs %ld\n", out, nEvents, nPairs);
}
