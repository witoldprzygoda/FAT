// export_extra.C -- additional inputs of the deck, written to data/deck_extra.json (no pass over the big ntuples).
// Sources (all from existing outputs): the simulation study, the mixing outputs of plots/k_mix_pools.C and
// plots/mixing_build_betanew.C, the data outputs of plots/isbest_fraction_betanew.C and plots/isbest_emulation_betanew.C,
// plots/isbest_choice_sim.C, and the simulation PID polygons v6.
// Run from this directory: root -l -b -q 'export_extra.C+'
#include <TCutG.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <cstdio>
#include <string>
#include <vector>

static FILE* J = nullptr; static bool first = true;
static void key(const std::string& k) { fprintf(J, "%s\"%s\":", first ? "" : ",\n", k.c_str()); first = false; }
static void arr(const std::string& k, const std::vector<double>& v) { key(k); fputc('[', J); for (size_t i = 0; i < v.size(); ++i) fprintf(J, "%s%.10g", i ? "," : "", v[i]); fputc(']', J); }
static void dump(TFile& f, const std::string& name, const std::string& as, bool withErr = false) {
    TH1* h = (TH1*)f.Get(name.c_str()); if (!h) { printf("MISSING %s in %s\n", name.c_str(), f.GetName()); return; }
    std::vector<double> v, e; for (int b = 1; b <= h->GetNbinsX(); ++b) { v.push_back(h->GetBinContent(b)); e.push_back(h->GetBinError(b)); }
    arr(as, v); if (withErr) arr(as + "_err", e);
}

static void dump2(TFile& f, const std::string& name, const std::string& as) {   // flattened: for each x bin all y bins
    TH2* h = (TH2*)f.Get(name.c_str()); if (!h) { printf("MISSING %s in %s\n", name.c_str(), f.GetName()); return; }
    std::vector<double> v; for (int i = 1; i <= h->GetNbinsX(); ++i) for (int j = 1; j <= h->GetNbinsY(); ++j) v.push_back(h->GetBinContent(i, j)); arr(as, v);
}

void export_extra(const char* out = "data/deck_extra.json")
{
    J = fopen(out, "w"); fputs("{\n", J);
    const std::string O = "../../plots/output/";
    {   TFile f((O + "cb_study_sim_std_full.root").c_str());
        for (const char* v : {"all", "PT3"}) { for (const char* p : {"pm", "pp", "mm"}) dump(f, Form("hmix_%s_oa0_%s", v, p), Form("simmix_%s_%s", v, p));
            for (const char* n : {"Btrue", "Npm", "Npp", "Nmm", "Sgen"}) dump(f, Form("%s_%s_all_oa0", n, v), Form("sim_%s_%s", n, v));
            dump(f, Form("CBhat_%s_all_oa0", v), Form("sim_CBhat_%s", v), true); }
        for (int c = 2; c <= 5; ++c) { for (const char* p : {"pm", "pp", "mm"}) dump(f, Form("hc_PT3_all_oa0_%s_c%d", p, c), Form("sc_%s_c%d", p, c)); dump(f, Form("hcB_PT3_all_oa0_c%d", c), Form("scB_c%d", c)); }
    }
    {   TFile f((O + "k_mix_pools_sim.root").c_str());
        for (const char* pool : {"all", "nge3"}) for (const char* s : {"all", "PT3"}) for (const char* p : {"pm", "pp", "mm"}) dump(f, Form("h_%s_%s_%s", pool, s, p), Form("kpool_sim_%s_%s_%s", pool, s, p));
        dump(f, "h_rpc_all", "rpc_sim_all"); dump(f, "h_rpc_PT3", "rpc_sim_PT3");
        for (const char* s : {"all", "PT3", "even", "odd"}) { for (const char* q : {"B", "S", "N", "pp", "mm"}) dump(f, Form("h_truth_%s_%s", s, q), Form("t_%s_%s", s, q));
            if (std::string(s) == "even" || std::string(s) == "odd") for (const char* pool : {"all", "nge3"}) for (const char* p : {"pm", "pp", "mm"}) dump(f, Form("h_%s_%s_%s", pool, s, p), Form("kpool_sim_%s_%s_%s", pool, s, p)); }
        for (const char* c : {"c2", "c3", "c4"}) for (const char* p : {"pm", "pp", "mm"}) dump(f, Form("h_%s_PT3_%s", c, p), Form("kclass_sim_%s_%s", c, p));
        for (const char* pool : {"all", "nge3"}) for (const char* p : {"pm", "pp", "mm"}) dump2(f, Form("h2_%s_%s", pool, p), Form("k2_sim_%s_%s", pool, p));
        for (const char* q : {"B", "S", "N", "pp", "mm"}) dump2(f, Form("h2_truth_%s", q), Form("t2_sim_%s", q));
    }
    {   TFile f((O + "k_mix_pools_data.root").c_str());
        for (const char* pool : {"all", "nge3"}) for (const char* p : {"pm", "pp", "mm"}) dump(f, Form("h_%s_PT3_%s", pool, p), Form("kpool_data_%s_%s", pool, p));
        dump(f, "h_rpc_PT3", "rpc_data_PT3");
        for (const char* pool : {"all", "nge3"}) for (const char* p : {"pm", "pp", "mm"}) dump(f, Form("h_%s_PT2_%s", pool, p), Form("kpool_data2_%s_%s", pool, p));
        for (const char* c : {"c2", "c3", "c4"}) for (const char* p : {"pm", "pp", "mm"}) dump(f, Form("h_%s_PT3_%s", c, p), Form("kclass_data_%s_%s", c, p));
        for (const char* pool : {"all", "nge3"}) for (const char* p : {"pm", "pp", "mm"}) dump2(f, Form("h2_%s_%s", pool, p), Form("k2_data_%s_%s", pool, p));
        for (const char* q : {"N", "pp", "mm"}) dump2(f, Form("h2_truth_%s", q), Form("t2_data_%s", q));
    }
    {   TFile f((O + "mixing_betanew.root").c_str());
        for (const char* p : {"pm", "pp", "mm"}) dump(f, Form("h_mixnc_%s_PT3", p), Form("datamix_%s", p));
        for (const char* p : {"pm", "pp", "mm"}) { dump(f, Form("h_se_%s_PT3_fine", p), Form("fine_all_%s", p)); dump(f, Form("h_se_%s_best_PT3_fine", p), Form("fine_best_%s", p)); }
    }
    {   TFile f((O + "isbest_fraction_betanew.root").c_str());
        for (const char* p : {"pm", "pp", "mm"}) for (const char* v : {"all", "best"}) dump(f, Form("h_%s_%s_PT3", v, p), Form("data_%s_%s", v, p));
    }
    {   TFile fp((O + "isbest_emu_vtx_pm.root").c_str()), fa((O + "isbest_emu_vtx_pp.root").c_str()), fm((O + "isbest_emu_vtx_mm.root").c_str());
        for (int c = 2; c <= 5; ++c) {
            dump(fp, Form("h_cross_pm_PT3_c%d", c), Form("emu_cross_c%d", c)); dump(fp, Form("h_crosswin_pm_PT3_c%d", c), Form("emu_win_c%d", c));
            dump(fp, Form("h_all_pm_PT3_c%d", c), Form("emu_pm_c%d", c)); dump(fa, Form("h_all_pp_PT3_c%d", c), Form("emu_pp_c%d", c)); dump(fm, Form("h_all_mm_PT3_c%d", c), Form("emu_mm_c%d", c));
            dump(fp, Form("h_cross_pm_PT3_c%d_fine", c), Form("fine_cross_c%d", c)); dump(fp, Form("h_crosswin_pm_PT3_c%d_fine", c), Form("fine_win_c%d", c)); dump(fa, Form("h_all_pp_PT3_c%d_fine", c), Form("fine_pp_c%d", c)); dump(fm, Form("h_all_mm_PT3_c%d_fine", c), Form("fine_mm_c%d", c));
        }
    }
    {   TFile f((O + "isbest_choice_sim.root").c_str());
        for (const char* t : {"gen", "lep", "fake"}) for (const char* c : {"n2", "n3"}) dump(f, Form("hD_%s_%s", t, c), Form("hD_%s_%s", t, c));
        arr("hD_edges", {-6, 6, 240}); dump(f, "h_win", "choice_win"); dump(f, "h_tot", "choice_tot"); dump(f, "h_rpc", "choice_rpc");
    }
    {   TFile f((O + "k_eps_projections.root").c_str());
        for (int w = 0; w < 2; ++w) for (const char* v : {"mm", "thcms", "cosh", "pt"}) {
            for (const char* p : {"all", "nge3"}) for (const char* q : {"pm", "pp", "mm"}) dump(f, Form("mix_w%d_%s_%s_%s", w, v, p, q), Form("mix_w%d_%s_%s_%s", w, v, p, q));
            for (const char* c : {"n3", "n4"}) for (const char* t : {"all", "best"}) dump(f, Form("eps_w%d_%s_%s_%s", w, v, c, t), Form("eps_w%d_%s_%s_%s", w, v, c, t)); }
    }
    {   TFile fp((O + "isbest_emu_vtx_pm.root").c_str()), fa((O + "isbest_emu_vtx_pp.root").c_str()), fm((O + "isbest_emu_vtx_mm.root").c_str());
        for (int c = 3; c <= 4; ++c) { dump(fp, Form("h_cross_pm_PT2_c%d", c), Form("emu2_cross_c%d", c)); dump(fp, Form("h_crosswin_pm_PT2_c%d", c), Form("emu2_win_c%d", c)); }
        for (int c = 3; c <= 4; ++c) for (const char* t : {"cross", "crosswin", "all", "best"}) { dump(fa, Form("h_%s_pp_PT3_c%d", t, c), Form("cl_pp_%s_c%d", t, c)); dump(fm, Form("h_%s_mm_PT3_c%d", t, c), Form("cl_mm_%s_c%d", t, c)); }
    }
    for (const char* m : {"data", "sim"}) { TFile f((O + "isbest_choice_tag_" + m + ".root").c_str()); dump(f, "h_win", std::string("tag_") + m + "_win"); dump(f, "h_tot", std::string("tag_") + m + "_tot"); }
    {   TFile f("/hdd3/przygoda/HADES/files/sim/pp45_pid_cuts_sim_betanew_v6.root");
        for (const char* q : {"ep", "em"}) for (const char* s : {"1sig", "3sig", "5sig"}) { TCutG* g = (TCutG*)f.Get(Form("%s_cut_%s", q, s)); if (!g) { printf("no polygon %s %s\n", q, s); continue; }
            arr(std::string("v6_") + q + "_" + s + "_x", std::vector<double>(g->GetX(), g->GetX() + g->GetN())); arr(std::string("v6_") + q + "_" + s + "_y", std::vector<double>(g->GetY(), g->GetY() + g->GetN())); }
    }
    fputs("\n}\n", J); fclose(J); printf("written %s\n", out);
}
