#include <TFile.h>
#include <TCutG.h>

void make_cutg()
{
    // Przybliżenie czerwonej obwódki ze slajdu
    // X = missing mass [GeV/c^2]
    // Y = invariant mass [GeV/c^2]

    const int n = 9;

    double x[n] = {
        1.790, 1.790, 1.940, 2.580, 2.820, 2.950, 2.950, 1.790, 1.790
    };

    double y[n] = {
        0.360, 1.490, 1.490, 0.720, 0.570, 0.570, 0.360, 0.360, 0.360
    };

    TCutG *cut_red = new TCutG("cut_2d", n, x, y);

    // Dostosuj do nazw swoich zmiennych:
    cut_red->SetVarX("mm_pippimepem");
    cut_red->SetVarY("m_pippimepem");
    cut_red->SetTitle("Approximate red graphical cut from slide (GeV)");

    TFile *fout = new TFile("cut_2d.root", "RECREATE");
    cut_red->Write();
    fout->Close();

    printf("Saved cut to file: cut_red_slide_GeV.root\n");
    printf("Cut name: cut_red_slide_GeV\n");
}
