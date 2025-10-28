#include "data.h"
#include "src/datamanager.h"

DataManager::DataManager(const char *name) : Manager(name)
{
    initData();
}

DataManager::~DataManager()
{
    finalizeData();
}

void DataManager::initData()
{
    // Beam parameters
    double proton_kinetic_energy = 1580; // MeV
    //double proton_kinetic_energy = 4500; // MeV
    double proton_mass = 938.27231;
    double proton_energy = proton_kinetic_energy + proton_mass;
    double proton_momentum = sqrt(proton_energy * proton_energy - proton_mass * proton_mass);

    // TLorentzVectors: beam is proj + targ
    proj = new TLorentzVector(0, 0, proton_momentum, proton_energy); // PROTON BEAM
    targ = new TLorentzVector(0, 0, 0, 938.27231);                   // PROTON TARG
    beam = new TLorentzVector(0, 0, 0, 0);                           // proj + targ
    *beam = *proj + *targ;

    // Dummy objects which are filled in the event loop
    pip = new TLorentzVector(0, 0, 0, 0);
    pip_LAB = new TLorentzVector(0, 0, 0, 0);
    p = new TLorentzVector(0, 0, 0, 0);
    p_LAB = new TLorentzVector(0, 0, 0, 0);
    n = new TLorentzVector(0, 0, 0, 0);
    n_LAB = new TLorentzVector(0, 0, 0, 0);
    deltaP = new TLorentzVector(0, 0, 0, 0);
    deltaP_LAB = new TLorentzVector(0, 0, 0, 0);
    deltaPP = new TLorentzVector(0, 0, 0, 0);
    deltaPP_LAB = new TLorentzVector(0, 0, 0, 0);

    p_pip = new TLorentzVector(0, 0, 0, 0);
    n_pip = new TLorentzVector(0, 0, 0, 0);
    pn = new TLorentzVector(0, 0, 0, 0);
    p_NPIP = new TLorentzVector(0, 0, 0, 0);
    n_PPIP = new TLorentzVector(0, 0, 0, 0);
    pip_PN = new TLorentzVector(0, 0, 0, 0);
    n_PN = new TLorentzVector(0, 0, 0, 0);
    pip_PPIP = new TLorentzVector(0, 0, 0, 0);
    pip_NPIP = new TLorentzVector(0, 0, 0, 0);

    proj_PPIP = new TLorentzVector(0, 0, 0, 0);
    proj_NPIP = new TLorentzVector(0, 0, 0, 0);
    proj_PN = new TLorentzVector(0, 0, 0, 0);

    // Histograms and ntuple
    //hp_p = createHistogram("hp_p", "Histogram of p_p", 1000, 0., 2500.);
    //hp_theta = createHistogram("hp_theta", "Histogram of p_theta", 180, 0., 180.);
    //hp_phi = createHistogram("hp_phi", "Histogram of p_phi", 36, 0., 360.);
    //hp_theta_phi = createHistogram("hp_theta_phi", "Histogram of p_theta_phi", 180, 0, 180, 36, 0, 360);

    //banana_p = createHistogram("banana_p", "banana_p", 1000, 0, 2, 2500, 0, 2500);
    //banana_pip = createHistogram("banana_pip", "banana_pip", 1000, 0, 2, 1500, 0, 1500);

    // COS theta - INV MASS resonance PLANE
    exp_sig = createHistogram("exp_sig", "exp_sig", 50, 0.8, 1.8, 20, -1, 1);

    // A
    pwa_pip_costh = createHistogram("pwa_pip_costh", "pwa_pip_costh", 40, -1., 1.);
    pwa_p_costh = createHistogram("pwa_p_costh", "pwa_p_costh", 40, -1., 1.);
    pwa_n_costh = createHistogram("pwa_n_costh", "pwa_n_costh", 40, -1., 1.);

    // B
    pwa_pip_p = createHistogram("pwa_pip_p", "pip_p", 50, 0., 1.0);
    pwa_p_p = createHistogram("pwa_p_p", "p_p", 50, 0., 2.0);
    pwa_n_p = createHistogram("pwa_n_p", "n_p", 50, 0., 2.0);

    // C
    pwa_ppip_m = createHistogram("pwa_ppip_m", "pwa_ppip_m", 60, 1.0, 1.6);
    pwa_npip_m = createHistogram("pwa_npip_m", "pwa_npip_m", 60, 0.9, 1.5);
    pwa_pn_m = createHistogram("pwa_pn_m", "pwa_pn_m", 60, 1.7, 2.3);

    // D
    pwa_pip_helicity = createHistogram("pwa_pip_helicity", "pwa_pip_helicity", 40, -1., 1.);
    pwa_pipn_helicity = createHistogram("pwa_pipn_helicity", "pwa_pipn_helicity", 40, -1., 1.);
    pwa_n_helicity = createHistogram("pwa_n_helicity", "pwa_n_helicity", 40, -1., 1.);

    // E
    pwa_pip_gj = createHistogram("pwa_pip_gj", "pwa_pip_gj", 40, -1., 1.);
    pwa_pipn_gj = createHistogram("pwa_pipn_gj", "pwa_pipn_gj", 40, -1., 1.);
    pwa_n_gj = createHistogram("pwa_n_gj", "pwa_n_gj", 40, -1., 1.);

    // ----------

    mass_deltaPP = createHistogram("mass_deltaPP", "mass_deltaPP", 50, 0.8, 1.8);
    mass_deltaP = createHistogram("mass_deltaP", "mass_deltaP", 50, 0.8, 1.8);
    cos_theta_deltaPP = createHistogram("cos_theta_deltaPP", "cos_theta_deltaPP", 20, -1, 1);
    cos_theta_deltaP = createHistogram("cos_theta_deltaP", "cos_theta_deltaP", 20, -1, 1);
    cos_theta_p = createHistogram("cos_theta_p", "cos_theta_p", 20, -1, 1);
    cos_theta_n = createHistogram("cos_theta_n", "cos_theta_n", 20, -1, 1);
    mass_n = createHistogram("mass_n", "mass_n", 1000, 0.5, 1.5);
    mass_p = createHistogram("mass_p", "mass_p", 1000, 0.5, 1.5);
    mass_pip = createHistogram("mass_pip", "mass_pip", 1000, 0.0, 0.5);

    char name[100];
    char namep[100];

    for (int i = 0; i < 40; ++i)
    {
        sprintf(name, "mass_n_tab_%d", i);
        mass_n_tab[i] = createHistogram(name, name, 100, 0.5, 1.5);
        sprintf(name, "mass_nn_tab_%d", i);
        mass_nn_tab[i] = createHistogram(name, name, 100, 0.5, 1.5);
    }

    for (int i = 0; i < 100; ++i)
    {
        sprintf(namep, "mass_p_tab_%d", i);
        mass_p_tab[i] = createHistogram(namep, namep, 100, 0.5, 1.5);
        sprintf(namep, "mass_pp_tab_%d", i);
        mass_pp_tab[i] = createHistogram(namep, namep, 100, 0.5, 1.5);
    }

    char name2dim[100];
    for (int i = 0; i < 20; ++i)
    {
        for (int j = 0; j < 25; ++j)
        {
            sprintf(name2dim, "miss_mass_invmass_costh_%d_%d", j + 1, i + 1);
            miss_mass_tab[j][i] = createHistogram(name2dim, name2dim, 100, 0.5, 1.5);
        }
    }



    eVertX = createHistogram("eVertX", "eVertX", 200,-100,100);
    eVertY = createHistogram("eVertY", "eVertY", 200,-100,100);
    eVertZ = createHistogram("eVertZ", "eVertZ", 500,-400,100);

    // ---

    nt = createNtuple("nt", "slimmed physics ntuple");
}

void DataManager::finalizeData()
{

    /****************************************************************************************/
}
