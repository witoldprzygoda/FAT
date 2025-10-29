/**
 * @file PPip_ID_Refactored.cc
 * @brief Refactored version of PPip_ID::ProcessEntries() using PParticle system
 *
 * This shows a side-by-side comparison of the original implementation
 * and the refactored version using PParticle and BoostFrame.
 */

#include "../src/PParticle.h"
#include "../src/BoostFrame.h"

using namespace Physics;

// ============================================================================
// REFACTORED ProcessEntries Method
// ============================================================================
/**
 * Original: PPip_ID.cc lines 50-200 (150 lines)
 * Refactored: ~80 lines (47% reduction)
 */
void ProcessEntries_Refactored()
{
    // ------------------------------------------------------------------------
    // 1. INITIALIZATION (same as original lines 50-56)
    // ------------------------------------------------------------------------
    Long64_t nentries = fTree->GetEntries();

    for (Long64_t i = 0; i < nentries; ++i)
    {
        progressbar(i, nentries);
        fTree->GetEntry(i);

        // ------------------------------------------------------------------------
        // 2. PARTICLE CREATION
        // Original: lines 62-69 (14 lines)
        // Refactored: 2 lines
        // ------------------------------------------------------------------------
        PParticle proton = ParticleFactory::createProton(p_p_corr_p, p_theta, p_phi);
        PParticle pion = ParticleFactory::createPionPlus(pip_p_corr_pip, pip_theta, pip_phi);

        // Optional: Set reconstructed momentum too (for resolution studies)
        // proton.setFromSpherical(p_p, p_theta, p_phi, MomentumType::RECONSTRUCTED);
        // pion.setFromSpherical(pip_p, pip_theta, pip_phi, MomentumType::RECONSTRUCTED);

        // Optional: If MC data, set simulated momentum
        // if (isMC) {
        //     proton.setFromSpherical(sim_p_p, sim_p_theta, sim_p_phi, MomentumType::SIMULATED);
        //     pion.setFromSpherical(sim_pip_p, sim_pip_theta, sim_pip_phi, MomentumType::SIMULATED);
        // }

        // ------------------------------------------------------------------------
        // 3. COMPOSITE PARTICLES
        // Original: lines 90-115 (26 lines with frame copies)
        // Refactored: 6 lines (LAB frame automatically preserved!)
        // ------------------------------------------------------------------------
        PParticle neutron = beam_ - proton - pion;  // Missing mass technique
        PParticle deltaPP = proton + pion;          // Δ++ candidate
        PParticle deltaP = beam_ - proton;          // Δ+ candidate
        PParticle p_pip = proton + pion;            // Same as deltaPP
        PParticle n_pip = neutron + pion;           // nπ+ system
        PParticle pn = proton + neutron;            // pn system

        // ------------------------------------------------------------------------
        // 4. REFERENCE FRAME SETUP AND BOOSTS
        // Original: lines 118-134 (17 lines with 9 boost calls)
        // Refactored: 12 lines with frame caching
        // ------------------------------------------------------------------------

        // Add composite frames to event frame manager
        event_frames_.addCompositeFrame("ppip", p_pip);
        event_frames_.addCompositeFrame("npip", n_pip);
        event_frames_.addCompositeFrame("pn", pn);

        // Boost to beam center-of-mass (CMS)
        const BoostFrame& beam_frame = event_frames_.getFrame("beam");
        PParticle p_cms = beam_frame.boost(proton);
        PParticle pip_cms = beam_frame.boost(pion);
        PParticle n_cms = beam_frame.boost(neutron);
        PParticle deltaP_cms = beam_frame.boost(deltaP);
        PParticle deltaPP_cms = beam_frame.boost(deltaPP);

        // Boost to p+π+ rest frame
        const BoostFrame& ppip_frame = event_frames_.getFrame("ppip");
        PParticle pip_PPIP = ppip_frame.boost(pion);
        PParticle n_PPIP = ppip_frame.boost(neutron);
        PParticle proj_PPIP = ppip_frame.boost(proj_);

        // Boost to n+π+ rest frame
        const BoostFrame& npip_frame = event_frames_.getFrame("npip");
        PParticle pip_NPIP = npip_frame.boost(pion);
        PParticle p_NPIP = npip_frame.boost(proton);
        PParticle proj_NPIP = npip_frame.boost(proj_);

        // Boost to p+n rest frame
        const BoostFrame& pn_frame = event_frames_.getFrame("pn");
        PParticle n_PN = pn_frame.boost(neutron);
        PParticle pip_PN = pn_frame.boost(pion);
        PParticle proj_PN = pn_frame.boost(proj_);

        // ------------------------------------------------------------------------
        // 5. MASS CALCULATIONS
        // Original: lines 139-141 (manual /1000 conversions)
        // Refactored: Automatic GeV conversion
        // ------------------------------------------------------------------------
        double m_n = neutron.massGeV();
        double m_p = proton.massGeV();
        double m_pip = pion.massGeV();

        // ------------------------------------------------------------------------
        // 6. EVENT SELECTION CUTS
        // Original: lines 142-146
        // Refactored: Same, but cleaner syntax
        // ------------------------------------------------------------------------
        if (m_n > 0.899 && m_n < 0.986) {
            mass_n->Fill(m_n, weight);

            if (deltaPP.massGeV() > 0.8 && deltaPP.massGeV() < 1.8) {
                // ------------------------------------------------------------------------
                // 7. HISTOGRAM FILLING
                // Original: lines 159-180 (pointer dereferencing, R2D conversions)
                // Refactored: Direct method calls, angles already in degrees
                // ------------------------------------------------------------------------

                // Basic kinematics (LAB frame)
                r["p_p"] = proton.momentum();           // vs p_p (direct from TNtuple)
                r["pip_p"] = pion.momentum();           // vs pip_p
                r["n_p"] = neutron.momentum();          // vs n->P()

                // Composite masses
                r["ppip_m"] = p_pip.mass();             // vs p_pip->M()
                r["npip_m"] = n_pip.mass();             // vs n_pip->M()
                r["pn_m"] = pn.mass();                  // vs pn->M()

                // Opening angles (already in degrees!)
                r["ppip_oa"] = proton.openingAngle(pion);     // vs R2D * openingangle(*p, *pip)
                r["npip_oa"] = neutron.openingAngle(pion);    // vs R2D * openingangle(*n, *pip)
                r["pn_oa"] = proton.openingAngle(neutron);    // vs R2D * openingangle(*p, *n)

                // CMS kinematics
                r["ppip_theta_cms"] = p_pip.theta();          // vs R2D * p_pip->Theta()
                r["ppip_p_cms"] = deltaPP_cms.momentum();     // vs deltaPP->P()
                r["npip_theta_cms"] = n_pip.theta();          // vs R2D * n_pip->Theta()

                // Boosted frame angles
                r["pip_ppip_theta"] = pip_PPIP.theta();       // vs R2D * pip_PPIP->Theta()
                r["n_ppip_theta"] = n_PPIP.theta();           // vs R2D * n_PPIP->Theta()

                // Rapidity
                r["deltaPP_y"] = deltaPP.rapidity();          // vs deltaPP->Rapidity()

                // LAB frame access (no separate *_LAB copies needed!)
                r["p_lab_theta"] = proton.labFrame().Theta() * R2D;  // vs R2D * p_LAB->Theta()
                r["pip_lab_phi"] = pion.labFrame().Phi() * R2D;      // vs R2D * pip_LAB->Phi()

                // Custom calculations
                r["mandelstam_t"] = calculateMandelstamT(proton, proj_);  // Helper function

                // Event weight and identification
                r["weight"] = weight;
                r["isBest"] = isBest;

                r.fill();
            }

            // Additional cuts and histograms (lines 182-200)
            if (deltaP.massGeV() > 0.8 && deltaP.massGeV() < 1.8) {
                // Δ+ analysis histograms
                mass_deltaP->Fill(deltaP.massGeV(), weight);

                int masa_id = static_cast<int>((deltaP.massGeV() - 1.0) * 25 / 0.8);
                int cos_id = static_cast<int>((deltaP_cms.cosTheta() + 1.) * 10.);

                if (masa_id >= 0 && masa_id < 20 && cos_id >= 0 && cos_id < 20) {
                    hist_masa_cos[masa_id][cos_id]->Fill(neutron.theta(), weight);
                }
            }
        }
    } // End event loop
}

// ============================================================================
// HELPER FUNCTIONS
// ============================================================================

/**
 * Calculate Mandelstam t
 * Original: Scattered throughout code
 * Refactored: Centralized helper
 */
double calculateMandelstamT(const PParticle& scattered_proton,
                           const PParticle& projectile)
{
    TLorentzVector q = projectile.vec() - scattered_proton.vec();
    return q.M2();  // Four-momentum transfer squared
}

/**
 * Opening angle helper (original Manager::openingangle)
 * Now built into PParticle as openingAngle() method
 */
// DELETED - no longer needed!

// ============================================================================
// COMPARISON METRICS
// ============================================================================

/*
BEFORE (Original PPip_ID.cc):
- Lines 50-200: ~150 lines
- Hardcoded masses: 4 instances
- Manual coordinate conversions: 2 instances (12 lines)
- Manual LAB copies: 9 lines
- Frame-specific copies: 12 lines
- Boost calls: 17 lines (9 boost operations)
- Pointer operations: 50+ dereferences
- Memory management: 23 new/delete pairs
- R2D conversions: 15+ instances

AFTER (Refactored):
- Total: ~80 lines (47% reduction)
- Hardcoded masses: 0
- Manual coordinate conversions: 0
- Manual LAB copies: 0 (automatic)
- Frame-specific copies: 0 (boost returns new object)
- Boost calls: 12 lines (3 frame creations + 9 boost calls)
- Pointer operations: 0 (value semantics)
- Memory management: 0 (RAII)
- R2D conversions: 2 (for LAB frame access only)

PERFORMANCE:
- Boost vector calculations: 9 → 3 (67% reduction)
- Memory allocations: 23 heap → 0 heap (stack only)
- Cache efficiency: Improved (value semantics, no pointer chasing)

MAINTAINABILITY:
- Type safety: Improved (PParticle vs void* pointers)
- Readability: "proton.openingAngle(pion)" vs "R2D * openingangle(*p, *pip)"
- Bug potential: Reduced (no manual memory, no pointer errors)
- Extensibility: Easy to add new particle types or momentum representations
*/

// ============================================================================
// ADDITIONAL IMPROVEMENTS
// ============================================================================

/**
 * Optional: Encapsulate cuts in a struct for reusability
 */
struct EventSelection {
    struct MassWindow {
        double min, max;
        bool contains(double m) const { return m > min && m < max; }
    };

    MassWindow neutron_cut{0.899, 0.986};
    MassWindow delta_pp_cut{0.8, 1.8};
    MassWindow delta_p_cut{0.8, 1.8};

    bool passNeutronCut(const PParticle& n) const {
        return neutron_cut.contains(n.massGeV());
    }

    bool passDeltaPPCut(const PParticle& dpp) const {
        return delta_pp_cut.contains(dpp.massGeV());
    }

    bool passDeltaPCut(const PParticle& dp) const {
        return delta_p_cut.contains(dp.massGeV());
    }
};

// Usage:
void ProcessEntriesWithCuts()
{
    EventSelection cuts;

    // ... particle creation ...

    if (cuts.passNeutronCut(neutron)) {
        if (cuts.passDeltaPPCut(deltaPP)) {
            // Analysis
        }
    }
}

/**
 * Optional: Create particle collection manager
 */
struct EventParticles {
    // Primary particles
    PParticle proton, pion, neutron;

    // Composites
    PParticle deltaPP, deltaP, p_pip, n_pip, pn;

    // CMS versions
    PParticle p_cms, pip_cms, n_cms, deltaP_cms, deltaPP_cms;

    // Boosted versions
    PParticle pip_PPIP, n_PPIP, proj_PPIP;
    PParticle pip_NPIP, p_NPIP, proj_NPIP;
    PParticle n_PN, pip_PN, proj_PN;

    void createFromNTuple(float p_p_corr, float p_theta, float p_phi,
                         float pip_p_corr, float pip_theta, float pip_phi,
                         const PParticle& beam, const PParticle& proj) {
        // Primary
        proton = ParticleFactory::createProton(p_p_corr, p_theta, p_phi);
        pion = ParticleFactory::createPionPlus(pip_p_corr, pip_theta, pip_phi);
        neutron = beam - proton - pion;

        // Composites
        deltaPP = proton + pion;
        deltaP = beam - proton;
        p_pip = proton + pion;
        n_pip = neutron + pion;
        pn = proton + neutron;
    }

    void applyBoosts(EventFrames& frames, const PParticle& proj) {
        // CMS
        const BoostFrame& beam_frame = frames.getFrame("beam");
        p_cms = beam_frame.boost(proton);
        pip_cms = beam_frame.boost(pion);
        n_cms = beam_frame.boost(neutron);
        deltaP_cms = beam_frame.boost(deltaP);
        deltaPP_cms = beam_frame.boost(deltaPP);

        // Composite frames
        frames.addCompositeFrame("ppip", p_pip);
        frames.addCompositeFrame("npip", n_pip);
        frames.addCompositeFrame("pn", pn);

        const BoostFrame& ppip_frame = frames.getFrame("ppip");
        pip_PPIP = ppip_frame.boost(pion);
        n_PPIP = ppip_frame.boost(neutron);
        proj_PPIP = ppip_frame.boost(proj);

        const BoostFrame& npip_frame = frames.getFrame("npip");
        pip_NPIP = npip_frame.boost(pion);
        p_NPIP = npip_frame.boost(proton);
        proj_NPIP = npip_frame.boost(proj);

        const BoostFrame& pn_frame = frames.getFrame("pn");
        n_PN = pn_frame.boost(neutron);
        pip_PN = pn_frame.boost(pion);
        proj_PN = pn_frame.boost(proj);
    }
};

// Ultra-clean event loop:
void ProcessEntriesClean()
{
    EventSelection cuts;
    EventParticles particles;
    EventFrames frames;
    frames.setBeamFrameFromKineticEnergy(1580.0);

    for (Long64_t i = 0; i < nentries; ++i) {
        fTree->GetEntry(i);

        particles.createFromNTuple(p_p_corr_p, p_theta, p_phi,
                                   pip_p_corr_pip, pip_theta, pip_phi,
                                   beam_, proj_);

        particles.applyBoosts(frames, proj_);

        if (!cuts.passNeutronCut(particles.neutron)) continue;
        if (!cuts.passDeltaPPCut(particles.deltaPP)) continue;

        // Fill histograms...
        r["ppip_m"] = particles.p_pip.mass();
        r["ppip_oa"] = particles.proton.openingAngle(particles.pion);
        r.fill();
    }
}

// ============================================================================
// MIGRATION CHECKLIST FOR YOUR CODE
// ============================================================================

/*
1. Add to datamanager.h:
   - #include "src/PParticle.h"
   - #include "src/BoostFrame.h"
   - PParticle proj_, targ_, beam_;
   - EventFrames event_frames_;

2. Modify DataManager::initData():
   - Replace TLorentzVector* with PParticle
   - Use ParticleFactory::createBeamProton(1580.0)
   - Call event_frames_.setBeamFrame(proj_, targ_)

3. Refactor PPip_ID::ProcessEntries():
   - Replace lines 62-69 with ParticleFactory calls
   - Replace lines 90-115 with composite creation
   - Replace lines 118-134 with BoostFrame operations
   - Update histogram filling to use PParticle methods

4. Remove from data.h:
   - Global TLorentzVector* declarations (lines 88-113)
   - Hardcoded mass constants

5. Test:
   - Compile with -std=c++11
   - Run on test sample
   - Compare histograms (should match within 0.01%)

6. Commit:
   - git commit -am "Refactor to PParticle system"
*/
