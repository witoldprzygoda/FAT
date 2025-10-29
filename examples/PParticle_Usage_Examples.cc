/**
 * @file PParticle_Usage_Examples.cc
 * @brief Comprehensive examples showing PParticle usage and migration from old code
 *
 * This file demonstrates:
 * 1. Before/After comparisons
 * 2. Common use cases
 * 3. Migration patterns
 */

#include "../src/PParticle.h"
#include "../src/BoostFrame.h"
#include "TLorentzVector.h"
#include <iostream>

using namespace Physics;

// ============================================================================
// EXAMPLE 1: Particle Creation - BEFORE vs AFTER
// ============================================================================

void example1_particle_creation() {
    std::cout << "\n=== EXAMPLE 1: Particle Creation ===" << std::endl;

    // -----------------------------------------------------------------------
    // OLD WAY (from PPip_ID.cc lines 62-69)
    // -----------------------------------------------------------------------
    std::cout << "\n--- OLD WAY ---" << std::endl;
    {
        // Hardcoded masses (error-prone!)
        const double proton_mass = 938.27231;
        const double pion_mass = 139.56995;
        const double D2R = 1.74532925199432955e-02;

        // Variables from TNtuple
        double p_p_corr_p = 1580.0;
        double p_theta = 45.0;
        double p_phi = 30.0;
        double pip_p_corr_pip = 850.0;
        double pip_theta = 60.0;
        double pip_phi = 120.0;

        // Manual spherical to Cartesian conversion
        TVector3 v0, v1;
        v0.SetXYZ(
            p_p_corr_p * sin(D2R * p_theta) * cos(D2R * p_phi),
            p_p_corr_p * sin(D2R * p_theta) * sin(D2R * p_phi),
            p_p_corr_p * cos(D2R * p_theta)
        );
        v1.SetXYZ(
            pip_p_corr_pip * sin(D2R * pip_theta) * cos(D2R * pip_phi),
            pip_p_corr_pip * sin(D2R * pip_theta) * sin(D2R * pip_phi),
            pip_p_corr_pip * cos(D2R * pip_theta)
        );

        // Create TLorentzVectors
        TLorentzVector* p = new TLorentzVector();
        TLorentzVector* pip = new TLorentzVector();
        p->SetVectM(v0, proton_mass);
        pip->SetVectM(v1, pion_mass);

        std::cout << "Proton:  E = " << p->E() << " MeV" << std::endl;
        std::cout << "Pion+:   E = " << pip->E() << " MeV" << std::endl;

        delete p;
        delete pip;
    }

    // -----------------------------------------------------------------------
    // NEW WAY (with PParticle)
    // -----------------------------------------------------------------------
    std::cout << "\n--- NEW WAY ---" << std::endl;
    {
        // Variables from TNtuple
        double p_p_corr_p = 1580.0;
        double p_theta = 45.0;
        double p_phi = 30.0;
        double pip_p_corr_pip = 850.0;
        double pip_theta = 60.0;
        double pip_phi = 120.0;

        // One-line creation with automatic mass and coordinate conversion
        PParticle proton = ParticleFactory::createProton(p_p_corr_p, p_theta, p_phi);
        PParticle pion = ParticleFactory::createPionPlus(pip_p_corr_pip, pip_theta, pip_phi);

        std::cout << "Proton:  E = " << proton.energy() << " MeV" << std::endl;
        std::cout << "Pion+:   E = " << pion.energy() << " MeV" << std::endl;
    }

    // BENEFITS:
    // ✓ No hardcoded masses
    // ✓ No manual coordinate conversion
    // ✓ No memory management (RAII)
    // ✓ Type-safe particle identification
    // ✓ 8 lines reduced to 2 lines
}

// ============================================================================
// EXAMPLE 2: Multiple Momentum Representations
// ============================================================================

void example2_momentum_variants() {
    std::cout << "\n=== EXAMPLE 2: Multiple Momentum Representations ===" << std::endl;

    // Create proton with three momentum types
    PParticle proton(MASS_PROTON, "p");

    // Set reconstructed momentum (raw detector data)
    proton.setFromSpherical(1550.0, 45.0, 30.0, MomentumType::RECONSTRUCTED);

    // Set corrected momentum (energy loss corrected)
    proton.setFromSpherical(1580.0, 45.0, 30.0, MomentumType::CORRECTED);

    // Set simulated momentum (MC truth, if available)
    proton.setFromSpherical(1575.0, 45.2, 30.1, MomentumType::SIMULATED);

    // Access any representation transparently
    std::cout << "Reconstructed E: " << proton.energy(MomentumType::RECONSTRUCTED) << " MeV" << std::endl;
    std::cout << "Corrected E:     " << proton.energy(MomentumType::CORRECTED) << " MeV" << std::endl;
    std::cout << "Simulated E:     " << proton.energy(MomentumType::SIMULATED) << " MeV" << std::endl;

    // Default is RECONSTRUCTED
    std::cout << "Default E:       " << proton.energy() << " MeV" << std::endl;

    // Use corrected for analysis (typical workflow)
    double analysis_momentum = proton.momentum(MomentumType::CORRECTED);
    std::cout << "Analysis momentum: " << analysis_momentum << " MeV/c" << std::endl;
}

// ============================================================================
// EXAMPLE 3: Composite Particles - BEFORE vs AFTER
// ============================================================================

void example3_composite_particles() {
    std::cout << "\n=== EXAMPLE 3: Composite Particles ===" << std::endl;

    // Create particles
    PParticle proton = ParticleFactory::createProton(1580.0, 45.0, 30.0);
    PParticle pion = ParticleFactory::createPionPlus(850.0, 60.0, 120.0);
    PParticle projectile = ParticleFactory::createBeamProton(1580.0);
    PParticle target = ParticleFactory::createTargetProton();

    // -----------------------------------------------------------------------
    // OLD WAY (from PPip_ID.cc lines 100-102)
    // -----------------------------------------------------------------------
    std::cout << "\n--- OLD WAY ---" << std::endl;
    {
        TLorentzVector* p = new TLorentzVector(proton.vec());
        TLorentzVector* pip = new TLorentzVector(pion.vec());
        TLorentzVector* beam = new TLorentzVector(projectile.vec() + target.vec());

        TLorentzVector* deltaPP = new TLorentzVector();
        TLorentzVector* n = new TLorentzVector();

        *deltaPP = *p + *pip;
        *n = *beam - *p - *pip;

        std::cout << "Delta++: M = " << deltaPP->M()/1000. << " GeV/c^2" << std::endl;
        std::cout << "Neutron: M = " << n->M()/1000. << " GeV/c^2" << std::endl;

        delete p; delete pip; delete beam; delete deltaPP; delete n;
    }

    // -----------------------------------------------------------------------
    // NEW WAY
    // -----------------------------------------------------------------------
    std::cout << "\n--- NEW WAY ---" << std::endl;
    {
        PParticle beam = projectile + target;
        PParticle deltaPP = proton + pion;
        PParticle neutron = beam - proton - pion;

        std::cout << "Delta++: M = " << deltaPP.massGeV() << " GeV/c^2" << std::endl;
        std::cout << "Neutron: M = " << neutron.massGeV() << " GeV/c^2" << std::endl;
    }

    // BENEFITS:
    // ✓ Natural mathematical syntax
    // ✓ No pointer dereferencing
    // ✓ No memory management
    // ✓ Automatic /1000 with massGeV()
    // ✓ 12 lines reduced to 3 lines
}

// ============================================================================
// EXAMPLE 4: Reference Frame Boosts - BEFORE vs AFTER
// ============================================================================

void example4_boosts() {
    std::cout << "\n=== EXAMPLE 4: Reference Frame Boosts ===" << std::endl;

    // Create particles
    PParticle proton = ParticleFactory::createProton(1580.0, 45.0, 30.0);
    PParticle pion = ParticleFactory::createPionPlus(850.0, 60.0, 120.0);
    PParticle projectile = ParticleFactory::createBeamProton(1580.0);
    PParticle target = ParticleFactory::createTargetProton();
    PParticle beam = projectile + target;

    // -----------------------------------------------------------------------
    // OLD WAY (from PPip_ID.cc lines 118-134)
    // -----------------------------------------------------------------------
    std::cout << "\n--- OLD WAY: Beam Rest Frame ---" << std::endl;
    {
        // Must copy before boosting
        TLorentzVector* p = new TLorentzVector(proton.vec());
        TLorentzVector* pip = new TLorentzVector(pion.vec());
        TLorentzVector* beam_copy = new TLorentzVector(beam.vec());

        // Manual boost call for each particle
        p->Boost(0.0, 0.0, -(*beam_copy).Beta());
        pip->Boost(0.0, 0.0, -(*beam_copy).Beta());

        std::cout << "Proton E_cms: " << p->E() << " MeV" << std::endl;
        std::cout << "Pion E_cms:   " << pip->E() << " MeV" << std::endl;

        delete p; delete pip; delete beam_copy;
    }

    // -----------------------------------------------------------------------
    // NEW WAY: Single Frame, Multiple Boosts
    // -----------------------------------------------------------------------
    std::cout << "\n--- NEW WAY: Beam Rest Frame ---" << std::endl;
    {
        BoostFrame beam_frame = BoostFrame::createBeamFrame(beam.beta());

        PParticle p_cms = beam_frame.boost(proton);
        PParticle pip_cms = beam_frame.boost(pion);

        std::cout << "Proton E_cms: " << p_cms.energy() << " MeV" << std::endl;
        std::cout << "Pion E_cms:   " << pip_cms.energy() << " MeV" << std::endl;
    }

    // -----------------------------------------------------------------------
    // OLD WAY: Composite Rest Frames (very verbose!)
    // -----------------------------------------------------------------------
    std::cout << "\n--- OLD WAY: Multiple Composite Frames ---" << std::endl;
    {
        // Create composite systems
        TLorentzVector* p_pip = new TLorentzVector(proton.vec() + pion.vec());
        TLorentzVector* neutron_local = new TLorentzVector((beam - proton - pion).vec());

        // Must copy before boosting
        TLorentzVector* pip_PPIP = new TLorentzVector(pion.vec());
        TLorentzVector* n_PPIP = new TLorentzVector(neutron_local->Px(),
                                                     neutron_local->Py(),
                                                     neutron_local->Pz(),
                                                     neutron_local->E());

        // Boost to p+pi+ frame (lines 125-127)
        pip_PPIP->Boost(-(*p_pip).BoostVector());
        n_PPIP->Boost(-(*p_pip).BoostVector());

        std::cout << "Pion in p+pi+ frame:   E = " << pip_PPIP->E() << " MeV" << std::endl;
        std::cout << "Neutron in p+pi+ frame: E = " << n_PPIP->E() << " MeV" << std::endl;

        delete p_pip; delete neutron_local; delete pip_PPIP; delete n_PPIP;
    }

    // -----------------------------------------------------------------------
    // NEW WAY: Composite Rest Frames
    // -----------------------------------------------------------------------
    std::cout << "\n--- NEW WAY: Multiple Composite Frames ---" << std::endl;
    {
        PParticle p_pip = proton + pion;
        PParticle neutron = beam - proton - pion;

        BoostFrame ppip_frame(p_pip);

        PParticle pip_boosted = ppip_frame.boost(pion);
        PParticle n_boosted = ppip_frame.boost(neutron);

        std::cout << "Pion in p+pi+ frame:   E = " << pip_boosted.energy() << " MeV" << std::endl;
        std::cout << "Neutron in p+pi+ frame: E = " << n_boosted.energy() << " MeV" << std::endl;
    }

    // BENEFITS:
    // ✓ Boost vector cached (calculated once)
    // ✓ Clear intent: "boost to ppip_frame"
    // ✓ No manual copying
    // ✓ 14 lines reduced to 5 lines
}

// ============================================================================
// EXAMPLE 5: EventFrames - Managing All Frames Together
// ============================================================================

void example5_event_frames() {
    std::cout << "\n=== EXAMPLE 5: EventFrames Manager ===" << std::endl;

    // Create particles
    PParticle proton = ParticleFactory::createProton(1580.0, 45.0, 30.0);
    PParticle pion = ParticleFactory::createPionPlus(850.0, 60.0, 120.0);
    PParticle projectile = ParticleFactory::createBeamProton(1580.0);
    PParticle target = ParticleFactory::createTargetProton();

    // Set up all reference frames at once
    EventFrames frames;
    frames.setBeamFrame(projectile, target);
    frames.addCompositeFrame("ppip", proton + pion);
    frames.addCompositeFrame("npip", (projectile + target - proton - pion) + pion);

    // Access any frame by name
    PParticle p_cms = frames.getFrame("beam").boost(proton);
    PParticle pip_ppip = frames.getFrame("ppip").boost(pion);

    std::cout << "Proton in CMS:      E = " << p_cms.energy() << " MeV" << std::endl;
    std::cout << "Pion in p+pi+ frame: E = " << pip_ppip.energy() << " MeV" << std::endl;

    // BENEFITS:
    // ✓ Centralized frame management
    // ✓ Named frames (self-documenting)
    // ✓ Easy to add new frames
}

// ============================================================================
// EXAMPLE 6: LAB Frame Preservation
// ============================================================================

void example6_lab_frame() {
    std::cout << "\n=== EXAMPLE 6: LAB Frame Preservation ===" << std::endl;

    PParticle proton = ParticleFactory::createProton(1580.0, 45.0, 30.0);

    std::cout << "LAB frame:  E = " << proton.energy() << " MeV" << std::endl;

    // Boost to another frame
    PParticle beam = ParticleFactory::createBeamProton(1580.0) +
                     ParticleFactory::createTargetProton();
    proton.boostZ(-beam.beta());

    std::cout << "After boost: E = " << proton.energy() << " MeV" << std::endl;

    // LAB frame still accessible
    std::cout << "LAB frame (preserved): E = " << proton.labFrame().E() << " MeV" << std::endl;

    // Reset to LAB
    proton.resetToLAB();
    std::cout << "After reset: E = " << proton.energy() << " MeV" << std::endl;

    // BENEFITS:
    // ✓ LAB frame always available
    // ✓ Easy frame switching
    // ✓ No manual LAB copies needed
}

// ============================================================================
// EXAMPLE 7: Complete Event Analysis Pattern
// ============================================================================

void example7_complete_event_analysis() {
    std::cout << "\n=== EXAMPLE 7: Complete Event Analysis ===" << std::endl;

    // -----------------------------------------------------------------------
    // Simulating TNtuple input
    // -----------------------------------------------------------------------
    struct NTupleData {
        float p_p_corr_p = 1580.0;
        float p_theta = 45.0;
        float p_phi = 30.0;
        float pip_p_corr_pip = 850.0;
        float pip_theta = 60.0;
        float pip_phi = 120.0;
    } ntuple;

    // -----------------------------------------------------------------------
    // Event reconstruction
    // -----------------------------------------------------------------------

    // Create particles from TNtuple
    PParticle proton = ParticleFactory::createProton(
        ntuple.p_p_corr_p, ntuple.p_theta, ntuple.p_phi);

    PParticle pion = ParticleFactory::createPionPlus(
        ntuple.pip_p_corr_pip, ntuple.pip_theta, ntuple.pip_phi);

    // Beam setup
    PParticle projectile = ParticleFactory::createBeamProton(1580.0);
    PParticle target = ParticleFactory::createTargetProton();
    PParticle beam = projectile + target;

    // Composite systems
    PParticle deltaPP = proton + pion;
    PParticle neutron = beam - proton - pion;

    // -----------------------------------------------------------------------
    // Apply cuts (natural syntax)
    // -----------------------------------------------------------------------
    bool pass_mass_cut = (neutron.massGeV() > 0.899 && neutron.massGeV() < 0.986);
    bool pass_deltapp_cut = (deltaPP.massGeV() > 0.8 && deltaPP.massGeV() < 1.8);

    std::cout << "Pass neutron mass cut: " << pass_mass_cut << std::endl;
    std::cout << "Pass Delta++ cut:      " << pass_deltapp_cut << std::endl;

    // -----------------------------------------------------------------------
    // Reference frame analysis
    // -----------------------------------------------------------------------
    EventFrames frames;
    frames.setBeamFrame(projectile, target);
    frames.addCompositeFrame("ppip", deltaPP);

    PParticle proton_cms = frames.getFrame("beam").boost(proton);
    PParticle pion_ppip = frames.getFrame("ppip").boost(pion);

    // -----------------------------------------------------------------------
    // Fill histograms (example)
    // -----------------------------------------------------------------------
    std::cout << "\nHistogram values:" << std::endl;
    std::cout << "  mass_n:       " << neutron.massGeV() << " GeV/c^2" << std::endl;
    std::cout << "  mass_deltaPP: " << deltaPP.massGeV() << " GeV/c^2" << std::endl;
    std::cout << "  cos_theta_cms: " << proton_cms.cosTheta() << std::endl;
    std::cout << "  opening_angle: " << proton.openingAngle(pion) << " deg" << std::endl;

    // BENEFITS: Complete workflow in ~30 clean, readable lines
}

// ============================================================================
// MAIN
// ============================================================================

int main() {
    std::cout << "=====================================================" << std::endl;
    std::cout << "    PParticle Usage Examples" << std::endl;
    std::cout << "    Before/After Migration Guide" << std::endl;
    std::cout << "=====================================================" << std::endl;

    example1_particle_creation();
    example2_momentum_variants();
    example3_composite_particles();
    example4_boosts();
    example5_event_frames();
    example6_lab_frame();
    example7_complete_event_analysis();

    std::cout << "\n=====================================================" << std::endl;
    std::cout << "Summary of Improvements:" << std::endl;
    std::cout << "  ✓ 70% reduction in boilerplate code" << std::endl;
    std::cout << "  ✓ Eliminated hardcoded masses" << std::endl;
    std::cout << "  ✓ Eliminated manual memory management" << std::endl;
    std::cout << "  ✓ Eliminated repetitive coordinate conversions" << std::endl;
    std::cout << "  ✓ Eliminated repetitive boost calls" << std::endl;
    std::cout << "  ✓ Type-safe particle handling" << std::endl;
    std::cout << "  ✓ Self-documenting code" << std::endl;
    std::cout << "=====================================================" << std::endl;

    return 0;
}
