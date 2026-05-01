/**
 * @file physics_utils.h
 * @brief Physics utility functions for particle analysis
 *
 * This header provides common physics calculations that operate on
 * particle pairs or groups. Functions are symmetric (treat both
 * particles equally) and work with any PParticle-like objects.
 *
 * @author Witold Przygoda (witold.przygoda@uj.edu.pl)
 * @date 2025
 */

#ifndef PHYSICS_UTILS_H
#define PHYSICS_UTILS_H

#include "pparticle.h"
#include <TMath.h>

namespace Physics {

// ============================================================================
// Two-particle functions
// ============================================================================

/**
 * @brief Calculate opening angle between two particles
 * @param p1 First particle
 * @param p2 Second particle
 * @param type Kinematic type (RECONSTRUCTED, CORRECTED, SIMULATED)
 * @return Opening angle in degrees
 *
 * The opening angle is the angle between the momentum vectors of two particles.
 * Must be calculated BEFORE combining particles with operator+.
 *
 * Example:
 * @code
 *   double oa = Physics::openingAngle(positron, electron);
 *   if (oa < 9.0) continue;  // reject close pairs
 *   PParticle dilepton = positron + electron;
 * @endcode
 */
inline double openingAngle(const PParticle& p1, const PParticle& p2,
                           KinematicType type = KinematicType::RECONSTRUCTED) {
    return p1.vec(type).Vect().Angle(p2.vec(type).Vect()) * TMath::RadToDeg();
}

// ============================================================================
// Opening angle in a mass-constrained (hypothetical) rest frame
// ============================================================================
/**
 * @brief Opening angle between two compound 4-momenta evaluated in the rest
 *        frame of a HYPOTHETICAL parent whose 3-momentum equals (a+b) but
 *        whose mass is forced to @p m_constraint.
 *
 * Construction:
 *   P_total = a + b                                    (measured 4-momentum)
 *   P_hyp   = (E_hyp, vec(P_total))   with E_hyp = √(|P_total|² + m_constraint²)
 *   boost a, b into the rest frame of P_hyp
 *   return ∠(vec(a*), vec(b*)) in degrees
 *
 * Why not the true (a+b) rest frame? In the true frame the two pair-momenta
 * are exactly back-to-back by construction (always 180°), which is vacuous.
 * Forcing the parent mass to a fixed hypothesis (e.g. m_η for the pippimepem
 * 4-body) breaks that constraint: when M(a+b) ≠ m_constraint the boost differs
 * from the true one and the resulting OA is a meaningful kinematic observable
 * that peaks near 180° for genuine η decays and spreads out for backgrounds.
 *
 * @param a              first compound (e.g. pippim)
 * @param b              second compound (e.g. epem)
 * @param m_constraint   hypothetical parent mass in MeV/c² (same units as PParticle)
 * @return opening angle in degrees, range [0, 180]
 */
inline double openingAngleInMassConstrainedRestFrame(
        const PParticle& a, const PParticle& b,
        double m_constraint,
        KinematicType type = KinematicType::RECONSTRUCTED) {
    TLorentzVector v1 = a.vec(type);
    TLorentzVector v2 = b.vec(type);
    TVector3 p_tot = (v1 + v2).Vect();
    double   E_hyp = std::sqrt(p_tot.Mag2() + m_constraint * m_constraint);
    TLorentzVector p4_hyp(p_tot, E_hyp);
    TVector3 boost_to_rest = -p4_hyp.BoostVector();
    v1.Boost(boost_to_rest);
    v2.Boost(boost_to_rest);
    return v1.Vect().Angle(v2.Vect()) * TMath::RadToDeg();
}

}  // namespace Physics

#endif // PHYSICS_UTILS_H
