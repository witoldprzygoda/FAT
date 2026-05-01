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
// Dihedral / plane angle between two pairs of particles
// ============================================================================
/**
 * @brief Angle between the planes spanned by two particle pairs.
 *
 * For pairs (a, b) and (c, d) — typically AFTER boosting all four particles
 * to a common rest frame — compute:
 *   n_ab = p_a × p_b   (normal to plane (a, b))
 *   n_cd = p_c × p_d   (normal to plane (c, d))
 *   φ    = angle(n_ab, n_cd)   ∈ [0, 180] deg
 *
 * Used for the pippimepem analysis as the angle between the (π+π-) and
 * (e+e-) decay planes in the pippimepem rest frame — a meaningful kinematic
 * observable for the 4-body decay topology. NOT the same as the (vacuous)
 * opening angle between the two pair-momenta in the same frame, which is
 * always 180° by 3-momentum conservation.
 *
 * @return plane-to-plane angle in degrees, range [0, 180]
 */
inline double planeAngle(const PParticle& a, const PParticle& b,
                         const PParticle& c, const PParticle& d,
                         KinematicType type = KinematicType::RECONSTRUCTED) {
    TVector3 n_ab = a.vec(type).Vect().Cross(b.vec(type).Vect());
    TVector3 n_cd = c.vec(type).Vect().Cross(d.vec(type).Vect());
    return n_ab.Angle(n_cd) * TMath::RadToDeg();
}

}  // namespace Physics

#endif // PHYSICS_UTILS_H
