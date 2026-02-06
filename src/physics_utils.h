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
    return p1.vec(type).Angle(p2.vec(type)) * TMath::RadToDeg();
}

}  // namespace Physics

#endif // PHYSICS_UTILS_H
