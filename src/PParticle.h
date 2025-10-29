#ifndef PPARTICLE_H
#define PPARTICLE_H

#include "TLorentzVector.h"
#include "TVector3.h"
#include <memory>
#include <string>
#include <stdexcept>

// ============================================================================
// Particle Physics Constants (PDG 2024 values)
// ============================================================================
namespace Physics {
    constexpr double MASS_PROTON   = 938.27231;   // MeV/c^2
    constexpr double MASS_NEUTRON  = 939.56542;   // MeV/c^2
    constexpr double MASS_PION_PLUS  = 139.56995; // MeV/c^2
    constexpr double MASS_PION_MINUS = 139.56995; // MeV/c^2
    constexpr double MASS_PION_ZERO  = 134.9768;  // MeV/c^2

    constexpr double D2R = 1.74532925199432955e-02; // Degrees to radians
    constexpr double R2D = 57.2957795130823229;     // Radians to degrees
}

// ============================================================================
// Momentum Type Enumeration
// ============================================================================
enum class MomentumType {
    RECONSTRUCTED,  // Default: measured momentum from detector
    CORRECTED,      // Energy loss corrected
    SIMULATED       // MC truth (if available)
};

// ============================================================================
// PParticle: Generic Particle Class with Multiple Momentum Representations
// ============================================================================
/**
 * @class PParticle
 * @brief Encapsulates a particle with reconstructed, corrected, and simulated momenta
 *
 * Key Features:
 * - Stores up to 3 momentum representations per particle
 * - Transparent boosting to any reference frame
 * - Composite particle creation via operator+
 * - Automatic mass assignment
 * - LAB frame preservation
 *
 * Usage Example:
 * @code
 *   // Create proton from spherical coordinates
 *   PParticle proton(Physics::MASS_PROTON);
 *   proton.setFromSpherical(1580, 45.0, 30.0, MomentumType::RECONSTRUCTED);
 *   proton.setFromSpherical(1590, 45.0, 30.0, MomentumType::CORRECTED);
 *
 *   // Boost to beam rest frame
 *   TVector3 beam_beta(0, 0, 0.85);
 *   proton.boost(-beam_beta);
 *
 *   // Create composite
 *   PParticle delta_pp = proton + pion;
 * @endcode
 */
class PParticle {
public:
    // ========================================================================
    // Constructors
    // ========================================================================

    /**
     * @brief Construct particle with given rest mass
     * @param mass Rest mass in MeV/c^2
     * @param name Optional particle name for debugging
     */
    explicit PParticle(double mass, const std::string& name = "")
        : mass_(mass), name_(name) {}

    /**
     * @brief Construct from existing TLorentzVector
     * @param p4 Four-momentum vector
     * @param name Optional particle name
     */
    PParticle(const TLorentzVector& p4, const std::string& name = "")
        : mass_(p4.M()), name_(name) {
        p4_reconstructed_ = p4;
        lab_frame_reconstructed_ = p4;
    }

    /**
     * @brief Copy constructor with deep copy of all momentum representations
     */
    PParticle(const PParticle& other)
        : mass_(other.mass_), name_(other.name_),
          p4_reconstructed_(other.p4_reconstructed_),
          p4_corrected_(other.p4_corrected_),
          p4_simulated_(other.p4_simulated_),
          lab_frame_reconstructed_(other.lab_frame_reconstructed_),
          lab_frame_corrected_(other.lab_frame_corrected_),
          lab_frame_simulated_(other.lab_frame_simulated_) {}

    // ========================================================================
    // Momentum Setters: Spherical Coordinates (p, theta, phi)
    // ========================================================================

    /**
     * @brief Set momentum from spherical coordinates
     * @param p Momentum magnitude in MeV/c
     * @param theta_deg Polar angle in degrees [0, 180]
     * @param phi_deg Azimuthal angle in degrees [0, 360]
     * @param type Which momentum representation to set
     */
    void setFromSpherical(double p, double theta_deg, double phi_deg,
                         MomentumType type = MomentumType::RECONSTRUCTED) {
        TVector3 p3 = sphericalToCartesian(p, theta_deg, phi_deg);
        setFromVector(p3, type);
    }

    /**
     * @brief Set momentum from 3-vector
     * @param p3 Three-momentum vector
     * @param type Which momentum representation to set
     */
    void setFromVector(const TVector3& p3,
                      MomentumType type = MomentumType::RECONSTRUCTED) {
        TLorentzVector p4;
        p4.SetVectM(p3, mass_);

        switch (type) {
            case MomentumType::RECONSTRUCTED:
                p4_reconstructed_ = p4;
                lab_frame_reconstructed_ = p4;
                break;
            case MomentumType::CORRECTED:
                p4_corrected_ = p4;
                lab_frame_corrected_ = p4;
                break;
            case MomentumType::SIMULATED:
                p4_simulated_ = p4;
                lab_frame_simulated_ = p4;
                break;
        }
    }

    /**
     * @brief Set momentum from Cartesian coordinates
     * @param px, py, pz Momentum components in MeV/c
     * @param type Which momentum representation to set
     */
    void setFromCartesian(double px, double py, double pz,
                         MomentumType type = MomentumType::RECONSTRUCTED) {
        TVector3 p3(px, py, pz);
        setFromVector(p3, type);
    }

    // ========================================================================
    // Momentum Getters
    // ========================================================================

    /**
     * @brief Get four-momentum for specified type
     * @param type Momentum representation to retrieve
     * @return Const reference to TLorentzVector
     * @throws std::runtime_error if requested type not set
     */
    const TLorentzVector& vec(MomentumType type = MomentumType::RECONSTRUCTED) const {
        switch (type) {
            case MomentumType::RECONSTRUCTED:
                return p4_reconstructed_;
            case MomentumType::CORRECTED:
                if (p4_corrected_.E() == 0)
                    throw std::runtime_error("Corrected momentum not set for " + name_);
                return p4_corrected_;
            case MomentumType::SIMULATED:
                if (p4_simulated_.E() == 0)
                    throw std::runtime_error("Simulated momentum not set for " + name_);
                return p4_simulated_;
        }
        return p4_reconstructed_; // Default fallback
    }

    /**
     * @brief Get mutable reference to four-momentum
     * @warning Use sparingly; prefer immutable interface
     */
    TLorentzVector& vecMutable(MomentumType type = MomentumType::RECONSTRUCTED) {
        return const_cast<TLorentzVector&>(vec(type));
    }

    /**
     * @brief Get LAB frame momentum (before any boosts)
     * @param type Momentum representation to retrieve
     * @return Const reference to LAB frame TLorentzVector
     */
    const TLorentzVector& labFrame(MomentumType type = MomentumType::RECONSTRUCTED) const {
        switch (type) {
            case MomentumType::RECONSTRUCTED:
                return lab_frame_reconstructed_;
            case MomentumType::CORRECTED:
                return lab_frame_corrected_;
            case MomentumType::SIMULATED:
                return lab_frame_simulated_;
        }
        return lab_frame_reconstructed_;
    }

    // ========================================================================
    // Reference Frame Transformations
    // ========================================================================

    /**
     * @brief Apply Lorentz boost to all momentum representations
     * @param beta_vector Boost vector (beta_x, beta_y, beta_z)
     *
     * This boosts ALL momentum types that have been set.
     * LAB frame copies remain unchanged for later reference.
     */
    void boost(const TVector3& beta_vector) {
        if (p4_reconstructed_.E() != 0)
            p4_reconstructed_.Boost(beta_vector);
        if (p4_corrected_.E() != 0)
            p4_corrected_.Boost(beta_vector);
        if (p4_simulated_.E() != 0)
            p4_simulated_.Boost(beta_vector);
    }

    /**
     * @brief Apply z-axis boost (for collinear beam geometry)
     * @param beta_z Boost velocity along z-axis [-1, 1]
     */
    void boostZ(double beta_z) {
        boost(TVector3(0, 0, beta_z));
    }

    /**
     * @brief Boost to rest frame of another particle/system
     * @param reference_system PParticle or composite defining rest frame
     * @param type Which momentum representation to use for boost
     */
    void boostToRestFrame(const PParticle& reference_system,
                         MomentumType type = MomentumType::RECONSTRUCTED) {
        TVector3 boost_vec = -reference_system.vec(type).BoostVector();
        boost(boost_vec);
    }

    /**
     * @brief Reset to LAB frame (undo all boosts)
     */
    void resetToLAB() {
        p4_reconstructed_ = lab_frame_reconstructed_;
        p4_corrected_ = lab_frame_corrected_;
        p4_simulated_ = lab_frame_simulated_;
    }

    // ========================================================================
    // Composite Particle Creation
    // ========================================================================

    /**
     * @brief Create composite particle by adding four-momenta
     * @param other Another particle to combine with
     * @return New PParticle representing composite system
     *
     * Combines all available momentum representations.
     * Composite mass is computed from invariant mass.
     */
    PParticle operator+(const PParticle& other) const {
        TLorentzVector composite_p4 = p4_reconstructed_ + other.p4_reconstructed_;
        PParticle composite(composite_p4.M(), name_ + "+" + other.name_);

        // Reconstructed
        composite.p4_reconstructed_ = composite_p4;
        composite.lab_frame_reconstructed_ = composite_p4;

        // Corrected (if both available)
        if (p4_corrected_.E() != 0 && other.p4_corrected_.E() != 0) {
            composite.p4_corrected_ = p4_corrected_ + other.p4_corrected_;
            composite.lab_frame_corrected_ = composite.p4_corrected_;
        }

        // Simulated (if both available)
        if (p4_simulated_.E() != 0 && other.p4_simulated_.E() != 0) {
            composite.p4_simulated_ = p4_simulated_ + other.p4_simulated_;
            composite.lab_frame_simulated_ = composite.p4_simulated_;
        }

        return composite;
    }

    /**
     * @brief Subtract four-momentum (for missing mass calculations)
     */
    PParticle operator-(const PParticle& other) const {
        TLorentzVector diff_p4 = p4_reconstructed_ - other.p4_reconstructed_;
        PParticle result(diff_p4.M(), name_ + "-" + other.name_);
        result.p4_reconstructed_ = diff_p4;
        result.lab_frame_reconstructed_ = diff_p4;
        return result;
    }

    // ========================================================================
    // Kinematic Accessors (shortcuts to most common quantities)
    // ========================================================================

    double mass(MomentumType type = MomentumType::RECONSTRUCTED) const {
        return vec(type).M();
    }

    double massGeV(MomentumType type = MomentumType::RECONSTRUCTED) const {
        return vec(type).M() / 1000.0;
    }

    double momentum(MomentumType type = MomentumType::RECONSTRUCTED) const {
        return vec(type).P();
    }

    double energy(MomentumType type = MomentumType::RECONSTRUCTED) const {
        return vec(type).E();
    }

    double theta(MomentumType type = MomentumType::RECONSTRUCTED) const {
        return vec(type).Theta() * Physics::R2D;
    }

    double phi(MomentumType type = MomentumType::RECONSTRUCTED) const {
        return vec(type).Phi() * Physics::R2D;
    }

    double cosTheta(MomentumType type = MomentumType::RECONSTRUCTED) const {
        return vec(type).CosTheta();
    }

    double rapidity(MomentumType type = MomentumType::RECONSTRUCTED) const {
        return vec(type).Rapidity();
    }

    double beta(MomentumType type = MomentumType::RECONSTRUCTED) const {
        return vec(type).Beta();
    }

    TVector3 boostVector(MomentumType type = MomentumType::RECONSTRUCTED) const {
        return vec(type).BoostVector();
    }

    // ========================================================================
    // Utility Methods
    // ========================================================================

    /**
     * @brief Calculate opening angle with another particle
     * @param other Another particle
     * @param type Momentum representation to use
     * @return Opening angle in degrees
     */
    double openingAngle(const PParticle& other,
                       MomentumType type = MomentumType::RECONSTRUCTED) const {
        return vec(type).Angle(other.vec(type).Vect()) * Physics::R2D;
    }

    /**
     * @brief Calculate relative angle in azimuthal plane
     * @param other Another particle
     * @param type Momentum representation to use
     * @return Delta phi in degrees
     */
    double deltaPhi(const PParticle& other,
                   MomentumType type = MomentumType::RECONSTRUCTED) const {
        double dphi = vec(type).Phi() - other.vec(type).Phi();
        // Wrap to [-pi, pi]
        while (dphi > TMath::Pi()) dphi -= 2*TMath::Pi();
        while (dphi < -TMath::Pi()) dphi += 2*TMath::Pi();
        return dphi * Physics::R2D;
    }

    const std::string& name() const { return name_; }
    void setName(const std::string& name) { name_ = name; }

    /**
     * @brief Print particle information (for debugging)
     */
    void print(MomentumType type = MomentumType::RECONSTRUCTED) const {
        const TLorentzVector& p4 = vec(type);
        std::cout << "PParticle: " << name_ << std::endl;
        std::cout << "  Mass: " << mass_ << " MeV/c^2" << std::endl;
        std::cout << "  (E, px, py, pz) = ("
                  << p4.E() << ", " << p4.Px() << ", "
                  << p4.Py() << ", " << p4.Pz() << ")" << std::endl;
        std::cout << "  (p, theta, phi) = ("
                  << momentum(type) << ", " << theta(type) << ", "
                  << phi(type) << ")" << std::endl;
    }

private:
    // ========================================================================
    // Internal Helper Methods
    // ========================================================================

    /**
     * @brief Convert spherical to Cartesian coordinates
     */
    static TVector3 sphericalToCartesian(double p, double theta_deg, double phi_deg) {
        double theta_rad = theta_deg * Physics::D2R;
        double phi_rad = phi_deg * Physics::D2R;
        return TVector3(
            p * sin(theta_rad) * cos(phi_rad),
            p * sin(theta_rad) * sin(phi_rad),
            p * cos(theta_rad)
        );
    }

    // ========================================================================
    // Data Members
    // ========================================================================

    double mass_;                           ///< Rest mass in MeV/c^2
    std::string name_;                      ///< Particle identifier

    // Current four-momenta (subject to boosts)
    TLorentzVector p4_reconstructed_;       ///< Measured momentum
    TLorentzVector p4_corrected_;           ///< Energy loss corrected
    TLorentzVector p4_simulated_;           ///< MC truth (if available)

    // LAB frame copies (immune to boosts)
    TLorentzVector lab_frame_reconstructed_;
    TLorentzVector lab_frame_corrected_;
    TLorentzVector lab_frame_simulated_;
};

// ============================================================================
// Convenience Factory Functions
// ============================================================================

namespace ParticleFactory {
    /**
     * @brief Create proton from TNtuple variables
     * @param p_corr Corrected momentum magnitude
     * @param theta Polar angle in degrees
     * @param phi Azimuthal angle in degrees
     * @return PParticle configured as proton
     */
    inline PParticle createProton(double p_corr, double theta, double phi) {
        PParticle proton(Physics::MASS_PROTON, "p");
        proton.setFromSpherical(p_corr, theta, phi, MomentumType::CORRECTED);
        return proton;
    }

    /**
     * @brief Create positive pion from TNtuple variables
     */
    inline PParticle createPionPlus(double p_corr, double theta, double phi) {
        PParticle pion(Physics::MASS_PION_PLUS, "pi+");
        pion.setFromSpherical(p_corr, theta, phi, MomentumType::CORRECTED);
        return pion;
    }

    /**
     * @brief Create beam proton from kinetic energy
     * @param T_kin Kinetic energy in MeV
     * @return PParticle representing beam along +z
     */
    inline PParticle createBeamProton(double T_kin) {
        double E = T_kin + Physics::MASS_PROTON;
        double p = sqrt(E*E - Physics::MASS_PROTON*Physics::MASS_PROTON);
        PParticle beam(Physics::MASS_PROTON, "beam");
        beam.setFromCartesian(0, 0, p);
        return beam;
    }

    /**
     * @brief Create target proton at rest
     */
    inline PParticle createTargetProton() {
        PParticle target(Physics::MASS_PROTON, "target");
        target.setFromCartesian(0, 0, 0);
        return target;
    }
}

#endif // PPARTICLE_H
