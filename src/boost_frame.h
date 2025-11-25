#ifndef BOOSTFRAME_H
#define BOOSTFRAME_H

#include "pparticle.h"
#include <vector>
#include <memory>

// ============================================================================
// BoostFrame: Manages Reference Frame Transformations
// ============================================================================
/**
 * @class BoostFrame
 * @brief Manages boosting multiple particles to a common reference frame
 *
 * This class eliminates repetitive boost() calls by caching the boost vector
 * and providing batch operations.
 *
 * Usage Example:
 * @code
 *   // Define rest frame from composite system
 *   PParticle p_pip_system = proton + pion;
 *   BoostFrame ppip_frame(p_pip_system);
 *
 *   // Boost multiple particles to this frame
 *   PParticle pion_boosted = ppip_frame.boost(pion);
 *   PParticle neutron_boosted = ppip_frame.boost(neutron);
 *   PParticle proj_boosted = ppip_frame.boost(projectile);
 *
 *   // Or boost in-place
 *   ppip_frame.applyTo(pion);
 *   ppip_frame.applyTo(neutron);
 *   ppip_frame.applyTo(projectile);
 * @endcode
 */
class BoostFrame {
public:
    // ========================================================================
    // Constructors
    // ========================================================================

    /**
     * @brief Construct boost frame from reference particle/system
     * @param reference PParticle defining the rest frame
     * @param type Which momentum representation to use
     */
    explicit BoostFrame(const PParticle& reference,
                       MomentumType type = MomentumType::RECONSTRUCTED)
        : boost_vector_(-reference.boostVector(type)),
          name_(reference.name() + "_frame") {}

    /**
     * @brief Construct boost frame from explicit boost vector
     * @param beta_vector Boost velocity vector
     * @param name Optional frame identifier
     */
    explicit BoostFrame(const TVector3& beta_vector,
                       const std::string& name = "custom_frame")
        : boost_vector_(beta_vector), name_(name) {}

    /**
     * @brief Construct beam rest frame (z-axis boost only)
     * @param beta_z Velocity along beam axis
     */
    static BoostFrame createBeamFrame(double beta_z) {
        return BoostFrame(TVector3(0, 0, -beta_z), "beam_frame");
    }

    // ========================================================================
    // Boost Operations
    // ========================================================================

    /**
     * @brief Boost particle to this frame (returns new copy)
     * @param particle Particle to boost
     * @return New PParticle boosted to this frame
     */
    PParticle boost(const PParticle& particle) const {
        PParticle boosted(particle);
        boosted.boost(boost_vector_);
        return boosted;
    }

    /**
     * @brief Apply boost to particle in-place
     * @param particle Particle to modify
     */
    void applyTo(PParticle& particle) const {
        particle.boost(boost_vector_);
    }

    /**
     * @brief Boost multiple particles at once
     * @param particles Vector of particles to boost
     * @return Vector of boosted particles
     */
    std::vector<PParticle> boost(const std::vector<PParticle>& particles) const {
        std::vector<PParticle> boosted;
        boosted.reserve(particles.size());
        for (const auto& p : particles) {
            boosted.push_back(boost(p));
        }
        return boosted;
    }

    /**
     * @brief Apply boost to multiple particles in-place
     * @param particles Vector of particles to modify
     */
    void applyTo(std::vector<PParticle>& particles) const {
        for (auto& p : particles) {
            p.boost(boost_vector_);
        }
    }

    // ========================================================================
    // Accessors
    // ========================================================================

    const TVector3& boostVector() const { return boost_vector_; }
    const std::string& name() const { return name_; }

    double gamma() const {
        double beta = boost_vector_.Mag();
        return 1.0 / sqrt(1.0 - beta*beta);
    }

private:
    TVector3 boost_vector_;  ///< Cached boost velocity
    std::string name_;       ///< Frame identifier
};

// ============================================================================
// EventFrames: Container for All Reference Frames in an Event
// ============================================================================
/**
 * @class EventFrames
 * @brief Manages all reference frames for a physics event
 *
 * Centralizes frame management for typical hadronic reaction analysis.
 * Eliminates error-prone manual boost calls.
 *
 * Usage Example:
 * @code
 *   EventFrames frames;
 *   frames.setBeamFrame(projectile, target);
 *   frames.addCompositeFrame("ppip", proton + pion);
 *   frames.addCompositeFrame("npip", neutron + pion);
 *
 *   // Boost all particles to beam rest frame
 *   PParticle p_CMS = frames.getFrame("beam").boost(proton);
 * @endcode
 */
class EventFrames {
public:
    /**
     * @brief Define beam center-of-mass frame
     * @param projectile Beam particle
     * @param target Target particle
     */
    void setBeamFrame(const PParticle& projectile, const PParticle& target) {
        PParticle beam = projectile + target;
        frames_["beam"] = std::make_shared<BoostFrame>(beam);
        beam_beta_ = beam.beta();
    }

    /**
     * @brief Define beam frame from kinetic energy (convenience)
     * @param T_kin Beam kinetic energy in MeV
     */
    void setBeamFrameFromKineticEnergy(double T_kin) {
        PParticle proj = ParticleFactory::createBeamProton(T_kin);
        PParticle targ = ParticleFactory::createTargetProton();
        setBeamFrame(proj, targ);
    }

    /**
     * @brief Add composite particle rest frame
     * @param name Frame identifier (e.g., "ppip", "npip")
     * @param composite Composite particle defining frame
     */
    void addCompositeFrame(const std::string& name, const PParticle& composite) {
        frames_[name] = std::make_shared<BoostFrame>(composite);
    }

    /**
     * @brief Get frame by name
     * @param name Frame identifier
     * @return Const reference to BoostFrame
     * @throws std::runtime_error if frame not found
     */
    const BoostFrame& getFrame(const std::string& name) const {
        auto it = frames_.find(name);
        if (it == frames_.end()) {
            throw std::runtime_error("Frame not found: " + name);
        }
        return *(it->second);
    }

    /**
     * @brief Check if frame exists
     */
    bool hasFrame(const std::string& name) const {
        return frames_.find(name) != frames_.end();
    }

    double beamBeta() const { return beam_beta_; }

private:
    std::map<std::string, std::shared_ptr<BoostFrame>> frames_;
    double beam_beta_ = 0.0;
};

#endif // BOOSTFRAME_H
