/**
 * @file pparticle_fwd.h
 * @brief PParticleFwd - Extended particle class for Forward Tracker (FT) data
 *
 * PParticleFwd inherits from PParticle and adds Forward Tracker specific fields.
 * It can be used in arithmetic operations with PParticle objects (e.g., fw_p1 + pion).
 *
 * Forward Tracker variables in ntuple:
 * - fwdet_beta_N, fwdet_p_N, fwdet_mass_N, fwdet_mass2_N, fwdet_q_N
 * - fwdet_tofrec_N, fwdet_phi_N, fwdet_theta_N, fwdet_r_N, fwdet_z_N
 * - fwdet_chi2_N, fwdet_phi2_N, fwdet_theta2_N, fwdet_r2_N, fwdet_z2_N
 * - fwdet_dist_tofhit_N, fwdet_tof_N, fwdet_ndf_N
 * - fwdet_mult, fwdet_counter
 *
 * Where N = 1, 2, 3 for the three FT hits.
 *
 * @author Witold Przygoda (witold.przygoda@uj.edu.pl)
 * @date 2025
 */

#ifndef PPARTICLE_FWD_H
#define PPARTICLE_FWD_H

#include "pparticle.h"
#include "ntuple_reader.h"
#include <string>

/**
 * @class PParticleFwd
 * @brief Extended PParticle for Forward Tracker data
 *
 * Inherits all PParticle functionality (momentum representations, boosts,
 * arithmetic operations) and adds FT-specific fields.
 *
 * Usage:
 * @code
 *   PParticleFwd fw_p1(Physics::MASS_PROTON, "fw_p1");
 *   fw_p1.setFromReader(reader, 1);  // Fill from fwdet_*_1 variables
 *   
 *   // Can combine with other particles
 *   PParticle composite = fw_p1 + pion;
 * @endcode
 */
class PParticleFwd : public PParticle {
public:
    // ========================================================================
    // Forward Tracker specific fields
    // ========================================================================
    
    // Measured quantities
    double fwd_beta;          ///< Measured velocity (beta = v/c)
    double fwd_mass;          ///< Reconstructed mass [MeV/c^2]
    double fwd_mass2;         ///< Mass squared [MeV^2/c^4]
    double fwd_charge;        ///< Charge
    double fwd_tofRec;        ///< Reconstructed time of flight
    double fwd_tof;           ///< Time of flight
    
    // Position/angles - primary reconstruction
    double fwd_phi;           ///< Azimuthal angle [deg]
    double fwd_theta;         ///< Polar angle [deg]
    double fwd_r;             ///< Radial position [mm]
    double fwd_z;             ///< Z position [mm]
    
    // Position/angles - secondary reconstruction
    double fwd_phi2;          ///< Azimuthal angle (2nd method) [deg]
    double fwd_theta2;        ///< Polar angle (2nd method) [deg]
    double fwd_r2;            ///< Radial position (2nd method) [mm]
    double fwd_z2;            ///< Z position (2nd method) [mm]
    
    // Quality parameters
    double fwd_chi2;          ///< Fit chi-squared
    int    fwd_ndf;           ///< Number of degrees of freedom
    double fwd_chi2ndf;       ///< chi2/ndf (computed)
    double fwd_distToRpc;     ///< Distance to RPC hit [mm]
    
    // Multiplicity info (same for all FT particles in event)
    int    fwd_mult;          ///< Total FT multiplicity in event
    int    fwd_counter;       ///< Counter (hit index)
    
    // Validity flag
    bool   fwd_valid;         ///< Whether this particle was successfully filled

    // ========================================================================
    // Constructors
    // ========================================================================
    
    /**
     * @brief Construct empty FT particle with given mass hypothesis
     * @param mass Mass hypothesis in MeV/c^2 (default: proton)
     * @param name Particle name (e.g., "fw_p1")
     */
    explicit PParticleFwd(double mass = Physics::MASS_PROTON, 
                          const std::string& name = "fw_p")
        : PParticle(mass, name), fwd_valid(false)
    {
        resetFwdFields();
    }
    
    /**
     * @brief Copy constructor
     */
    PParticleFwd(const PParticleFwd& other)
        : PParticle(other),
          fwd_beta(other.fwd_beta), fwd_mass(other.fwd_mass), 
          fwd_mass2(other.fwd_mass2), fwd_charge(other.fwd_charge),
          fwd_tofRec(other.fwd_tofRec), fwd_tof(other.fwd_tof),
          fwd_phi(other.fwd_phi), fwd_theta(other.fwd_theta),
          fwd_r(other.fwd_r), fwd_z(other.fwd_z),
          fwd_phi2(other.fwd_phi2), fwd_theta2(other.fwd_theta2),
          fwd_r2(other.fwd_r2), fwd_z2(other.fwd_z2),
          fwd_chi2(other.fwd_chi2), fwd_ndf(other.fwd_ndf),
          fwd_chi2ndf(other.fwd_chi2ndf), fwd_distToRpc(other.fwd_distToRpc),
          fwd_mult(other.fwd_mult), fwd_counter(other.fwd_counter),
          fwd_valid(other.fwd_valid) {}

    // ========================================================================
    // Data Loading from NTupleReader
    // ========================================================================
    
    /**
     * @brief Fill particle from ntuple reader
     * @param reader NTupleReader with loaded event
     * @param index FT hit index (1, 2, or 3)
     * @param type Kinematic type to set (default: RECONSTRUCTED, for future: CORRECTED, SIMULATED)
     * @return true if successfully filled, false if variables missing
     *
     * Reads all fwdet_*_N variables where N = index.
     * Sets momentum using fwdet_p_N, fwdet_theta_N, fwdet_phi_N.
     */
    bool setFromReader(NTupleReader& reader, int index,
                       KinematicType type = KinematicType::RECONSTRUCTED) {
        
        // Build variable name suffix
        std::string idx = "_" + std::to_string(index);
        
        // Check if primary variables exist
        std::string p_var = "fwdet_p" + idx;
        std::string theta_var = "fwdet_theta" + idx;
        std::string phi_var = "fwdet_phi" + idx;
        
        if (!reader.hasVariable(p_var) || 
            !reader.hasVariable(theta_var) || 
            !reader.hasVariable(phi_var)) {
            fwd_valid = false;
            return false;
        }
        
        // Read kinematic variables
        double p_val = reader[p_var];
        double theta_val = reader[theta_var];
        double phi_val = reader[phi_var];
        
        // Set momentum (inherited from PParticle)
        setFromSpherical(p_val, theta_val, phi_val, type);
        
        // Read FT-specific fields
        fwd_phi = phi_val;
        fwd_theta = theta_val;
        
        fwd_beta      = readIfExists(reader, "fwdet_beta" + idx, -1.0);
        fwd_mass      = readIfExists(reader, "fwdet_mass" + idx, -1.0);
        fwd_mass2     = readIfExists(reader, "fwdet_mass2" + idx, -1.0);
        fwd_charge    = readIfExists(reader, "fwdet_q" + idx, 0.0);
        fwd_tofRec    = readIfExists(reader, "fwdet_tofrec" + idx, -1.0);
        fwd_tof       = readIfExists(reader, "fwdet_tof" + idx, -1.0);
        
        fwd_r         = readIfExists(reader, "fwdet_r" + idx, -1.0);
        fwd_z         = readIfExists(reader, "fwdet_z" + idx, -1.0);
        fwd_phi2      = readIfExists(reader, "fwdet_phi2" + idx, -1.0);
        fwd_theta2    = readIfExists(reader, "fwdet_theta2" + idx, -1.0);
        fwd_r2        = readIfExists(reader, "fwdet_r2" + idx, -1.0);
        fwd_z2        = readIfExists(reader, "fwdet_z2" + idx, -1.0);
        
        fwd_chi2      = readIfExists(reader, "fwdet_chi2" + idx, -1.0);
        fwd_ndf       = static_cast<int>(readIfExists(reader, "fwdet_ndf" + idx, 0.0));
        fwd_chi2ndf   = (fwd_ndf > 0) ? fwd_chi2 / fwd_ndf : 1e9;
        fwd_distToRpc = readIfExists(reader, "fwdet_dist_tofhit" + idx, -1.0);
        
        // Multiplicity (same for all hits, no index suffix)
        fwd_mult      = static_cast<int>(readIfExists(reader, "fwdet_mult", 0.0));
        fwd_counter   = static_cast<int>(readIfExists(reader, "fwdet_counter", 0.0));
        
        fwd_valid = true;
        return true;
    }
    
    /**
     * @brief Reset FT-specific fields to default values
     */
    void resetFwdFields() {
        fwd_beta = -1.0;
        fwd_mass = -1.0;
        fwd_mass2 = -1.0;
        fwd_charge = 0.0;
        fwd_tofRec = -1.0;
        fwd_tof = -1.0;
        fwd_phi = -1.0;
        fwd_theta = -1.0;
        fwd_r = -1.0;
        fwd_z = -1.0;
        fwd_phi2 = -1.0;
        fwd_theta2 = -1.0;
        fwd_r2 = -1.0;
        fwd_z2 = -1.0;
        fwd_chi2 = -1.0;
        fwd_ndf = 0;
        fwd_chi2ndf = 1e9;
        fwd_distToRpc = -1.0;
        fwd_mult = 0;
        fwd_counter = 0;
        fwd_valid = false;
    }
    
    // ========================================================================
    // Quality Accessors
    // ========================================================================
    
    /**
     * @brief Check if particle is valid (was successfully filled)
     */
    bool isValid() const { return fwd_valid; }
    
    /**
     * @brief Get chi2/ndf quality
     */
    double getChi2NDF() const { return fwd_chi2ndf; }
    
    /**
     * @brief Get distance to RPC hit
     */
    double getDistToRpc() const { return fwd_distToRpc; }
    
    /**
     * @brief Get measured beta
     */
    double getMeasuredBeta() const { return fwd_beta; }
    
    /**
     * @brief Check if passes basic quality cuts
     * @param maxChi2ndf Maximum chi2/ndf (default: 10)
     * @param maxDistRpc Maximum distance to RPC (default: no cut = 1e9)
     */
    bool passesQuality(double maxChi2ndf = 10.0, double maxDistRpc = 1e9) const {
        if (!fwd_valid) return false;
        if (fwd_chi2ndf > maxChi2ndf) return false;
        if (fwd_distToRpc > maxDistRpc && fwd_distToRpc > 0) return false;
        return true;
    }
    
    /**
     * @brief Print FT-specific info
     */
    void printFwd() const {
        print();  // Base class print
        std::cout << "  FT Info:" << std::endl;
        std::cout << "    beta: " << fwd_beta << ", mass: " << fwd_mass << " MeV/c^2" << std::endl;
        std::cout << "    chi2/ndf: " << fwd_chi2ndf << " (" << fwd_chi2 << "/" << fwd_ndf << ")" << std::endl;
        std::cout << "    r: " << fwd_r << " mm, z: " << fwd_z << " mm" << std::endl;
        std::cout << "    distToRpc: " << fwd_distToRpc << " mm" << std::endl;
    }

private:
    /**
     * @brief Helper to read variable if exists, else return default
     */
    double readIfExists(NTupleReader& reader, const std::string& name, double def) {
        return reader.hasVariable(name) ? reader[name] : def;
    }
};

#endif // PPARTICLE_FWD_H
