/**
 * @file pparticle_ecal.h
 * @brief PParticleEcal - Extended particle class for Electromagnetic Calorimeter (ECAL) data
 *
 * PParticleEcal inherits from PParticle and adds ECAL/EMC specific fields.
 * It can be used in arithmetic operations with PParticle objects (e.g., ecal_n1 + proton).
 *
 * ECAL neutral candidate variables in ntuple:
 * - neutr_pid_N, neutr_tof_N, neutr_dist_N, neutr_clusterid_N
 * - neutr_beta_N, neutr_p_N, neutr_p_pid_N, neutr_mass_N, neutr_mass2_N
 * - neutr_q_N, neutr_tofrec_N, neutr_phi_N, neutr_theta_N
 * - neutr_r_N, neutr_z_N, neutr_chi2_N
 * - neutr_phi2_N, neutr_theta2_N, neutr_r2_N, neutr_z2_N
 * - neutr_energy_N
 * - neutr_cluster_energy_N, neutr_cluster_theta_N, neutr_cluster_phi_N (cluster reco)
 * - neutr_mult, neutr_counter
 *
 * Where N = 1, 2, 3, 4, 5 for the five ECAL hits.
 *
 * The cluster variables (neutr_cluster_*) are used for PParticle kinematics.
 *
 * Typical use cases:
 * - Photon (gamma): mass = 0
 * - Electron/Positron: mass = MASS_ELECTRON
 * - Neutral pion (pi0): mass = MASS_PION_ZERO
 *
 * @author Witold Przygoda (witold.przygoda@uj.edu.pl)
 * @date 2025
 */

#ifndef PPARTICLE_ECAL_H
#define PPARTICLE_ECAL_H

#include "pparticle.h"
#include "ntuple_reader.h"
#include <string>

/**
 * @class PParticleEcal
 * @brief Extended PParticle for Electromagnetic Calorimeter data
 *
 * Inherits all PParticle functionality (momentum representations, boosts,
 * arithmetic operations) and adds ECAL-specific fields.
 *
 * Usage:
 * @code
 *   // Photon hypothesis (mass = 0)
 *   PParticleEcal ecal_n1(0.0, "ecal_n1");
 *   ecal_n1.setFromReader(reader, 1);  // Fill from neutr_*_1 variables
 *   
 *   // Electron hypothesis
 *   PParticleEcal ecal_e1(Physics::MASS_ELECTRON, "ecal_e1");
 *   ecal_e1.setFromReader(reader, 1);
 *   
 *   // Can combine with other particles
 *   PParticle composite = ecal_n1 + proton;
 * @endcode
 */
class PParticleEcal : public PParticle {
public:
    // ========================================================================
    // ECAL specific fields
    // ========================================================================

    // Identification
    int    ecal_pid;          ///< Particle ID
    int    ecal_clusterId;    ///< EMC cluster index
    double ecal_charge;       ///< Charge

    // Energy/momentum
    double ecal_energy;       ///< Energy from calorimeter [MeV]
    double ecal_p_pid;        ///< Momentum from PID [MeV/c]

    // Mass
    double ecal_mass;         ///< Reconstructed mass [MeV/c^2]
    double ecal_mass2;        ///< Mass squared [MeV^2/c^4]

    // Timing
    double ecal_beta;         ///< Measured velocity (beta = v/c)
    double ecal_tof;          ///< Time of flight
    double ecal_tofRec;       ///< Reconstructed TOF

    // Distance
    double ecal_dist;         ///< Distance to EMC [mm]

    // Position/angles - primary reconstruction
    double ecal_phi;          ///< Azimuthal angle [deg]
    double ecal_theta;        ///< Polar angle [deg]
    double ecal_r;            ///< Radial position [mm]
    double ecal_z;            ///< Z position [mm]

    // Position/angles - secondary reconstruction
    double ecal_phi2;         ///< Azimuthal angle (2nd method) [deg]
    double ecal_theta2;       ///< Polar angle (2nd method) [deg]
    double ecal_r2;           ///< Radial position (2nd method) [mm]
    double ecal_z2;           ///< Z position (2nd method) [mm]

    // ========================================================================
    // ECAL cluster reconstruction (used for PParticle kinematics)
    // ========================================================================
    double cluster_energy;    ///< Cluster energy [MeV] - used for momentum
    double cluster_theta;     ///< Cluster polar angle [deg] - used for kinematics
    double cluster_phi;       ///< Cluster azimuthal angle [deg] - used for kinematics
    int    cluster_ncells;    ///< Number of cells in cluster (cluster size)

    // Quality
    double ecal_chi2;         ///< Fit chi-squared

    // Multiplicity info (same for all ECAL particles in event)
    int    ecal_mult;         ///< Total ECAL multiplicity in event
    int    ecal_counter;      ///< Counter (hit index)

    // Validity flag
    bool   ecal_valid;        ///< Whether this particle was successfully filled

    // ========================================================================
    // Constructors
    // ========================================================================
    
    /**
     * @brief Construct empty ECAL particle with given mass hypothesis
     * @param mass Mass hypothesis in MeV/c^2 (default: 0 for photon)
     * @param name Particle name (e.g., "ecal_n1", "gamma", "e-")
     *
     * Common mass hypotheses:
     * - 0.0: Photon (gamma)
     * - Physics::MASS_ELECTRON: Electron/Positron
     * - Physics::MASS_PION_ZERO: Neutral pion
     */
    explicit PParticleEcal(double mass = 0.0, 
                           const std::string& name = "ecal_n")
        : PParticle(mass, name), ecal_valid(false)
    {
        resetEcalFields();
    }
    
    /**
     * @brief Copy constructor
     */
    PParticleEcal(const PParticleEcal& other)
        : PParticle(other),
          ecal_pid(other.ecal_pid), ecal_clusterId(other.ecal_clusterId),
          ecal_charge(other.ecal_charge),
          ecal_energy(other.ecal_energy), ecal_p_pid(other.ecal_p_pid),
          ecal_mass(other.ecal_mass), ecal_mass2(other.ecal_mass2),
          ecal_beta(other.ecal_beta), ecal_tof(other.ecal_tof),
          ecal_tofRec(other.ecal_tofRec), ecal_dist(other.ecal_dist),
          ecal_phi(other.ecal_phi), ecal_theta(other.ecal_theta),
          ecal_r(other.ecal_r), ecal_z(other.ecal_z),
          ecal_phi2(other.ecal_phi2), ecal_theta2(other.ecal_theta2),
          ecal_r2(other.ecal_r2), ecal_z2(other.ecal_z2),
          cluster_energy(other.cluster_energy),
          cluster_theta(other.cluster_theta), cluster_phi(other.cluster_phi),
          cluster_ncells(other.cluster_ncells),
          ecal_chi2(other.ecal_chi2),
          ecal_mult(other.ecal_mult), ecal_counter(other.ecal_counter),
          ecal_valid(other.ecal_valid) {}

    // ========================================================================
    // Data Loading from NTupleReader
    // ========================================================================
    
    /**
     * @brief Fill particle from ntuple reader
     * @param reader NTupleReader with loaded event
     * @param index ECAL hit index (1, 2, 3, 4, or 5)
     * @param type Kinematic type to set (default: RECONSTRUCTED)
     * @return true if successfully filled, false if variables missing
     *
     * Reads all neutr_*_N variables where N = index.
     * Uses cluster reconstruction variables (neutr_cluster_*) for PParticle kinematics.
     * For photons (mass=0), uses cluster_energy as momentum.
     * For massive particles, uses neutr_p_N.
     */
    bool setFromReader(NTupleReader& reader, int index,
                       KinematicType type = KinematicType::RECONSTRUCTED) {

        // Build variable name suffix
        std::string idx = "_" + std::to_string(index);

        // Check if cluster variables exist (required for kinematics)
        std::string cluster_theta_var = "neutr_cluster_theta" + idx;
        std::string cluster_phi_var = "neutr_cluster_phi" + idx;
        std::string cluster_energy_var = "neutr_cluster_energy" + idx;
        std::string p_var = "neutr_p" + idx;

        if (!reader.hasVariable(cluster_theta_var) ||
            !reader.hasVariable(cluster_phi_var)) {
            ecal_valid = false;
            return false;
        }

        // Read cluster reconstruction variables (for PParticle kinematics)
        cluster_theta = reader[cluster_theta_var];
        cluster_phi = reader[cluster_phi_var];
        cluster_energy = readIfExists(reader, cluster_energy_var, -1.0);
        cluster_ncells = static_cast<int>(readIfExists(reader, "neutr_cluster_ncells" + idx, -1.0));

        // For momentum: use cluster_energy for photons (mass=0), p for massive
        double p_val;
        double mass_hypo = vec().M();  // Get mass from base class

        if (mass_hypo == 0.0 && cluster_energy > 0) {
            // Photon: p = E (massless)
            p_val = cluster_energy;
        } else if (reader.hasVariable(p_var)) {
            // Massive particle: use momentum
            p_val = reader[p_var];
        } else if (cluster_energy > 0) {
            // Fallback: use cluster energy
            p_val = cluster_energy;
        } else {
            ecal_valid = false;
            return false;
        }

        // Set momentum using cluster angles (inherited from PParticle)
        setFromSpherical(p_val, cluster_theta, cluster_phi, type);

        // Read ECAL-specific fields (primary reconstruction)
        ecal_theta     = readIfExists(reader, "neutr_theta" + idx, -1.0);
        ecal_phi       = readIfExists(reader, "neutr_phi" + idx, -1.0);

        ecal_pid       = static_cast<int>(readIfExists(reader, "neutr_pid" + idx, 0.0));
        ecal_clusterId = static_cast<int>(readIfExists(reader, "neutr_clusterid" + idx, -1.0));
        ecal_charge    = readIfExists(reader, "neutr_q" + idx, 0.0);

        ecal_energy    = readIfExists(reader, "neutr_energy" + idx, -1.0);
        ecal_p_pid     = readIfExists(reader, "neutr_p_pid" + idx, -1.0);

        ecal_mass      = readIfExists(reader, "neutr_mass" + idx, -1.0);
        ecal_mass2     = readIfExists(reader, "neutr_mass2" + idx, -1.0);

        ecal_beta      = readIfExists(reader, "neutr_beta" + idx, -1.0);
        ecal_tof       = readIfExists(reader, "neutr_tof" + idx, -1.0);
        ecal_tofRec    = readIfExists(reader, "neutr_tofrec" + idx, -1.0);
        ecal_dist      = readIfExists(reader, "neutr_dist" + idx, -1.0);

        ecal_r         = readIfExists(reader, "neutr_r" + idx, -1.0);
        ecal_z         = readIfExists(reader, "neutr_z" + idx, -1.0);
        ecal_chi2      = readIfExists(reader, "neutr_chi2" + idx, -1.0);

        ecal_phi2      = readIfExists(reader, "neutr_phi2" + idx, -1.0);
        ecal_theta2    = readIfExists(reader, "neutr_theta2" + idx, -1.0);
        ecal_r2        = readIfExists(reader, "neutr_r2" + idx, -1.0);
        ecal_z2        = readIfExists(reader, "neutr_z2" + idx, -1.0);

        // Multiplicity (same for all hits, no index suffix)
        ecal_mult      = static_cast<int>(readIfExists(reader, "neutr_mult", 0.0));
        ecal_counter   = static_cast<int>(readIfExists(reader, "neutr_counter", 0.0));

        ecal_valid = true;
        return true;
    }
    
    /**
     * @brief Reset ECAL-specific fields to default values
     */
    void resetEcalFields() {
        ecal_pid = 0;
        ecal_clusterId = -1;
        ecal_charge = 0.0;
        ecal_energy = -1.0;
        ecal_p_pid = -1.0;
        ecal_mass = -1.0;
        ecal_mass2 = -1.0;
        ecal_beta = -1.0;
        ecal_tof = -1.0;
        ecal_tofRec = -1.0;
        ecal_dist = -1.0;
        ecal_phi = -1.0;
        ecal_theta = -1.0;
        ecal_r = -1.0;
        ecal_z = -1.0;
        ecal_phi2 = -1.0;
        ecal_theta2 = -1.0;
        ecal_r2 = -1.0;
        ecal_z2 = -1.0;
        cluster_energy = -1.0;
        cluster_theta = -1.0;
        cluster_phi = -1.0;
        cluster_ncells = -1;
        ecal_chi2 = -1.0;
        ecal_mult = 0;
        ecal_counter = 0;
        ecal_valid = false;
    }
    
    // ========================================================================
    // Quality Accessors
    // ========================================================================
    
    /**
     * @brief Check if particle is valid (was successfully filled)
     */
    bool isValid() const { return ecal_valid; }
    
    /**
     * @brief Get cluster ID
     */
    int getClusterId() const { return ecal_clusterId; }
    
    /**
     * @brief Get chi2 quality
     */
    double getChi2() const { return ecal_chi2; }
    
    /**
     * @brief Get distance to EMC
     */
    double getDistToEmc() const { return ecal_dist; }
    
    /**
     * @brief Get measured energy from calorimeter
     */
    double getEnergy() const { return ecal_energy; }
    
    /**
     * @brief Get measured beta
     */
    double getMeasuredBeta() const { return ecal_beta; }
    
    /**
     * @brief Check if passes basic quality cuts
     * @param minEnergy Minimum energy (default: 0)
     * @param maxChi2 Maximum chi2 (default: no cut = 1e9)
     */
    bool passesQuality(double minEnergy = 0.0, double maxChi2 = 1e9) const {
        if (!ecal_valid) return false;
        if (ecal_energy < minEnergy && ecal_energy > 0) return false;
        if (ecal_chi2 > maxChi2 && ecal_chi2 > 0) return false;
        return true;
    }
    
    /**
     * @brief Print ECAL-specific info
     */
    void printEcal() const {
        print();  // Base class print
        std::cout << "  ECAL Info:" << std::endl;
        std::cout << "    pid: " << ecal_pid << ", clusterId: " << ecal_clusterId << std::endl;
        std::cout << "    energy: " << ecal_energy << " MeV" << std::endl;
        std::cout << "    beta: " << ecal_beta << ", chi2: " << ecal_chi2 << std::endl;
        std::cout << "    r: " << ecal_r << " mm, z: " << ecal_z << " mm" << std::endl;
        std::cout << "    distToEmc: " << ecal_dist << " mm" << std::endl;
    }

private:
    /**
     * @brief Helper to read variable if exists, else return default
     */
    double readIfExists(NTupleReader& reader, const std::string& name, double def) {
        return reader.hasVariable(name) ? reader[name] : def;
    }
};

#endif // PPARTICLE_ECAL_H
