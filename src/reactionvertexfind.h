//////////////////////////////////////////////////////////////////////////////
// ReactionVertexFind - Flexible vertex finder for N tracks (N >= 2)
//
// Supports: 2-track (p, pip), 3-track, 4-track, or more
//////////////////////////////////////////////////////////////////////////////

#ifndef REACTIONVERTEXFIND_H
#define REACTIONVERTEXFIND_H

#include "TMath.h"
#include "Math/Vector3D.h"
#include <vector>

using ROOT::Math::XYZVector;

class ReactionVertexFind {

public:
    // Track info structure
    struct TrackInfo {
        Double_t r, z, theta, phi;
        Double_t rkchi2, mdcinnerchi2, beta;
        XYZVector pos;
        XYZVector dir;
        Double_t weight;
        Double_t distance;
        Bool_t   isGood;
    };
    
    // Quality result structure
    struct VertexQuality {
        Float_t chi2;
        Float_t rVertex;
        Float_t zSpread;
        Float_t sumOfWeights;
        Int_t   nTracks;
        Int_t   nGoodTracks;
        Int_t   nConvergence;
        std::vector<Float_t> trackDCA;
        std::vector<Float_t> trackWeight;
    };

private:
    std::vector<TrackInfo> fTracks;
    XYZVector fVertex;
    VertexQuality fQuality;
    
    // Cut parameters
    Float_t fMinZ, fMaxZ, fMaxR;
    Float_t fMaxRKChi2, fMaxSegChi2;
    Float_t fMinBeta, fMaxBeta;
    
    // Fit parameters
    Bool_t  fUseTukey;
    Float_t fTukeyConst;
    Float_t fEpsilon;
    Int_t   fMaxIterations;

public:
    ReactionVertexFind() :
        fVertex(-1000., -1000., -1000.),
        fMinZ(-120.), fMaxZ(20.), fMaxR(15.),
        fMaxRKChi2(100.), fMaxSegChi2(6.),
        fMinBeta(0.), fMaxBeta(1.2),
        fUseTukey(true), fTukeyConst(6.0),
        fEpsilon(0.3), fMaxIterations(100)
    {
        reset();
    }
    
    //========================================================================
    // Reset - clear all tracks
    //========================================================================
    void reset() {
        fTracks.clear();
        fVertex.SetXYZ(-1000., -1000., -1000.);
        fQuality.chi2 = -1.;
        fQuality.rVertex = -1.;
        fQuality.zSpread = -1.;
        fQuality.sumOfWeights = 0.;
        fQuality.nTracks = 0;
        fQuality.nGoodTracks = 0;
        fQuality.nConvergence = -1;
        fQuality.trackDCA.clear();
        fQuality.trackWeight.clear();
    }
    
    //========================================================================
    // Add a single track
    //========================================================================
    void addTrack(Float_t r, Float_t z, Float_t theta, Float_t phi,
                  Float_t rkchi2, Float_t mdcinnerchi2, Float_t beta) {
        TrackInfo track;
        track.r = r;
        track.z = z;
        track.theta = theta;
        track.phi = phi;
        track.rkchi2 = rkchi2;
        track.mdcinnerchi2 = mdcinnerchi2;
        track.beta = beta;
        track.weight = 1.0;
        track.distance = 0.;
        track.isGood = checkTrackQuality(track);
        
        if (track.isGood) {
            transformTrack(track);
        }
        
        fTracks.push_back(track);
    }
    
    //========================================================================
    // Convenience: set 2 tracks (p, pip)
    //========================================================================
    void setTracks2(
        Float_t p_r, Float_t p_z, Float_t p_theta, Float_t p_phi,
        Float_t p_rkchi2, Float_t p_mdcinnerchi2, Float_t p_beta,
        Float_t pip_r, Float_t pip_z, Float_t pip_theta, Float_t pip_phi,
        Float_t pip_rkchi2, Float_t pip_mdcinnerchi2, Float_t pip_beta
    ) {
        reset();
        addTrack(p_r, p_z, p_theta, p_phi, p_rkchi2, p_mdcinnerchi2, p_beta);
        addTrack(pip_r, pip_z, pip_theta, pip_phi, pip_rkchi2, pip_mdcinnerchi2, pip_beta);
    }
    
    //========================================================================
    // Convenience: set 3 tracks
    //========================================================================
    void setTracks3(
        Float_t t1_r, Float_t t1_z, Float_t t1_theta, Float_t t1_phi,
        Float_t t1_rkchi2, Float_t t1_mdcinnerchi2, Float_t t1_beta,
        Float_t t2_r, Float_t t2_z, Float_t t2_theta, Float_t t2_phi,
        Float_t t2_rkchi2, Float_t t2_mdcinnerchi2, Float_t t2_beta,
        Float_t t3_r, Float_t t3_z, Float_t t3_theta, Float_t t3_phi,
        Float_t t3_rkchi2, Float_t t3_mdcinnerchi2, Float_t t3_beta
    ) {
        reset();
        addTrack(t1_r, t1_z, t1_theta, t1_phi, t1_rkchi2, t1_mdcinnerchi2, t1_beta);
        addTrack(t2_r, t2_z, t2_theta, t2_phi, t2_rkchi2, t2_mdcinnerchi2, t2_beta);
        addTrack(t3_r, t3_z, t3_theta, t3_phi, t3_rkchi2, t3_mdcinnerchi2, t3_beta);
    }
    
    //========================================================================
    // Convenience: set 4 tracks (pim, pip, em, ep)
    //========================================================================
    void setTracks4(
        Float_t t1_r, Float_t t1_z, Float_t t1_theta, Float_t t1_phi,
        Float_t t1_rkchi2, Float_t t1_mdcinnerchi2, Float_t t1_beta,
        Float_t t2_r, Float_t t2_z, Float_t t2_theta, Float_t t2_phi,
        Float_t t2_rkchi2, Float_t t2_mdcinnerchi2, Float_t t2_beta,
        Float_t t3_r, Float_t t3_z, Float_t t3_theta, Float_t t3_phi,
        Float_t t3_rkchi2, Float_t t3_mdcinnerchi2, Float_t t3_beta,
        Float_t t4_r, Float_t t4_z, Float_t t4_theta, Float_t t4_phi,
        Float_t t4_rkchi2, Float_t t4_mdcinnerchi2, Float_t t4_beta
    ) {
        reset();
        addTrack(t1_r, t1_z, t1_theta, t1_phi, t1_rkchi2, t1_mdcinnerchi2, t1_beta);
        addTrack(t2_r, t2_z, t2_theta, t2_phi, t2_rkchi2, t2_mdcinnerchi2, t2_beta);
        addTrack(t3_r, t3_z, t3_theta, t3_phi, t3_rkchi2, t3_mdcinnerchi2, t3_beta);
        addTrack(t4_r, t4_z, t4_theta, t4_phi, t4_rkchi2, t4_mdcinnerchi2, t4_beta);
    }
    
    //========================================================================
    // Configure cuts
    //========================================================================
    void setCuts(Float_t minZ, Float_t maxZ, Float_t maxR) {
        fMinZ = minZ; fMaxZ = maxZ; fMaxR = maxR;
    }
    
    void setQualityCuts(Float_t maxRKChi2, Float_t maxSegChi2,
                        Float_t minBeta, Float_t maxBeta) {
        fMaxRKChi2 = maxRKChi2; fMaxSegChi2 = maxSegChi2;
        fMinBeta = minBeta; fMaxBeta = maxBeta;
    }
    
    void useTukeyWeights(Bool_t use) { fUseTukey = use; }
    void setTukeyConst(Float_t c)    { fTukeyConst = c; }
    void setEpsilon(Float_t eps)     { fEpsilon = eps; }
    void setMaxIterations(Int_t n)   { fMaxIterations = n; }
    
    //========================================================================
    // Main vertex finding
    //========================================================================
    Bool_t findVertex() {
        // Reset quality
        fVertex.SetXYZ(-1000., -1000., -1000.);
        fQuality.chi2 = -1.;
        fQuality.nConvergence = -1;
        fQuality.nTracks = fTracks.size();
        fQuality.nGoodTracks = 0;
        fQuality.trackDCA.resize(fTracks.size(), -1.);
        fQuality.trackWeight.resize(fTracks.size(), 0.);
        
        // Count good tracks
        for (size_t i = 0; i < fTracks.size(); i++) {
            if (fTracks[i].isGood) fQuality.nGoodTracks++;
        }
        
        // Need at least 2 good tracks
        if (fQuality.nGoodTracks < 2) return false;
        
        // Initial fit
        if (!solveLSM()) return false;
        
        if (!fUseTukey) {
            calcQualityParameters();
            fQuality.nConvergence = 1;
            return true;
        }
        
        return doTukeyIteration();
    }
    
    //========================================================================
    // Simple Z vertex (weighted average, fast)
    //========================================================================
    Float_t calcSimpleZVertex() const {
        Double_t sumWZ = 0., sumW = 0.;
        for (size_t i = 0; i < fTracks.size(); i++) {
            if (!fTracks[i].isGood) continue;
            Double_t tanTh = TMath::Tan(fTracks[i].theta * TMath::DegToRad());
            Double_t w = tanTh * tanTh;
            sumWZ += w * fTracks[i].z;
            sumW += w;
        }
        return (sumW > 0) ? sumWZ / sumW : -1000.;
    }
    
    //========================================================================
    // Get results
    //========================================================================
    void getVertex(Float_t& vx, Float_t& vy, Float_t& vz) const {
        vx = fVertex.X(); vy = fVertex.Y(); vz = fVertex.Z();
    }
    Float_t getVx() const { return fVertex.X(); }
    Float_t getVy() const { return fVertex.Y(); }
    Float_t getVz() const { return fVertex.Z(); }
    const XYZVector& getVertexVector() const { return fVertex; }
    
    const VertexQuality& getQuality() const { return fQuality; }
    Float_t getChi2()         const { return fQuality.chi2; }
    Float_t getRVertex()      const { return fQuality.rVertex; }
    Float_t getZSpread()      const { return fQuality.zSpread; }
    Float_t getSumOfWeights() const { return fQuality.sumOfWeights; }
    Int_t   getNTracks()      const { return fQuality.nTracks; }
    Int_t   getNGoodTracks()  const { return fQuality.nGoodTracks; }
    Int_t   getIterations()   const { return fQuality.nConvergence; }
    
    Float_t getTrackDCA(Int_t i) const { 
        return (i >= 0 && i < (Int_t)fQuality.trackDCA.size()) ? fQuality.trackDCA[i] : -1; 
    }
    Float_t getTrackWeight(Int_t i) const { 
        return (i >= 0 && i < (Int_t)fQuality.trackWeight.size()) ? fQuality.trackWeight[i] : -1; 
    }
    Bool_t isTrackGood(Int_t i) const { 
        return (i >= 0 && i < (Int_t)fTracks.size()) ? fTracks[i].isGood : false; 
    }
    
    Bool_t isGoodVertex(Float_t maxChi2 = 2., Float_t maxRVertex = 3., 
                        Float_t minWeightSum = -1.) const {
        if (fVertex.X() < -999.) return false;
        if (fQuality.chi2 < 0 || fQuality.chi2 > maxChi2) return false;
        if (fQuality.rVertex > maxRVertex) return false;
        // For 2 tracks, expect sum ~2; for 4 tracks, expect sum ~4
        Float_t expectedMinWeight = (minWeightSum < 0) ? fQuality.nGoodTracks * 0.75 : minWeightSum;
        if (fQuality.sumOfWeights < expectedMinWeight) return false;
        return true;
    }

private:
    Bool_t checkTrackQuality(const TrackInfo& track) const {
        if (track.z < fMinZ || track.z > fMaxZ) return false;
        if (track.r > fMaxR) return false;
        if (track.rkchi2 > fMaxRKChi2 && track.rkchi2 > 0) return false;
        if (track.mdcinnerchi2 > fMaxSegChi2 && track.mdcinnerchi2 > 0) return false;
        if (track.beta < fMinBeta || track.beta > fMaxBeta) return false;
        return true;
    }
    
    void transformTrack(TrackInfo& track) {
        Double_t phi_rad = track.phi * TMath::DegToRad();
        Double_t theta_rad = track.theta * TMath::DegToRad();
        Double_t pi2 = TMath::Pi() / 2.;
        
        track.pos.SetXYZ(
            track.r * TMath::Cos(phi_rad + pi2),
            track.r * TMath::Sin(phi_rad + pi2),
            track.z
        );
        
        track.dir.SetXYZ(
            TMath::Sin(theta_rad) * TMath::Cos(phi_rad),
            TMath::Sin(theta_rad) * TMath::Sin(phi_rad),
            TMath::Cos(theta_rad)
        );
    }
    
    Double_t calcTrackDistance(const TrackInfo& track, const XYZVector& vertex) const {
        XYZVector delta = track.pos - vertex;
        XYZVector cross = delta.Cross(track.dir);
        return TMath::Sqrt(cross.Mag2());
    }
    
    Bool_t solveLSM() {
        // Build 3x3 matrix A and vector b using direct arrays
        Double_t A[3][3] = {{0,0,0}, {0,0,0}, {0,0,0}};
        Double_t b[3] = {0, 0, 0};
        
        Int_t nContributing = 0;
        
        for (size_t i = 0; i < fTracks.size(); i++) {
            if (!fTracks[i].isGood || fTracks[i].weight <= 0) continue;
            
            Double_t w  = fTracks[i].weight;
            Double_t ax = fTracks[i].dir.X();
            Double_t ay = fTracks[i].dir.Y();
            Double_t az = fTracks[i].dir.Z();
            Double_t rx = fTracks[i].pos.X();
            Double_t ry = fTracks[i].pos.Y();
            Double_t rz = fTracks[i].pos.Z();
            
            Double_t M00 = w * (ay*ay + az*az);
            Double_t M11 = w * (ax*ax + az*az);
            Double_t M22 = w * (ax*ax + ay*ay);
            Double_t M01 = w * (-ax*ay);
            Double_t M02 = w * (-ax*az);
            Double_t M12 = w * (-ay*az);
            
            A[0][0] += M00;  A[0][1] += M01;  A[0][2] += M02;
            A[1][0] += M01;  A[1][1] += M11;  A[1][2] += M12;
            A[2][0] += M02;  A[2][1] += M12;  A[2][2] += M22;
            
            b[0] += M00*rx + M01*ry + M02*rz;
            b[1] += M01*rx + M11*ry + M12*rz;
            b[2] += M02*rx + M12*ry + M22*rz;
            
            nContributing++;
        }
        
        if (nContributing < 2) {
            fVertex.SetXYZ(-1000., -1000., -1000.);
            return false;
        }
        
        // Add small regularization (Tikhonov) to prevent singularity
        Double_t reg = 1e-6;
        A[0][0] += reg;
        A[1][1] += reg;
        A[2][2] += reg;
        
        // Calculate determinant
        Double_t det = A[0][0] * (A[1][1]*A[2][2] - A[1][2]*A[2][1])
                     - A[0][1] * (A[1][0]*A[2][2] - A[1][2]*A[2][0])
                     + A[0][2] * (A[1][0]*A[2][1] - A[1][1]*A[2][0]);
        
        if (TMath::Abs(det) < 1e-10) {
            fVertex.SetXYZ(-1000., -1000., -1000.);
            return false;
        }
        
        // Direct 3x3 inversion (Cramer's rule)
        Double_t invDet = 1.0 / det;
        
        Double_t inv[3][3];
        inv[0][0] = (A[1][1]*A[2][2] - A[1][2]*A[2][1]) * invDet;
        inv[0][1] = (A[0][2]*A[2][1] - A[0][1]*A[2][2]) * invDet;
        inv[0][2] = (A[0][1]*A[1][2] - A[0][2]*A[1][1]) * invDet;
        inv[1][0] = (A[1][2]*A[2][0] - A[1][0]*A[2][2]) * invDet;
        inv[1][1] = (A[0][0]*A[2][2] - A[0][2]*A[2][0]) * invDet;
        inv[1][2] = (A[0][2]*A[1][0] - A[0][0]*A[1][2]) * invDet;
        inv[2][0] = (A[1][0]*A[2][1] - A[1][1]*A[2][0]) * invDet;
        inv[2][1] = (A[0][1]*A[2][0] - A[0][0]*A[2][1]) * invDet;
        inv[2][2] = (A[0][0]*A[1][1] - A[0][1]*A[1][0]) * invDet;
        
        Double_t vx = inv[0][0]*b[0] + inv[0][1]*b[1] + inv[0][2]*b[2];
        Double_t vy = inv[1][0]*b[0] + inv[1][1]*b[1] + inv[1][2]*b[2];
        Double_t vz = inv[2][0]*b[0] + inv[2][1]*b[1] + inv[2][2]*b[2];
        
        fVertex.SetXYZ(vx, vy, vz);
        return true;
    }
    
    Bool_t doTukeyIteration() {
        Double_t sumResiduals = 0.;
        Int_t nTracks = 0;
        
        for (size_t i = 0; i < fTracks.size(); i++) {
            if (!fTracks[i].isGood) continue;
            Double_t d = calcTrackDistance(fTracks[i], fVertex);
            fTracks[i].distance = d;
            sumResiduals += d * d;
            nTracks++;
        }
        
        if (nTracks == 0) return false;
        
        Double_t width = fTukeyConst * TMath::Sqrt(sumResiduals / nTracks);
        XYZVector oldVertex;
        Int_t iteration = 0;
        
        do {
            oldVertex = fVertex;
            iteration++;
            
            sumResiduals = 0.;
            Double_t sumWeights = 0.;
            
            for (size_t i = 0; i < fTracks.size(); i++) {
                if (!fTracks[i].isGood) {
                    fTracks[i].weight = 0.;
                    continue;
                }
                
                Double_t d = calcTrackDistance(fTracks[i], oldVertex);
                fTracks[i].distance = d;
                
                Double_t w = 0.;
                if (d < width && width > 0) {
                    Double_t u = d / width;
                    w = (1. - u*u) * (1. - u*u);
                }
                
                fTracks[i].weight = w;
                sumResiduals += w * d * d;
                sumWeights += w;
            }
            
            if (sumWeights > 0) {
                width = fTukeyConst * TMath::Sqrt(sumResiduals / sumWeights);
            }
            
            if (!solveLSM()) {
                fQuality.nConvergence = -1;
                return false;
            }
            
            XYZVector dV = fVertex - oldVertex;
            if (TMath::Sqrt(dV.Mag2()) < fEpsilon) {
                fQuality.nConvergence = iteration;
                break;
            }
            
        } while (iteration < fMaxIterations);
        
        if (iteration >= fMaxIterations) {
            fQuality.nConvergence = -2;
        }
        
        calcQualityParameters();
        return true;
    }
    
    void calcQualityParameters() {
        Double_t sumWD2 = 0., sumW = 0.;
        Double_t sumZ = 0., sumZ2 = 0.;
        Int_t n = 0;
        
        fQuality.trackDCA.resize(fTracks.size());
        fQuality.trackWeight.resize(fTracks.size());
        
        for (size_t i = 0; i < fTracks.size(); i++) {
            Double_t d = calcTrackDistance(fTracks[i], fVertex);
            fTracks[i].distance = d;
            
            fQuality.trackDCA[i] = d;
            fQuality.trackWeight[i] = fTracks[i].weight;
            
            if (fTracks[i].isGood) {
                sumWD2 += fTracks[i].weight * d * d;
                sumW += fTracks[i].weight;
                sumZ += fTracks[i].z;
                sumZ2 += fTracks[i].z * fTracks[i].z;
                n++;
            }
        }
        
        fQuality.chi2 = (sumW > 0) ? sumWD2 / sumW : -1.;
        fQuality.sumOfWeights = sumW;
        fQuality.rVertex = TMath::Sqrt(fVertex.X()*fVertex.X() + 
                                        fVertex.Y()*fVertex.Y());
        
        if (n > 1) {
            Double_t meanZ = sumZ / n;
            fQuality.zSpread = TMath::Sqrt(sumZ2/n - meanZ*meanZ);
        } else {
            fQuality.zSpread = 0.;
        }
    }
};

#endif // REACTIONVERTEXFIND_H
