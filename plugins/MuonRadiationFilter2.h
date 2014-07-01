#ifndef TauAnalysis_MCEmbeddingTools_MuonRadiationFilter2_h
#define TauAnalysis_MCEmbeddingTools_MuonRadiationFilter2_h

/** \class MuonRadiationFilter2
 *
 * Veto events in which a muon from Z --> mu+ mu- decay radiates a photon:
 *  muon -> muon + photon
 * 
 * \author Christian Veelken, LLR
 *
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/Common/interface/View.h"

class MuonRadiationFilter2 : public edm::EDFilter 
{
 public:
  explicit MuonRadiationFilter2(const edm::ParameterSet&);
  ~MuonRadiationFilter2() {}

 private:
  bool filter(edm::Event&, const edm::EventSetup&);

  typedef edm::View<reco::PFCandidate> PFCandidateView;

  double compCaloEnECAL(const reco::Candidate::LorentzVector&, const PFCandidateView&);
  void compPFIso_raw(const reco::Candidate::LorentzVector&, const PFCandidateView&, 
		     const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&, 
		     double&, double&, double&);
  double compPFIso_puCorr(const reco::Candidate::LorentzVector&, const PFCandidateView&, const PFCandidateView&, 
			  const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&);
  bool checkMuonRadiation(const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector*, double, const PFCandidateView&, const PFCandidateView&, 
			  const reco::Candidate&, const reco::Candidate&);

  edm::InputTag srcSelectedMuons_;
  edm::InputTag srcPFCandidates_;

  // parameters for reconstruction of eta x phi strips                                
  double seedPtECAL_;
  double seedDeltaR_;
  typedef std::vector<int> vint;
  vint stripCandidatesParticleIds_;
  double stripEtaAssociationDistance_;
  double stripPhiAssociationDistance_;

  // selection of muon FSR events                    
  struct stripSelectionEntryType
  {
    double minPt_;
    double maxDeltaR_;
    double maxHoverE_;
    bool applyMassWindowSelection_;
  };
  std::vector<stripSelectionEntryType> stripSelection_;

  // track isolation parameters
  double dRvetoCone_;
  double dRisoCone_;
  double minTrackPt_;
  double maxTransverseImpactParameter_;
  double maxDeltaZ_;
  double maximumSumPtCut_;

  // global flags
  bool invert_;
  bool filter_;

  int verbosity_;
};

#endif
