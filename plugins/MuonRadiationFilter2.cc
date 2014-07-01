#include "Validation/MCEmbedding/plugins/MuonRadiationFilter2.h"

#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <TMath.h>

// The code in this anonymous namespace is taken from TauAnalysis/MCEmbeddingTools/src/embeddingAuxFunctions.cc.
// It would be better for this package to depend on it, however, the code is only available since CMSSW_7_0_X.
// In order for this package to be able to be used also with CMSSW_5_3_X, the code is duplicated here. Once
// TODO: Once CMSSW_5_3_X becomes irrelevant, we should instead use the code from TauAnalysis/MCEmbeddingTools
// instead of this duplication.
//#include "TauAnalysis/MCEmbeddingTools/interface/embeddingAuxFunctions.h"
namespace
{

bool higherPt(const reco::CandidateBaseRef& muon1, const reco::CandidateBaseRef& muon2)
{
  return (muon1->pt() > muon2->pt());
}

std::vector<reco::CandidateBaseRef> getSelMuons(const edm::Event& evt, const edm::InputTag& srcSelMuons)
{
  std::vector<reco::CandidateBaseRef> selMuons;

  edm::Handle<reco::CompositeCandidateCollection> combCandidatesHandle;
  if ( evt.getByLabel(srcSelMuons, combCandidatesHandle) ) {
    if ( combCandidatesHandle->size() >= 1 ) {
      const reco::CompositeCandidate& combCandidate = combCandidatesHandle->at(0); // TF: use only the first combined candidate
      for ( size_t idx = 0; idx < combCandidate.numberOfDaughters(); ++idx ) { 
	const reco::Candidate* daughter = combCandidate.daughter(idx);
	reco::CandidateBaseRef selMuon;
	if ( daughter->hasMasterClone() ) {
	  selMuon = daughter->masterClone();
	} 
	if ( selMuon.isNull() ) 
	  throw cms::Exception("Configuration") 
	    << "Collection 'selectedMuons' = " << srcSelMuons.label() << " of CompositeCandidates does not refer to daughters of valid type !!\n";
	selMuons.push_back(selMuon);
      }
    }
  } else {
    typedef edm::View<reco::Candidate> CandidateView;
    edm::Handle<CandidateView> selMuonsHandle;
    if ( evt.getByLabel(srcSelMuons, selMuonsHandle) ) {
      for ( size_t idx = 0; idx < selMuonsHandle->size(); ++idx ) {
	selMuons.push_back(reco::CandidateBaseRef(selMuonsHandle->refAt(idx)));
      }
    } else {
      throw cms::Exception("Configuration") 
	<< "Invalid input collection 'selectedMuons' = " << srcSelMuons.label() << " !!\n";
    }
  }

  // sort collection of selected muons by decreasing Pt
  std::sort(selMuons.begin(), selMuons.end(), higherPt);

  return selMuons;
}

reco::CandidateBaseRef getTheMuPlus(const std::vector<reco::CandidateBaseRef>& selMuons)
{
//--- return highest Pt muon of positive charge
//
//    NOTE: function assumes that collection of muons passed as function argument is sorted by decreasing Pt
//         (= as returned by 'getSelMuons' function)
  
  for ( std::vector<reco::CandidateBaseRef>::const_iterator selMuon = selMuons.begin();
	selMuon != selMuons.end(); ++selMuon ) {
    if ( (*selMuon)->charge() > +0.5 ) return (*selMuon);
  }

  // no muon of positive charge found
  return reco::CandidateBaseRef();
}

reco::CandidateBaseRef getTheMuMinus(const std::vector<reco::CandidateBaseRef>& selMuons)
{
//--- return highest Pt muon of negative charge
//
//    NOTE: function assumes that collection of muons passed as function argument is sorted by decreasing Pt
//         (= as returned by 'getSelMuons' function)

  for ( std::vector<reco::CandidateBaseRef>::const_iterator selMuon = selMuons.begin();
	selMuon != selMuons.end(); ++selMuon ) {
    if ( (*selMuon)->charge() < -0.5 ) return (*selMuon);
  }

  // no muon of negative charge found
  return reco::CandidateBaseRef();
}

} // anonymous namespace

const double nomMassZ = 91.1876;

typedef std::pair<const reco::PFCandidate*, const reco::Track*> pfCand_track_pair;

MuonRadiationFilter2::MuonRadiationFilter2(const edm::ParameterSet& cfg)
{
  srcSelectedMuons_ = cfg.getParameter<edm::InputTag>("srcSelectedMuons");
  srcPFCandidates_ = cfg.getParameter<edm::InputTag>("srcPFCandidates");

  seedPtECAL_ = cfg.getParameter<double>("seedPtECAL");
  seedDeltaR_ = cfg.getParameter<double>("seedDeltaR");
  stripCandidatesParticleIds_ = cfg.getParameter<vint>("stripCandidatesParticleIds");
  stripEtaAssociationDistance_ = cfg.getParameter<double>("stripEtaAssociationDistance");
  stripPhiAssociationDistance_ = cfg.getParameter<double>("stripPhiAssociationDistance");

  edm::VParameterSet cfgStripSelection = cfg.getParameter<edm::VParameterSet>("stripSelection");
  for ( edm::VParameterSet::const_iterator cfgStripSelectionEntry = cfgStripSelection.begin();
	cfgStripSelectionEntry != cfgStripSelection.end(); ++cfgStripSelectionEntry ) {
    stripSelectionEntryType stripSelectionEntry;
    stripSelectionEntry.minPt_ = cfgStripSelectionEntry->getParameter<double>("minPt");
    stripSelectionEntry.maxDeltaR_ = cfgStripSelectionEntry->getParameter<double>("maxDeltaR");
    stripSelectionEntry.maxHoverE_ = cfgStripSelectionEntry->getParameter<double>("maxHoverE");
    stripSelectionEntry.applyMassWindowSelection_ = cfgStripSelectionEntry->getParameter<bool>("applyMassWindowSelection");
    stripSelection_.push_back(stripSelectionEntry);    
  }

  dRvetoCone_ = cfg.getParameter<double>("dRvetoCone");
  dRisoCone_ = cfg.getParameter<double>("dRisoCone");
  minTrackPt_ = cfg.getParameter<double>("minTrackPt");
  maxTransverseImpactParameter_ = cfg.getParameter<double>("maxTransverseImpactParameter");
  maxDeltaZ_ = cfg.getParameter<double>("maxDeltaZ");
  maximumSumPtCut_ = cfg.getParameter<double>("maximumSumPtCut");

  invert_ = cfg.getParameter<bool>("invert");
  filter_ = cfg.getParameter<bool>("filter");
  if ( !filter_ ) {
    produces<bool>();
  }

  verbosity_ = ( cfg.exists("verbosity") ) ?
    cfg.getParameter<int>("verbosity") : 0;
}

namespace
{
  const reco::Track* getMuonTrack(const reco::Muon& muon)
  {
    const reco::Track* track = 0;
    if      ( muon.innerTrack().isNonnull()  ) track = muon.innerTrack().get();
    else if ( muon.globalTrack().isNonnull() ) track = muon.globalTrack().get();
    else if ( muon.outerTrack().isNonnull()  ) track = muon.outerTrack().get();
    return track;
  }

  std::string getPFCandidateType(reco::PFCandidate::ParticleType pfCandidateType)
  {
    if      ( pfCandidateType == reco::PFCandidate::X         ) return "undefined";
    else if ( pfCandidateType == reco::PFCandidate::h         ) return "PFChargedHadron";
    else if ( pfCandidateType == reco::PFCandidate::e         ) return "PFElectron";
    else if ( pfCandidateType == reco::PFCandidate::mu        ) return "PFMuon";
    else if ( pfCandidateType == reco::PFCandidate::gamma     ) return "PFGamma";
    else if ( pfCandidateType == reco::PFCandidate::h0        ) return "PFNeutralHadron";
    else if ( pfCandidateType == reco::PFCandidate::h_HF      ) return "HF_had";
    else if ( pfCandidateType == reco::PFCandidate::egamma_HF ) return "HF_em";
    else assert(0);
  }

  void printPFCandidate(const std::string& label, const reco::PFCandidate& pfCandidate, int idx, double* dR, const reco::Candidate::Point* vertex)
  {
    std::cout << label << " #" << idx << " (dR = ";
    if ( dR ) std::cout << (*dR);
    else std::cout << "N/A";
    std::cout << "): Pt = " << pfCandidate.pt() << ","
	      << " eta = " << pfCandidate.eta() << ", phi = " << pfCandidate.phi() << ", mass = " << pfCandidate.mass()
	      << " (type = "<< getPFCandidateType(pfCandidate.particleId()) << ")" << std::endl;
    std::cout << " Calorimeter energy: ECAL = " << pfCandidate.ecalEnergy() << ", HCAL = " << pfCandidate.hcalEnergy() << std::endl;
    const reco::Track* track = 0;
    if ( pfCandidate.trackRef().isNonnull() ) track = &(*pfCandidate.trackRef());
    else if ( pfCandidate.gsfTrackRef().isNonnull() ) track = &(*pfCandidate.gsfTrackRef());
    if ( track ) {
      std::cout << " Track #" << idx << ": Pt = " << track->pt() << ", eta = " << track->eta() << ", phi = " << track->phi() << ", charge = " << track->charge() << std::endl;
      if ( vertex ) {
	std::cout << " (dZ = " << track->dz(*vertex) << ", dXY = " << track->dxy(*vertex) << "," 
		  << " numHits = " << track->hitPattern().numberOfValidTrackerHits() << ", numPxlHits = " << track->hitPattern().numberOfValidPixelHits() << "," 
		  << " chi2 = " << track->normalizedChi2() << ", dPt/Pt = " << (track->ptError()/track->pt()) << ")" << std::endl;
      }
    }
  }

  reco::Candidate::LorentzVector makeCaloEnP4(const reco::Candidate::LorentzVector& refP4, double caloEn)
  {
    double ux = TMath::Cos(refP4.phi())*TMath::Sin(refP4.theta());
    double uy = TMath::Sin(refP4.phi())*TMath::Sin(refP4.theta());
    double uz = TMath::Cos(refP4.theta());
    reco::Candidate::LorentzVector p4(ux*caloEn, uy*caloEn, uz*caloEn, caloEn);
    return p4;
  }

  std::vector<const reco::PFCandidate*> makeStripPFCandidates(const reco::PFCandidate* seed, const std::vector<const reco::PFCandidate*>& pfCandidates,  
							      double stripEtaAssociationDistance, double stripPhiAssociationDistance)
  {
    reco::Candidate::LorentzVector stripP4 = seed->p4();
    std::vector<const reco::PFCandidate*> stripPFCandidates;
    stripPFCandidates.push_back(seed);
    std::set<size_t> pfCandidateIdsCurrentStrip;
    bool isPFCandidateAdded = false;
    int stripBuildIteration = 0;
    const int maxStripBuildIterations = 10;
    do {
      isPFCandidateAdded = false;
      
      size_t numPFCandidates = pfCandidates.size();
      for ( size_t pfCandidateId = 0; pfCandidateId < numPFCandidates; ++pfCandidateId ) {
	const reco::PFCandidate* pfCandidate = pfCandidates[pfCandidateId];
	if ( pfCandidate == seed || pfCandidateIdsCurrentStrip.find(pfCandidateId) != pfCandidateIdsCurrentStrip.end() ) continue; // do not include same pfCandidate twice
	if ( TMath::Abs(stripP4.eta() - pfCandidate->eta()) < stripEtaAssociationDistance &&
	     TMath::Abs(stripP4.phi() - pfCandidate->phi()) < stripPhiAssociationDistance ) {
	  stripP4 += makeCaloEnP4(pfCandidate->p4(), pfCandidate->ecalEnergy());
	  stripPFCandidates.push_back(pfCandidate);
	  isPFCandidateAdded = true;
	  pfCandidateIdsCurrentStrip.insert(pfCandidateId);
	}
      }

      ++stripBuildIteration;
    } while ( isPFCandidateAdded && (stripBuildIteration < maxStripBuildIterations || maxStripBuildIterations == -1) );
    return stripPFCandidates;
  }
  
  reco::Candidate::LorentzVector compStripP4(std::vector<const reco::PFCandidate*>& stripPFCandidates, double& stripPt_ECAL, double& stripPt_HCAL, double& stripEta, double& stripPhi)
  {
    reco::Candidate::LorentzVector stripP4_ECAL;
    reco::Candidate::LorentzVector stripP4_HCAL;
    for ( std::vector<const reco::PFCandidate*>::const_iterator stripPFCandidate = stripPFCandidates.begin();
	  stripPFCandidate != stripPFCandidates.end(); ++stripPFCandidate ) {
      stripP4_ECAL += makeCaloEnP4((*stripPFCandidate)->p4(), (*stripPFCandidate)->ecalEnergy());
      stripP4_HCAL += makeCaloEnP4((*stripPFCandidate)->p4(), (*stripPFCandidate)->hcalEnergy());
    }
    stripPt_ECAL = stripP4_ECAL.pt();
    stripPt_HCAL = stripP4_HCAL.pt();
    stripEta     = stripP4_ECAL.eta();
    stripPhi     = stripP4_ECAL.phi();
    return stripP4_ECAL;
  }

  double compTrackIsoPtSum(std::vector<const reco::PFCandidate*>& stripPFCandidates, double stripEta, double stripPhi, 
			   const std::vector<pfCand_track_pair>& selTracks, double dRvetoCone, double dRisoCone)
  {
    double trackIsoPtSum = 0.;
    for ( std::vector<pfCand_track_pair>::const_iterator selTrack = selTracks.begin();
	  selTrack != selTracks.end(); ++selTrack ) {
      double dR = deltaR(selTrack->first->eta(), selTrack->first->phi(), stripEta, stripPhi);
      if ( dR < dRisoCone ) {
	bool isStrip = false;
	for ( std::vector<const reco::PFCandidate*>::const_iterator stripPFCandidate = stripPFCandidates.begin();
	      stripPFCandidate != stripPFCandidates.end(); ++stripPFCandidate ) {
	  double dR = deltaR(selTrack->first->p4(), (*stripPFCandidate)->p4());
	  if ( dR < dRvetoCone ) {
	    isStrip = true;
	    break;
	  }
	}
	if ( !isStrip ) {
	  trackIsoPtSum += selTrack->second->pt();
	}
      }
    }
    return trackIsoPtSum;
  }
}

bool MuonRadiationFilter2::filter(edm::Event& evt, const edm::EventSetup& es)
{
  if ( verbosity_ >= 1 ) {
    std::cout << "<MuonRadiationFilter2::filter>:" << std::endl;
  }

  std::vector<reco::CandidateBaseRef> selMuons = getSelMuons(evt, srcSelectedMuons_);
  const reco::CandidateBaseRef muPlus  = getTheMuPlus(selMuons);
  const reco::CandidateBaseRef muMinus = getTheMuMinus(selMuons);

  if ( muPlus.isNull() || muMinus.isNull() ) return false; // not selected Z --> mu+ mu- event: reject event

  const reco::Track* muPlusTrack = 0;
  const reco::Track* muMinusTrack = 0;
  if ( dynamic_cast<const reco::Muon*>(muPlus.get()) && dynamic_cast<const reco::Muon*>(muMinus.get()) ) {
    muPlusTrack = getMuonTrack(*dynamic_cast<const reco::Muon*>(muPlus.get()));
    muMinusTrack = getMuonTrack(*dynamic_cast<const reco::Muon*>(muMinus.get()));
  }
  reco::Vertex::Point muPlusVtx = muPlus->vertex();
  reco::Vertex::Point muMinusVtx = muMinus->vertex();
  if ( verbosity_ >= 2 ) {
    std::cout << "mu+: Pt = " << muPlus->pt() << ", eta = " << muPlus->eta() << ", phi = " << muPlus->phi() << std::endl;
    std::cout << " track: Pt = " << muPlusTrack->pt() << ", eta = " << muPlusTrack->eta() << ", phi = " << muPlusTrack->phi() << std::endl;
    std::cout << " vtx: x = " << muPlusVtx.x() << ", y = " << muPlusVtx.y() << ", z = " << muPlusVtx.z() << std::endl;
    std::cout << "mu-: Pt = " << muMinus->pt() << ", eta = " << muMinus->eta() << ", phi = " << muMinus->phi() << std::endl;
    std::cout << " track: Pt = " << muMinusTrack->pt() << ", eta = " << muMinusTrack->eta() << ", phi = " << muMinusTrack->phi() << std::endl;
    std::cout << " vtx: x = " << muMinusVtx.x() << ", y = " << muMinusVtx.y() << ", z = " << muMinusVtx.z() << std::endl;
  }

  double massWithoutStrip = (muPlus->p4() + muMinus->p4()).mass();
  if ( verbosity_ >= 2 ) {
    std::cout << " massWithoutStrip = " << massWithoutStrip << std::endl;
  }

  edm::Handle<PFCandidateView> pfCandidates;
  evt.getByLabel(srcPFCandidates_, pfCandidates);

  std::vector<const reco::PFCandidate*> pfCandidates_passingParticleId;
  int pfCandidateIdx = 0;
  for ( PFCandidateView::const_iterator pfCandidate = pfCandidates->begin();
	pfCandidate != pfCandidates->end(); ++pfCandidate, ++pfCandidateIdx ) {
    reco::PFCandidate::ParticleType particleId = pfCandidate->particleId();
    bool passesParticleId = false;    
    for ( vint::const_iterator stripCandidatesParticleId = stripCandidatesParticleIds_.begin();
	  stripCandidatesParticleId != stripCandidatesParticleIds_.end(); ++stripCandidatesParticleId ) {
      if ( particleId == (*stripCandidatesParticleId) ) {
	passesParticleId = true;
	break;
      }
    }
    if ( verbosity_ >= 3 && pfCandidate->pt() > 1.0 ) {
      printPFCandidate("PFCandidate", *pfCandidate, pfCandidateIdx, 0, 0);
      std::cout << "particleId = " << particleId << " (passesParticleId = " << passesParticleId << ")" << std::endl;
    }
    if ( passesParticleId ) {
      if ( verbosity_ >= 3 ) {
	std::cout << "--> storing PFCandidate #" << pfCandidateIdx << " as seed #" << pfCandidates_passingParticleId.size() << "." << std::endl;
      }
      pfCandidates_passingParticleId.push_back(&(*pfCandidate));
    }
  }

  std::vector<pfCand_track_pair> selTracks;
  for ( PFCandidateView::const_iterator pfCandidate = pfCandidates->begin();
	pfCandidate != pfCandidates->end(); ++pfCandidate ) {
    const reco::Track* track = 0;
    if ( pfCandidate->muonRef().isNonnull() ) {
      const reco::Muon* muon = pfCandidate->muonRef().get();
      track = getMuonTrack(*muon);
    }
    if ( !track ) {
      if      ( pfCandidate->gsfTrackRef().isNonnull() ) track = pfCandidate->gsfTrackRef().get();
      else if ( pfCandidate->trackRef().isNonnull()    ) track = pfCandidate->trackRef().get();
    }
    if ( !track ) continue;

    if ( !(track->pt() > minTrackPt_) ) continue;
    
    double dRmuPlus  = ( muPlusTrack  ) ? deltaR(track->eta(), track->phi(), muPlusTrack->eta(),  muPlusTrack->phi())  : 1.e+3;
    double dRmuMinus = ( muMinusTrack ) ? deltaR(track->eta(), track->phi(), muMinusTrack->eta(), muMinusTrack->phi()) : 1.e+3;
    if ( dRmuPlus < dRvetoCone_ || dRmuMinus < dRvetoCone_ ) continue;
    
    double dZmuPlus   = track->dz(muPlusVtx);
    double dXYmuPlus  = track->dxy(muPlusVtx);
    double dZmuMinus  = track->dz(muMinusVtx);
    double dXYmuMinus = track->dxy(muMinusVtx);
    if ( (TMath::Abs(dZmuPlus)  < maxDeltaZ_ && TMath::Abs(dXYmuPlus)  < maxTransverseImpactParameter_) ||
	 (TMath::Abs(dZmuMinus) < maxDeltaZ_ && TMath::Abs(dXYmuMinus) < maxTransverseImpactParameter_) ) {
      selTracks.push_back(pfCand_track_pair(&(*pfCandidate), track));
    }
  }

  bool isMuonRadiation = false;

  int seedIdx = 0;
  for ( std::vector<const reco::PFCandidate*>::const_iterator seed = pfCandidates_passingParticleId.begin();
	seed != pfCandidates_passingParticleId.end(); ++seed, ++seedIdx ) {
    double dRmuPlus_seed  = deltaR((*seed)->p4(), muPlus->p4());
    double dRmuMinus_seed = deltaR((*seed)->p4(), muMinus->p4());
    if ( verbosity_ >= 2 && (*seed)->pt() > 1.0 ) {
      double dR = 1.e+3;
      reco::Candidate::Point vertex;
      if ( dRmuPlus_seed < dRmuMinus_seed ) {
	dR = dRmuPlus_seed;
	vertex = muPlusVtx;
      } else {
	dR = dRmuMinus_seed;
	vertex = muMinusVtx;
      }
      printPFCandidate("seed", **seed, seedIdx, &dR, &vertex);
      std::cout << "dR(mu+) = " << dRmuPlus_seed << ", dR(mu-) = " << dRmuMinus_seed << std::endl;
    }
    if ( !(dRmuPlus_seed < seedDeltaR_ || dRmuMinus_seed < seedDeltaR_) ) continue;
    
    double seedPt_ECAL = TMath::Sin((*seed)->theta())*(*seed)->ecalEnergy();
    if ( verbosity_ >= 2 ) {
      std::cout << "seedPt_ECAL = " << seedPt_ECAL << std::endl;
    }
    if ( !(seedPt_ECAL > seedPtECAL_) ) continue;
      
    std::vector<const reco::PFCandidate*> stripPFCandidates = makeStripPFCandidates(*seed, pfCandidates_passingParticleId, stripEtaAssociationDistance_, stripPhiAssociationDistance_);   
    double stripPt_ECAL, stripPt_HCAL, stripEta, stripPhi;
    reco::Candidate::LorentzVector stripP4 = compStripP4(stripPFCandidates, stripPt_ECAL, stripPt_HCAL, stripEta, stripPhi);
    double dRmuPlus_strip  = deltaR(stripEta, stripPhi, muPlus->eta(), muPlus->phi());
    double dRmuMinus_strip = deltaR(stripEta, stripPhi, muMinus->eta(), muMinus->phi());
    if ( verbosity_ >= 2 ) {
      std::cout << "strip: Pt(ECAL) = " << stripPt_ECAL << ", Pt(HCAL) = " << stripPt_HCAL << ", eta = " << stripEta << ", phi = " << stripPhi << std::endl;
      std::cout << " dR(mu+) = " << dRmuPlus_strip << ", dR(mu-) = " << dRmuMinus_strip << std::endl;
    }

    double trackIsoPtSum = compTrackIsoPtSum(stripPFCandidates, stripEta, stripPhi, selTracks, dRvetoCone_, dRisoCone_);
    bool isIsolated = (trackIsoPtSum < maximumSumPtCut_);
    if ( verbosity_ >= 2 ) {
      std::cout << "trackIsoPtSum = " << trackIsoPtSum << " (isIsolated = " << isIsolated << ")" << std::endl;
    }

    for ( std::vector<stripSelectionEntryType>::const_iterator stripSelectionEntry = stripSelection_.begin();
	  stripSelectionEntry != stripSelection_.end(); ++stripSelectionEntry ) {
      if ( stripPt_ECAL > stripSelectionEntry->minPt_ && 
	   (dRmuPlus_strip < stripSelectionEntry->maxDeltaR_ || dRmuMinus_strip < stripSelectionEntry->maxDeltaR_) &&
	   stripPt_HCAL < (stripSelectionEntry->maxHoverE_*stripPt_ECAL) ) {
	double massWithStrip = (muPlus->p4() + muMinus->p4() + stripP4).mass();
	if ( verbosity_ >= 2 ) {
	  std::cout << "massWithStrip = " << massWithStrip << std::endl;
	}
	if ( TMath::Abs(massWithStrip - nomMassZ) < TMath::Abs(massWithoutStrip - nomMassZ) || !stripSelectionEntry->applyMassWindowSelection_ ) {
	  if ( verbosity_ >= 2 ) {
	    std::cout << "--> setting isMuonRadiation = true." << std::endl;
	  }
	  isMuonRadiation = true;
	}
      }	   
    }
  }

  if ( verbosity_ >= 1 ) {
    std::cout << "isMuonRadiation = " << isMuonRadiation << std::endl;
  }

  if ( filter_ ) {
    if ( invert_ != isMuonRadiation ) return false; // reject events with muon -> muon + photon radiation
    else return true;
  } else {
    std::auto_ptr<bool> filter_result(new bool(invert_ != !isMuonRadiation));
    evt.put(filter_result);
    return true;
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(MuonRadiationFilter2);
