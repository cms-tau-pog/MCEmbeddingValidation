// -*- C++ -*-
//
// Package:    MuMuNTupleProducer
// Class:      MuMuNTupleProducer
// 
/**\class MuMuNTupleProducer MuMuNTupleProducer.cc DesyHTauTau/MuMuNTupleProducer/src/MuMuNTupleProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Armin Burgmeier,,,KIT
//         Created:  Mon Jan  9 11:36:09 CET 2012
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "DataFormats/RecoCandidate/interface/IsoDepositVetos.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "SimDataFormats/GeneratorProducts/interface/GenFilterInfo.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "TMath.h"
#include "TTree.h"
#include "TH1F.h"
#include "Math/VectorUtil.h"

#include <TVector3.h>

//
// class declaration
//

static const unsigned int NO_TRIGGER = UINT_MAX;

class MuMuNTupleProducer : public edm::EDAnalyzer
{
public:
  explicit MuMuNTupleProducer(const edm::ParameterSet&);
  ~MuMuNTupleProducer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  // ----------member data ---------------------------
  bool writeAllEvents;
  bool isEmbedded;
  bool isRHEmbedded;
  bool isGenEmbedded;
  bool isData;
  bool isSS;

  edm::InputTag triggerSource;
  edm::InputTag triggerEventSource;
  edm::InputTag origTriggerSource;

  edm::InputTag jetSource;
  edm::InputTag btagSource;
  edm::InputTag trackSource;

  edm::InputTag muonIsoDepsChargedParticlesSource;
  edm::InputTag muonIsoDepsChargedHadronsSource;
  edm::InputTag muonIsoDepsNeutralHadronsSource;
  edm::InputTag muonIsoDepsPhotonsSource;
  edm::InputTag muonIsoDepsPUSource;

  edm::InputTag muonSource;
  edm::InputTag origMuonSource;
  edm::InputTag rhoNeutralSource;
  edm::InputTag origRhoNeutralSource;

  edm::InputTag caloMetSource;
  edm::InputTag pfMetSource;
  edm::InputTag pfType1CorrectedMetSource;
  edm::InputTag mvaPfMetSource;
  edm::InputTag origCaloMetSource;
  edm::InputTag origPfMetSource;

  edm::InputTag vertexSource;
  edm::InputTag pileupSummaryInfoSource;
  edm::InputTag origVertexSource;

  edm::InputTag genParticlesSource;
  edm::InputTag origGenParticlesSource;

  HLTConfigProvider hltConfiguration;

  unsigned int hltIsoMu17Index;
  unsigned int hltIsoMu24Index;
  unsigned int hltIsoMu24Eta2p1Index;
  unsigned int hltDoubleMu7Index;
  unsigned int hltMu13Mu8Index;
  unsigned int hltMu17Mu8Index;

  std::string hltIsoMu17Filter;
  std::string hltIsoMu24Filter;
  std::string hltIsoMu24Eta2p1Filter;
  std::string hltDoubleMu7Filter;
  std::string hltMu13Mu8Filter8;
  std::string hltMu13Mu8Filter13;
  std::string hltMu17Mu8Filter8;
  std::string hltMu17Mu8Filter17;

  TTree* tree;
  TH1F* h_nPV;
  TH1F* h_nTrueInteractions;
  TH1F* h_nOrigPV;
  TH1F* h_dilepmass;
  TH1F* h_origdilepmass;

  unsigned int run;
  unsigned int lumi;
  unsigned int event;

  float minVisPtFilterWeight;
  float zmumuEvtSelEffCorrWeight;
  bool muonRadiationFilter;
  bool muonRadiationFilter2Sel1;
  bool muonRadiationFilter2Sel2;
  bool muonRadiationFilter2Sel3;

  unsigned int nPV;
  float nTrueInteractions;
  unsigned int nGenMEJets;
  unsigned int nOrigPV;

  bool hltIsoMu17;
  bool hltIsoMu24;
  bool hltIsoMu24Eta2p1;
  bool hltDoubleMu7;
  bool hltMu13Mu8;
  bool hltMu17Mu8;
  bool origHltIsoMu17;
  bool origHltDoubleMu7;
  bool origHltMu13Mu8;
  bool origHltMu17Mu8;

  float rho;
  float origRho;

  float caloMet;
  float caloMetPhi;
  float caloMetSigXX;
  float caloMetSigXY;
  float caloMetSigYX;
  float caloMetSigYY;
  float pfMet;
  float pfMetPhi;
  float pfMetSigXX;
  float pfMetSigXY;
  float pfMetSigYX;
  float pfMetSigYY;
  float pfType1CorrectedMet;
  float pfType1CorrectedMetPhi;
  float pfType1CorrectedMetSigXX;
  float pfType1CorrectedMetSigXY;
  float pfType1CorrectedMetSigYX;
  float pfType1CorrectedMetSigYY;
  float mvaPfMet;
  float mvaPfMetPhi;
  float mvaPfMetSigXX;
  float mvaPfMetSigXY;
  float mvaPfMetSigYX;
  float mvaPfMetSigYY;

  float origCaloMet;
  float origCaloMetPhi;
  float origPfMet;
  float origPfMetPhi;

  float posMuonPt;
  float posMuonEta;
  float posMuonPhi;
  float posMuonE;
  float posMuonTrackerIso03;
  float posMuonEcalIso03;
  float posMuonHcalIso03;
  float posMuonChargedHadronPfIso04;
  float posMuonChargedParticlePtPfIso04;
  float posMuonNeutralHadronEtPfIso04;
  float posMuonPhotonEtPfIso04;
  float posMuonPUPtPfIso04;
  bool posMuonQuality;
  bool posMuonHltIsoMu17;
  bool posMuonHltIsoMu24;
  bool posMuonHltIsoMu24Eta2p1;
  bool posMuonHltDoubleMu7;
  bool posMuonHltMu13Mu8Leg8;
  bool posMuonHltMu13Mu8Leg13;
  bool posMuonHltMu17Mu8Leg8;
  bool posMuonHltMu17Mu8Leg17;

  float negMuonPt;
  float negMuonEta;
  float negMuonPhi;
  float negMuonE;
  float negMuonTrackerIso03;
  float negMuonEcalIso03;
  float negMuonHcalIso03;
  float negMuonChargedHadronPfIso04;
  float negMuonChargedParticlePtPfIso04;
  float negMuonNeutralHadronEtPfIso04;
  float negMuonPhotonEtPfIso04;
  float negMuonPUPtPfIso04;
  bool negMuonQuality;
  bool negMuonHltIsoMu17;
  bool negMuonHltIsoMu24;
  bool negMuonHltIsoMu24Eta2p1;
  bool negMuonHltDoubleMu7;
  bool negMuonHltMu13Mu8Leg8;
  bool negMuonHltMu13Mu8Leg13;
  bool negMuonHltMu17Mu8Leg8;
  bool negMuonHltMu17Mu8Leg17;

  float genPosChargedLepPt;
  float genPosChargedLepEta;
  float genPosChargedLepPhi;
  float genPosChargedLepE;
  float genPosChargedLepFinalPt;
  float genPosChargedLepFinalEta;
  float genPosChargedLepFinalPhi;
  float genPosChargedLepFinalE;
  int genPosChargedLepPDG;

  float genNegChargedLepPt;
  float genNegChargedLepEta;
  float genNegChargedLepPhi;
  float genNegChargedLepE;
  float genNegChargedLepFinalPt;
  float genNegChargedLepFinalEta;
  float genNegChargedLepFinalPhi;
  float genNegChargedLepFinalE;
  int genNegChargedLepPDG;

  unsigned int nJets;
  float jetPt[20];
  float jetEta[20];
  float jetPhi[20];
  float jetEnergy[20];
  float jetBTag[20];

  unsigned int nTracks5;
  unsigned int nTracks10;
  unsigned int nTracks20;
  unsigned int nTracks30;
  unsigned int nTracks40;

  unsigned int nGlobalMuons;
  unsigned int nStandaloneMuons;
  unsigned int nPFMuons;

  float origPosMuonPt;
  float origPosMuonEta;
  float origPosMuonPhi;
  float origPosMuonE;
  float origPosMuonTrackerIso03;
  float origPosMuonEcalIso03;
  float origPosMuonHcalIso03;
  float origPosMuonChargedHadronPfIso04;
  float origPosMuonChargedParticlePtPfIso04;
  float origPosMuonNeutralHadronEtPfIso04;
  float origPosMuonPhotonEtPfIso04;
  float origPosMuonPUPtPfIso04;
  bool origPosMuonQuality;

  float origNegMuonPt;
  float origNegMuonEta;
  float origNegMuonPhi;
  float origNegMuonE;
  float origNegMuonTrackerIso03;
  float origNegMuonEcalIso03;
  float origNegMuonHcalIso03;
  float origNegMuonChargedHadronPfIso04;
  float origNegMuonChargedParticlePtPfIso04;
  float origNegMuonNeutralHadronEtPfIso04;
  float origNegMuonPhotonEtPfIso04;
  float origNegMuonPUPtPfIso04;
  bool origNegMuonQuality;

  float origGenPosChargedLepPt;
  float origGenPosChargedLepEta;
  float origGenPosChargedLepPhi;
  float origGenPosChargedLepE;
  float origGenPosChargedLepFinalPt;
  float origGenPosChargedLepFinalEta;
  float origGenPosChargedLepFinalPhi;
  float origGenPosChargedLepFinalE;
  int origGenPosChargedLepPDG;
  float origGenNegChargedLepPt;
  float origGenNegChargedLepEta;
  float origGenNegChargedLepPhi;
  float origGenNegChargedLepE;
  float origGenNegChargedLepFinalPt;
  float origGenNegChargedLepFinalEta;
  float origGenNegChargedLepFinalPhi;
  float origGenNegChargedLepFinalE;
  int origGenNegChargedLepPDG;
};

namespace
{

static bool isQualityMuon(const reco::Muon& muon, const edm::Handle<std::vector<reco::Vertex> >& vertexPtr)
{
  if(vertexPtr->empty()) return false;
  if(!muon.isGlobalMuon()) return false;
  if(!muon.isPFMuon()) return false;
  if(muon.globalTrack()->normalizedChi2() >= 10.0) return false;
  if(muon.globalTrack()->hitPattern().numberOfValidMuonHits() == 0) return false;
  if(muon.numberOfMatchedStations() <= 1) return false;
  if(fabs(muon.muonBestTrack()->dxy((*vertexPtr)[0].position())) >= 0.2) return false;
  if(fabs(muon.muonBestTrack()->dz((*vertexPtr)[0].position())) >= 0.5) return false;
  if(muon.innerTrack()->hitPattern().numberOfValidPixelHits() == 0) return false;
  if(muon.track()->hitPattern().trackerLayersWithMeasurement() <= 5) return false;
  return true;
}

bool isSuperior(const reco::Muon& candidate, const reco::Muon& to, const edm::Handle<std::vector<reco::Vertex> >& vertexPtr)
{
  const bool candidateQuality = isQualityMuon(candidate, vertexPtr);
  const bool toQuality = isQualityMuon(to, vertexPtr);

  if(candidateQuality != toQuality)
    return candidateQuality == true;
  else
    return candidate.pt() > to.pt();
}

unsigned int getTrigger(const HLTConfigProvider& hltConfig,
                        const std::string& triggerPattern)
{
  for(std::size_t i = 0; i < hltConfig.size(); ++i)
    if(hltConfig.triggerName(i).find(triggerPattern) == 0)
      return i;
  return NO_TRIGGER;
}

std::string getTriggerFilter(const HLTConfigProvider& hltConfig,
                             unsigned int triggerIndex,
                             const std::vector<std::string>& filterModules)
{
  const std::vector<std::string> saveTagsModules = hltConfig.saveTagsModules(triggerIndex);
  for(std::size_t j = 0; j < filterModules.size(); ++j)
    for(std::size_t k = 0; k < saveTagsModules.size(); ++k)
      if(filterModules[j] == saveTagsModules[k])
        return filterModules[j];

  std::cout << "Save tags modules:" << std::endl;
  for(std::size_t i = 0; i < saveTagsModules.size(); ++i)
    std::cout << "  " << saveTagsModules[i] << std::endl;

  std::cout << "Registered Filter Modules: " << std::endl;
  for(std::size_t i = 0; i < filterModules.size(); ++i)
    std::cout << "  " << filterModules[i] << std::endl;

  throw cms::Exception("MuMuNTupleProducer") << "Did not find filter for trigger " << hltConfig.triggerName(triggerIndex);
}

bool getTriggerMatch(const reco::Muon& muon,
                     const trigger::TriggerEvent& triggerEvent,
                     unsigned int triggerIndex,
                     const std::string& triggerProcess,
                     const std::string& filterName)
{
  if(triggerIndex != NO_TRIGGER)
  {
    const trigger::TriggerObjectCollection& triggerObjects = triggerEvent.getObjects();
    const unsigned int filterIndex = triggerEvent.filterIndex(edm::InputTag(filterName, "", triggerProcess));
    if(filterIndex < triggerEvent.sizeFilters())
    {
      const trigger::Keys& keys = triggerEvent.filterKeys(filterIndex);
      for(unsigned int i = 0; i < keys.size(); ++i)
      {
        if(ROOT::Math::VectorUtil::DeltaR(muon.p4(), triggerObjects[keys[i]].particle().p4()) < 0.3)
          return true;
      }
    }
  }

  return false;
}

std::pair<const reco::Muon*, const reco::Muon*> getMuons(const std::vector<reco::Muon>& muons, const edm::Handle<std::vector<reco::Vertex> >& vertexPtr, bool isSS)
{
   const reco::Muon* posMuon = NULL;
   const reco::Muon* negMuon = NULL;

   if(!isSS)
   {
      for(std::vector<reco::Muon>::const_iterator iter = muons.begin(); iter != muons.end(); ++iter)
      {
         if(iter->charge() > 0 && (!posMuon || isSuperior(*iter, *posMuon, vertexPtr))) posMuon = &*iter;
         if(iter->charge() < 0 && (!negMuon || isSuperior(*iter, *negMuon, vertexPtr))) negMuon = &*iter;
      }
   }
   else
   {
      const reco::Muon* posMuon1 = NULL;
      const reco::Muon* posMuon2 = NULL;
      const reco::Muon* negMuon1 = NULL;
      const reco::Muon* negMuon2 = NULL;
      for(std::vector<reco::Muon>::const_iterator iter = muons.begin(); iter != muons.end(); ++iter)
      {
         if(iter->charge() > 0 && (!posMuon1 || isSuperior(*iter, *posMuon1, vertexPtr))) { posMuon2 = posMuon1; posMuon1 = &*iter; }
         else if(iter->charge() > 0 && (!posMuon2 || isSuperior(*iter, *posMuon2, vertexPtr))) { posMuon2 = &*iter; }

         if(iter->charge() < 0 && (!negMuon1 || isSuperior(*iter, *negMuon1, vertexPtr))) { negMuon2 = negMuon1; negMuon1 = &*iter; }
         else if(iter->charge() < 0 && (!negMuon2 || isSuperior(*iter, *negMuon2, vertexPtr))) { negMuon2 = &*iter; }
      }

      // Decide whether to use the positive SS pair or the negative one
      if(!negMuon2) { posMuon = posMuon1; negMuon = posMuon2; }
      else if(!posMuon2) { posMuon = negMuon1; negMuon = negMuon2; }
      else if(isSuperior(*posMuon2, *negMuon2, vertexPtr)) { posMuon = posMuon1; negMuon = posMuon2; }
      else { posMuon = negMuon1; negMuon = negMuon2; }
   }

   return std::make_pair(posMuon, negMuon);
}

const reco::Candidate* getStableParticle(const reco::Candidate* part)
{
  const reco::Candidate* test = part;
  while(true)
  {
    const reco::Candidate* orig = test;
    for(unsigned int i = 0; i < test->numberOfDaughters(); ++i)
    {
      const reco::Candidate* cand = test->daughter(i);
      if(cand->pdgId() == test->pdgId())
      {
        test = cand;
	break;
      }
    }

    if(orig == test)
      return orig;
  }

  return NULL;
}

} // anonymous namespace

//
// constructors and destructor
//
MuMuNTupleProducer::MuMuNTupleProducer(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
   writeAllEvents = iConfig.getParameter<bool>("writeAllEvents");
   isEmbedded = iConfig.getParameter<bool>("isEmbedded");
   isRHEmbedded = iConfig.getParameter<bool>("isRHEmbedded");
   isGenEmbedded = iConfig.getParameter<bool>("isGenEmbedded");
   isData = iConfig.getParameter<bool>("isData");
   isSS = iConfig.getParameter<bool>("isSS");

   triggerSource = iConfig.getParameter<edm::InputTag>("TriggerSource");
   triggerEventSource = iConfig.getParameter<edm::InputTag>("TriggerEventSource");
   if(isEmbedded) origTriggerSource = iConfig.getParameter<edm::InputTag>("OrigTriggerSource");

   jetSource = iConfig.getParameter<edm::InputTag>("JetSource");
   btagSource = iConfig.getParameter<edm::InputTag>("BTagSource");
   trackSource = iConfig.getParameter<edm::InputTag>("TrackSource");

   muonIsoDepsChargedParticlesSource = iConfig.getParameter<edm::InputTag>("MuonIsoDepsChargedParticlesSource");
   muonIsoDepsChargedHadronsSource = iConfig.getParameter<edm::InputTag>("MuonIsoDepsChargedHadronsSource");
   muonIsoDepsNeutralHadronsSource = iConfig.getParameter<edm::InputTag>("MuonIsoDepsNeutralHadronsSource");
   muonIsoDepsPhotonsSource = iConfig.getParameter<edm::InputTag>("MuonIsoDepsPhotonsSource");
   muonIsoDepsPUSource = iConfig.getParameter<edm::InputTag>("MuonIsoDepsPUSource");

   muonSource = iConfig.getParameter<edm::InputTag>("MuonSource");
   if(isEmbedded) origMuonSource = iConfig.getParameter<edm::InputTag>("OrigMuonSource");
   rhoNeutralSource = iConfig.getParameter<edm::InputTag>("RhoNeutralSource");
   if(isEmbedded) origRhoNeutralSource = iConfig.getParameter<edm::InputTag>("OrigRhoNeutralSource");
   caloMetSource = iConfig.getParameter<edm::InputTag>("CaloMetSource");
   pfMetSource = iConfig.getParameter<edm::InputTag>("PfMetSource");
   pfType1CorrectedMetSource = iConfig.getParameter<edm::InputTag>("PfType1CorrectedMetSource");
   mvaPfMetSource = iConfig.getParameter<edm::InputTag>("MvaPfMetSource");
   if(isEmbedded) origCaloMetSource = iConfig.getParameter<edm::InputTag>("OrigCaloMetSource");
   if(isEmbedded) origPfMetSource = iConfig.getParameter<edm::InputTag>("OrigPfMetSource");
   vertexSource = iConfig.getParameter<edm::InputTag>("VertexSource");
   pileupSummaryInfoSource = iConfig.getParameter<edm::InputTag>("PileupSummaryInfoSource");
   if(isEmbedded) origVertexSource = iConfig.getParameter<edm::InputTag>("OrigVertexSource");

   if(isEmbedded || !isData) genParticlesSource = iConfig.getParameter<edm::InputTag>("GenParticlesSource");;
   if(isEmbedded && !isData) origGenParticlesSource = iConfig.getParameter<edm::InputTag>("OrigGenParticlesSource");;

   hltIsoMu17Index = hltIsoMu24Index = hltIsoMu24Eta2p1Index = hltDoubleMu7Index = hltMu13Mu8Index = hltMu17Mu8Index = NO_TRIGGER;
}


MuMuNTupleProducer::~MuMuNTupleProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

// ------------ method called for each event  ------------
void
MuMuNTupleProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   run = iEvent.id().run();
   lumi = iEvent.getLuminosityBlock().luminosityBlock();
   event = iEvent.id().event();
   //if(isData && !isCompatibleLS(run, lumi)) return;

   // Trigger
   edm::Handle<edm::TriggerResults> triggerResultsPtr, origTriggerResultsPtr;
   edm::Handle<trigger::TriggerEvent> triggerEventPtr;
   iEvent.getByLabel(triggerSource, triggerResultsPtr);
   if(isEmbedded && !isSS) iEvent.getByLabel(origTriggerSource, origTriggerResultsPtr);
   iEvent.getByLabel(triggerEventSource, triggerEventPtr);

   hltIsoMu17 = hltIsoMu24 = hltIsoMu24Eta2p1 = hltDoubleMu7 = hltMu13Mu8 = hltMu17Mu8 = false;
   if(hltIsoMu17Index != NO_TRIGGER) hltIsoMu17 = triggerResultsPtr->accept(hltIsoMu17Index);
   if(hltIsoMu24Index != NO_TRIGGER) hltIsoMu24 = triggerResultsPtr->accept(hltIsoMu24Index);
   if(hltIsoMu24Eta2p1Index != NO_TRIGGER) hltIsoMu24Eta2p1 = triggerResultsPtr->accept(hltIsoMu24Eta2p1Index);
   if(hltDoubleMu7Index != NO_TRIGGER) hltDoubleMu7 = triggerResultsPtr->accept(hltDoubleMu7Index);
   if(hltMu13Mu8Index != NO_TRIGGER) hltMu13Mu8 = triggerResultsPtr->accept(hltMu13Mu8Index);
   if(hltMu17Mu8Index != NO_TRIGGER) hltMu17Mu8 = triggerResultsPtr->accept(hltMu17Mu8Index);
  
   if(isEmbedded && !isGenEmbedded && !isSS)
   {
      origHltIsoMu17 = origHltDoubleMu7 = origHltMu13Mu8 = origHltMu17Mu8 = false;
      const edm::TriggerNames& origTriggerNames = iEvent.triggerNames(*origTriggerResultsPtr);
      for(unsigned int i = 0; i < origTriggerResultsPtr->size(); ++i)
      {
         if(origTriggerNames.triggerName(i).find("HLT_IsoMu17_v") != std::string::npos)
           origHltIsoMu17 = origTriggerResultsPtr->accept(i);
         if(origTriggerNames.triggerName(i).find("HLT_DoubleMu7_v") != std::string::npos)
           origHltDoubleMu7 = origTriggerResultsPtr->accept(i);
         if(origTriggerNames.triggerName(i).find("HLT_Mu13_Mu8_v") != std::string::npos)
           origHltMu13Mu8 = origTriggerResultsPtr->accept(i);
         if(origTriggerNames.triggerName(i).find("HLT_Mu17_Mu8_v") != std::string::npos)
           origHltMu17Mu8 = origTriggerResultsPtr->accept(i);
      }
   }

   // Rho
   edm::Handle<double> rhoPtr, origRhoPtr;
   iEvent.getByLabel(rhoNeutralSource, rhoPtr);
   if(isEmbedded && !isSS) iEvent.getByLabel(origRhoNeutralSource, origRhoPtr);

   rho = *rhoPtr;
   if(isEmbedded && !isSS) origRho = *origRhoPtr;
   else origRho = -1.0f;

   // MET
   edm::Handle<std::vector<reco::CaloMET> > caloMetPtr;
   edm::Handle<std::vector<reco::PFMET> > pfMetPtr, pfType1CorrectedMetPtr, mvaPfMetPtr;
   iEvent.getByLabel(caloMetSource, caloMetPtr);
   iEvent.getByLabel(pfMetSource, pfMetPtr);
   iEvent.getByLabel(pfType1CorrectedMetSource, pfType1CorrectedMetPtr);
   if(!mvaPfMetSource.label().empty()) iEvent.getByLabel(mvaPfMetSource, mvaPfMetPtr);

   caloMet = (*caloMetPtr)[0].pt();
   caloMetPhi = (*caloMetPtr)[0].phi();
   caloMetSigXX = (*caloMetPtr)[0].getSignificanceMatrix()(0,0);
   caloMetSigXY = (*caloMetPtr)[0].getSignificanceMatrix()(0,1);
   caloMetSigYX = (*caloMetPtr)[0].getSignificanceMatrix()(1,0);
   caloMetSigYY = (*caloMetPtr)[0].getSignificanceMatrix()(1,1);

   pfMet = (*pfMetPtr)[0].pt();
   pfMetPhi = (*pfMetPtr)[0].phi();
   pfMetSigXX = (*pfMetPtr)[0].getSignificanceMatrix()(0,0);
   pfMetSigXY = (*pfMetPtr)[0].getSignificanceMatrix()(0,1);
   pfMetSigYX = (*pfMetPtr)[0].getSignificanceMatrix()(1,0);
   pfMetSigYY = (*pfMetPtr)[0].getSignificanceMatrix()(1,1);

   pfType1CorrectedMet = (*pfType1CorrectedMetPtr)[0].pt();
   pfType1CorrectedMetPhi = (*pfType1CorrectedMetPtr)[0].phi();
   pfType1CorrectedMetSigXX = (*pfType1CorrectedMetPtr)[0].getSignificanceMatrix()(0,0);
   pfType1CorrectedMetSigXY = (*pfType1CorrectedMetPtr)[0].getSignificanceMatrix()(0,1);
   pfType1CorrectedMetSigYX = (*pfType1CorrectedMetPtr)[0].getSignificanceMatrix()(1,0);
   pfType1CorrectedMetSigYY = (*pfType1CorrectedMetPtr)[0].getSignificanceMatrix()(1,1);

   if(!mvaPfMetSource.label().empty())
   {
      mvaPfMet = (*mvaPfMetPtr)[0].pt();
      mvaPfMetPhi = (*mvaPfMetPtr)[0].phi();
      mvaPfMetSigXX = (*mvaPfMetPtr)[0].getSignificanceMatrix()(0,0);
      mvaPfMetSigXY = (*mvaPfMetPtr)[0].getSignificanceMatrix()(0,1);
      mvaPfMetSigYX = (*mvaPfMetPtr)[0].getSignificanceMatrix()(1,0);
      mvaPfMetSigYY = (*mvaPfMetPtr)[0].getSignificanceMatrix()(1,1);
   }

   if(isEmbedded && !isGenEmbedded)
   {
      edm::Handle<std::vector<reco::CaloMET> > origCaloMetPtr;
      edm::Handle<std::vector<reco::PFMET> > origPfMetPtr;
      iEvent.getByLabel(origCaloMetSource, origCaloMetPtr);
      iEvent.getByLabel(origPfMetSource, origPfMetPtr);

      origCaloMet = (*origCaloMetPtr)[0].pt();
      origCaloMetPhi = (*origCaloMetPtr)[0].phi();

      origPfMet = (*origPfMetPtr)[0].pt();
      origPfMetPhi = (*origPfMetPtr)[0].phi();
   }

   // Primary vertices
   edm::Handle<std::vector<reco::Vertex> > vertexPtr, origVertexPtr;
   iEvent.getByLabel(vertexSource, vertexPtr);
   if(isEmbedded && !isSS) iEvent.getByLabel(origVertexSource, origVertexPtr);
   nPV = vertexPtr->size();
   if(isEmbedded && !isSS) nOrigPV = origVertexPtr->size();
   else nOrigPV = 0;

   if(!isData)
   {
      edm::Handle<std::vector<PileupSummaryInfo> > puInfoPtr;
      iEvent.getByLabel(pileupSummaryInfoSource, puInfoPtr);

      nTrueInteractions = -1.f;
      for(std::vector<PileupSummaryInfo>::const_iterator iter = puInfoPtr->begin(); iter != puInfoPtr->end(); ++iter)
        if(iter->getBunchCrossing() == 0)
          nTrueInteractions = iter->getTrueNumInteractions();

      if(nTrueInteractions < 0.f)
        throw cms::Exception("MuMuNTupleProducer") << "No Pileup information for BX = 0 present";
   }
   else
   {
      nTrueInteractions = -1.f;
   }

   // Reconstructed Muons
   edm::Handle<std::vector<reco::Muon> > muonPtr, origMuonPtr;
   iEvent.getByLabel(muonSource, muonPtr);
   if(isEmbedded && !isSS) iEvent.getByLabel(origMuonSource, origMuonPtr);

   const std::pair<const reco::Muon*, const reco::Muon*> muons = getMuons(*muonPtr, vertexPtr, isSS);
   const reco::Muon* posMuon = muons.first;
   const reco::Muon* negMuon = muons.second;

   if(posMuon != NULL && negMuon != NULL)
   {
      reco::MuonRef posRef;
      reco::MuonRef negRef;

      for(unsigned int i = 0; i < muonPtr->size(); ++i)
      {
         if( &(*muonPtr)[i] == posMuon)
	   posRef = reco::MuonRef(muonPtr, i);
         if( &(*muonPtr)[i] == negMuon)
	   negRef = reco::MuonRef(muonPtr, i);
      }
      assert(posRef.isNonnull() && negRef.isNonnull());

      edm::Handle<edm::ValueMap<reco::IsoDeposit> > muonIsoDepsChargedParticlesPtr, muonIsoDepsChargedHadronsPtr, muonIsoDepsNeutralHadronsPtr, muonIsoDepsPhotonsPtr, muonIsoDepsPUPtr;
      iEvent.getByLabel(muonIsoDepsChargedParticlesSource, muonIsoDepsChargedParticlesPtr);
      iEvent.getByLabel(muonIsoDepsChargedHadronsSource, muonIsoDepsChargedHadronsPtr);
      iEvent.getByLabel(muonIsoDepsNeutralHadronsSource, muonIsoDepsNeutralHadronsPtr);
      iEvent.getByLabel(muonIsoDepsPhotonsSource, muonIsoDepsPhotonsPtr);
      iEvent.getByLabel(muonIsoDepsPUSource, muonIsoDepsPUPtr);

      reco::isodeposit::Direction posDir = reco::isodeposit::Direction(posRef->eta(), posRef->phi());
      reco::isodeposit::Direction negDir = reco::isodeposit::Direction(negRef->eta(), negRef->phi());
      reco::isodeposit::ConeVeto pos_pf_cone_veto_charged(posDir, 0.0001);
      reco::isodeposit::ConeVeto neg_pf_cone_veto_charged(negDir, 0.0001);
      reco::isodeposit::ThresholdVeto pf_threshold_veto_charged(0.0);
      reco::isodeposit::ConeVeto pos_pf_cone_veto(posDir, 0.01);
      reco::isodeposit::ConeVeto neg_pf_cone_veto(negDir, 0.01);
      reco::isodeposit::ThresholdVeto pf_threshold_veto(0.5);

      std::vector<reco::isodeposit::AbsVeto*> posVetosPFCharged;
      posVetosPFCharged.push_back(&pos_pf_cone_veto_charged);
      posVetosPFCharged.push_back(&pf_threshold_veto_charged);
      std::vector<reco::isodeposit::AbsVeto*> posVetosPF;
      posVetosPF.push_back(&pos_pf_cone_veto);
      posVetosPF.push_back(&pf_threshold_veto);

      std::vector<reco::isodeposit::AbsVeto*> negVetosPFCharged;
      negVetosPFCharged.push_back(&neg_pf_cone_veto_charged);
      negVetosPFCharged.push_back(&pf_threshold_veto_charged);
      std::vector<reco::isodeposit::AbsVeto*> negVetosPF;
      negVetosPF.push_back(&neg_pf_cone_veto);
      negVetosPF.push_back(&pf_threshold_veto);

      posMuonPt = posMuon->pt();
      posMuonEta = posMuon->eta();
      posMuonPhi = posMuon->phi();
      posMuonE = posMuon->energy();
      posMuonTrackerIso03 = posMuon->isolationR03().sumPt;
      posMuonEcalIso03 = posMuon->isolationR03().emEt;
      posMuonHcalIso03 = posMuon->isolationR03().hadEt;
      posMuonChargedHadronPfIso04 = (*muonIsoDepsChargedHadronsPtr)[posRef].depositWithin(0.4, posVetosPFCharged);
      posMuonChargedParticlePtPfIso04 = (*muonIsoDepsChargedParticlesPtr)[posRef].depositWithin(0.4, posVetosPFCharged);
      posMuonNeutralHadronEtPfIso04 = (*muonIsoDepsNeutralHadronsPtr)[posRef].depositWithin(0.4, posVetosPF);
      posMuonPhotonEtPfIso04 = (*muonIsoDepsPhotonsPtr)[posRef].depositWithin(0.4, posVetosPF);
      posMuonPUPtPfIso04 = (*muonIsoDepsPUPtr)[posRef].depositWithin(0.4, posVetosPF);
      posMuonQuality = isQualityMuon(*posMuon, vertexPtr);
      posMuonHltIsoMu17 = getTriggerMatch(*posMuon, *triggerEventPtr, hltIsoMu17Index, triggerEventSource.process(), hltIsoMu17Filter);
      posMuonHltIsoMu24 = getTriggerMatch(*posMuon, *triggerEventPtr, hltIsoMu24Index, triggerEventSource.process(), hltIsoMu24Filter);
      posMuonHltIsoMu24Eta2p1 = getTriggerMatch(*posMuon, *triggerEventPtr, hltIsoMu24Eta2p1Index, triggerEventSource.process(), hltIsoMu24Eta2p1Filter);
      posMuonHltDoubleMu7 = getTriggerMatch(*posMuon, *triggerEventPtr, hltDoubleMu7Index, triggerEventSource.process(), hltDoubleMu7Filter);
      posMuonHltMu13Mu8Leg8 = getTriggerMatch(*posMuon, *triggerEventPtr, hltMu13Mu8Index, triggerEventSource.process(), hltMu13Mu8Filter8);
      posMuonHltMu13Mu8Leg13 = getTriggerMatch(*posMuon, *triggerEventPtr, hltMu13Mu8Index, triggerEventSource.process(), hltMu13Mu8Filter13);
      posMuonHltMu17Mu8Leg8 = getTriggerMatch(*posMuon, *triggerEventPtr, hltMu17Mu8Index, triggerEventSource.process(), hltMu17Mu8Filter8);
      posMuonHltMu17Mu8Leg17 = getTriggerMatch(*posMuon, *triggerEventPtr, hltMu17Mu8Index, triggerEventSource.process(), hltMu17Mu8Filter17);

      negMuonPt = negMuon->pt();
      negMuonEta = negMuon->eta();
      negMuonPhi = negMuon->phi();
      negMuonE = negMuon->energy();
      negMuonTrackerIso03 = negMuon->isolationR03().sumPt;
      negMuonEcalIso03 = negMuon->isolationR03().emEt;
      negMuonHcalIso03 = negMuon->isolationR03().hadEt;
      negMuonChargedHadronPfIso04 = (*muonIsoDepsChargedHadronsPtr)[negRef].depositWithin(0.4, negVetosPFCharged);
      negMuonChargedParticlePtPfIso04 = (*muonIsoDepsChargedParticlesPtr)[negRef].depositWithin(0.4, negVetosPFCharged);
      negMuonNeutralHadronEtPfIso04 = (*muonIsoDepsNeutralHadronsPtr)[negRef].depositWithin(0.4, negVetosPF);
      negMuonPhotonEtPfIso04 = (*muonIsoDepsPhotonsPtr)[negRef].depositWithin(0.4, negVetosPF);
      negMuonPUPtPfIso04 = (*muonIsoDepsPUPtr)[negRef].depositWithin(0.4, negVetosPF);
      negMuonQuality = isQualityMuon(*negMuon, vertexPtr);
      negMuonHltIsoMu17 = getTriggerMatch(*negMuon, *triggerEventPtr, hltIsoMu17Index, triggerEventSource.process(), hltIsoMu17Filter);
      negMuonHltIsoMu24 = getTriggerMatch(*negMuon, *triggerEventPtr, hltIsoMu24Index, triggerEventSource.process(), hltIsoMu24Filter);
      negMuonHltIsoMu24Eta2p1 = getTriggerMatch(*negMuon, *triggerEventPtr, hltIsoMu24Eta2p1Index, triggerEventSource.process(), hltIsoMu24Eta2p1Filter);
      negMuonHltDoubleMu7 = getTriggerMatch(*negMuon, *triggerEventPtr, hltDoubleMu7Index, triggerEventSource.process(), hltDoubleMu7Filter);
      negMuonHltMu13Mu8Leg8 = getTriggerMatch(*negMuon, *triggerEventPtr, hltMu13Mu8Index, triggerEventSource.process(), hltMu13Mu8Filter8);
      negMuonHltMu13Mu8Leg13 = getTriggerMatch(*negMuon, *triggerEventPtr, hltMu13Mu8Index, triggerEventSource.process(), hltMu13Mu8Filter13);
      negMuonHltMu17Mu8Leg8 = getTriggerMatch(*negMuon, *triggerEventPtr, hltMu17Mu8Index, triggerEventSource.process(), hltMu17Mu8Filter8);
      negMuonHltMu17Mu8Leg17 = getTriggerMatch(*negMuon, *triggerEventPtr, hltMu17Mu8Index, triggerEventSource.process(), hltMu17Mu8Filter17);

      h_dilepmass->Fill( (posMuon->p4() + negMuon->p4()).M());

      // Jets
      edm::Handle<std::vector<reco::PFJet> > jetPtr;
      edm::Handle<reco::JetTagCollection> bTagPtr;
      iEvent.getByLabel(jetSource, jetPtr);
      iEvent.getByLabel(btagSource, bTagPtr);

      nJets = 0;
      for(unsigned int i = 0; i < jetPtr->size() && nJets < 20; ++i)
      {
         const reco::Jet* jet = &(*jetPtr)[i];
         if(ROOT::Math::VectorUtil::DeltaR(jet->p4(), negMuon->p4()) < 0.3 || ROOT::Math::VectorUtil::DeltaR(jet->p4(), posMuon->p4()) < 0.3)
           continue;
         if(jet->pt() < 20.0)
           continue;

         jetPt[nJets] = jet->pt();
         jetEta[nJets] = jet->eta();
         jetPhi[nJets] = jet->phi();
         jetEnergy[nJets] = jet->energy();

         unsigned int bTagIndex;
         for(bTagIndex = 0; bTagIndex < bTagPtr->size(); ++bTagIndex)
            if(ROOT::Math::VectorUtil::DeltaR(jet->p4(), (*bTagPtr)[bTagIndex].first->p4()) < 0.001)
               { jetBTag[nJets] = (*bTagPtr)[bTagIndex].second; break; }
         if(bTagIndex == bTagPtr->size()) throw cms::Exception("MuMuNTupleProducer") << "No BTag found for jet!";

         ++nJets;
      }

      // Tracks
      edm::Handle<reco::TrackCollection> trackPtr;
      iEvent.getByLabel(trackSource, trackPtr);

      nTracks5 = nTracks10 = nTracks20 = nTracks30 = nTracks40 = 0;
      for(unsigned int i = 0; i < trackPtr->size(); ++i)
      {
        const reco::Track& track = (*trackPtr)[i];
        if(track.pt() > 5) ++nTracks5;
        if(track.pt() > 10) ++nTracks10;
        if(track.pt() > 20) ++nTracks20;
        if(track.pt() > 30) ++nTracks30;
        if(track.pt() > 40) ++nTracks40;
      }

      // Muons
      nGlobalMuons = nStandaloneMuons = nPFMuons = 0;
      for(unsigned int i = 0; i < muonPtr->size(); ++i)
      {
        const reco::Muon& muon = (*muonPtr)[i];
        if(muon.isGlobalMuon()) ++nGlobalMuons;
        if(muon.isStandAloneMuon()) ++nStandaloneMuons;
        if(muon.isPFMuon()) ++nPFMuons;
      }
   }
   else
   {
      // no muon pair reconstructed. Set all variables to -1 to indicate this.
      posMuonPt = posMuonEta = posMuonPhi = posMuonE = posMuonTrackerIso03 = posMuonEcalIso03 = posMuonHcalIso03 = posMuonChargedHadronPfIso04 = posMuonChargedParticlePtPfIso04 = posMuonNeutralHadronEtPfIso04 = posMuonPhotonEtPfIso04 = posMuonPUPtPfIso04 = -1.0f;
      negMuonPt = negMuonEta = negMuonPhi = negMuonE = negMuonTrackerIso03 = negMuonEcalIso03 = negMuonHcalIso03 = negMuonChargedHadronPfIso04 = negMuonChargedParticlePtPfIso04 = negMuonNeutralHadronEtPfIso04 = negMuonPhotonEtPfIso04 = negMuonPUPtPfIso04 = -1.0f;
      posMuonHltIsoMu17 = posMuonHltIsoMu24 = posMuonHltIsoMu24Eta2p1 = posMuonHltDoubleMu7 = posMuonHltMu13Mu8Leg8 = posMuonHltMu13Mu8Leg13 = posMuonHltMu17Mu8Leg8 = posMuonHltMu17Mu8Leg17 = false;
      negMuonHltIsoMu17 = negMuonHltIsoMu24 = negMuonHltIsoMu24Eta2p1 = negMuonHltDoubleMu7 = negMuonHltMu13Mu8Leg8 = negMuonHltMu13Mu8Leg13 = negMuonHltMu17Mu8Leg8 = negMuonHltMu17Mu8Leg17 = false;

      posMuonQuality = false;
      negMuonQuality = false;
      nJets = 0;
      nTracks5 = nTracks10 = nTracks20 = nTracks30 = nTracks40 = 0;
      nGlobalMuons = nStandaloneMuons = nPFMuons = 0;
   }

   if(isEmbedded && !isGenEmbedded && !isSS)
   {
      const reco::Muon* posMuon = NULL;
      const reco::Muon* negMuon = NULL;

      // Use the actual muons that were used for the embedding, if they are available
      edm::Handle<std::vector<reco::CompositeCandidate> > compositePtr;
      iEvent.getByLabel("goldenZmumuCandidatesGe2IsoMuons", compositePtr);
      if(compositePtr.isValid() && !compositePtr->empty())
      {
         const reco::CompositeCandidate& comp = (*compositePtr)[0];
	 assert(comp.numberOfDaughters() == 2);
	 pat::MuonRef muon1 = comp.daughter(0)->masterClone().castTo<pat::MuonRef>();
	 pat::MuonRef muon2 = comp.daughter(1)->masterClone().castTo<pat::MuonRef>();
	 assert(muon1.isNonnull() && muon2.isNonnull());

	 if(muon1->charge() > 0)
	   { posMuon = &*muon1; negMuon = &*muon2; }
	 else
	   { negMuon = &*muon1; posMuon = &*muon2; }
      }
      else
      {
         const std::pair<const reco::Muon*, const reco::Muon*> muons = getMuons(*origMuonPtr, origVertexPtr, isSS);
         posMuon = muons.first;
         negMuon = muons.second;
      }

      // origMuons must always be available for embedded samples
      assert(posMuon != NULL && negMuon != NULL);

      origPosMuonPt = posMuon->pt();
      origPosMuonEta = posMuon->eta();
      origPosMuonPhi = posMuon->phi();
      origPosMuonE = posMuon->energy();
      origPosMuonTrackerIso03 = posMuon->isolationR03().sumPt;
      origPosMuonEcalIso03 = posMuon->isolationR03().emEt;
      origPosMuonHcalIso03 = posMuon->isolationR03().hadEt;
      origPosMuonChargedHadronPfIso04 = posMuon->pfIsolationR04().sumChargedHadronPt;
      origPosMuonChargedParticlePtPfIso04 = posMuon->pfIsolationR04().sumChargedParticlePt;
      origPosMuonNeutralHadronEtPfIso04 = posMuon->pfIsolationR04().sumNeutralHadronEt;
      origPosMuonPhotonEtPfIso04 = posMuon->pfIsolationR04().sumPhotonEt;
      origPosMuonPUPtPfIso04 = posMuon->pfIsolationR04().sumPUPt;
      origPosMuonQuality = isQualityMuon(*posMuon, origVertexPtr);

      origNegMuonPt = negMuon->pt();
      origNegMuonEta = negMuon->eta();
      origNegMuonPhi = negMuon->phi();
      origNegMuonE = negMuon->energy();
      origNegMuonTrackerIso03 = negMuon->isolationR03().sumPt;
      origNegMuonEcalIso03 = negMuon->isolationR03().emEt;
      origNegMuonHcalIso03 = negMuon->isolationR03().hadEt;
      origNegMuonChargedHadronPfIso04 = negMuon->pfIsolationR04().sumChargedHadronPt;
      origNegMuonChargedParticlePtPfIso04 = negMuon->pfIsolationR04().sumChargedParticlePt;
      origNegMuonNeutralHadronEtPfIso04 = negMuon->pfIsolationR04().sumNeutralHadronEt;
      origNegMuonPhotonEtPfIso04 = negMuon->pfIsolationR04().sumPhotonEt;
      origNegMuonPUPtPfIso04 = negMuon->pfIsolationR04().sumPUPt;
      origNegMuonQuality = isQualityMuon(*negMuon, origVertexPtr);

      h_origdilepmass->Fill( (posMuon->p4() + negMuon->p4()).M());
   }
   else
   {
      // Not an embedded sample. Set all variables to -1 to indicate this.
      origPosMuonPt = origPosMuonEta = origPosMuonPhi = origPosMuonE = origPosMuonTrackerIso03 = origPosMuonEcalIso03 = origPosMuonHcalIso03 = origPosMuonChargedHadronPfIso04 = origPosMuonChargedParticlePtPfIso04 = origPosMuonNeutralHadronEtPfIso04 = origPosMuonPhotonEtPfIso04 = origPosMuonPUPtPfIso04 = -1.0f;
      origNegMuonPt = origNegMuonEta = origNegMuonPhi = origNegMuonE = origNegMuonTrackerIso03 = origNegMuonEcalIso03 = origNegMuonHcalIso03 = origNegMuonChargedHadronPfIso04 = origNegMuonChargedParticlePtPfIso04 = origNegMuonNeutralHadronEtPfIso04 = origNegMuonPhotonEtPfIso04 = origNegMuonPUPtPfIso04 = -1.0f;

      origPosMuonQuality = false;
      origNegMuonQuality = false;
   }

   // Generated leptons
   genPosChargedLepPt = genPosChargedLepEta = genPosChargedLepPhi = genPosChargedLepE = -1.0f;
   genPosChargedLepFinalPt = genPosChargedLepFinalEta = genPosChargedLepFinalPhi = genPosChargedLepFinalE = -1.0f;
   genPosChargedLepPDG = 0;
   genNegChargedLepPt = genNegChargedLepEta = genNegChargedLepPhi = genNegChargedLepE = -1.0f;
   genNegChargedLepFinalPt = genNegChargedLepFinalEta = genNegChargedLepFinalPhi = genNegChargedLepFinalE = -1.0f;
   genNegChargedLepPDG = 0;
   if((isEmbedded && !isSS) || !isData)
   {
      edm::Handle<std::vector<reco::GenParticle> > genParticlesPtr;
      iEvent.getByLabel(genParticlesSource, genParticlesPtr);

      // Find exactly two opposite sign muons from Z
      for(unsigned int i = 0; i < genParticlesPtr->size(); ++i)
      {
         const reco::Candidate* cand = &(*genParticlesPtr)[i];
         if(cand->numberOfMothers() != 1 || (cand->mother(0)->pdgId() != 23 && cand->mother(0)->pdgId() != 25)) continue;
	 if(abs(cand->pdgId()) != 11 && abs(cand->pdgId()) != 13 && abs(cand->pdgId()) != 15)
            continue;

         const reco::Candidate* stable = getStableParticle(cand);

         if(cand->pdgId() < 0)
         {
            genPosChargedLepPt = cand->pt();
            genPosChargedLepEta = cand->eta();
            genPosChargedLepPhi = cand->phi();
            genPosChargedLepE = cand->energy();
            genPosChargedLepFinalPt = stable->pt();
            genPosChargedLepFinalEta = stable->eta();
            genPosChargedLepFinalPhi = stable->phi();
            genPosChargedLepFinalE = stable->energy();
            genPosChargedLepPDG = cand->pdgId();
         }
         else
         {
            genNegChargedLepPt = cand->pt();
            genNegChargedLepEta = cand->eta();
            genNegChargedLepPhi = cand->phi();
            genNegChargedLepE = cand->energy();
            genNegChargedLepFinalPt = stable->pt();
            genNegChargedLepFinalEta = stable->eta();
            genNegChargedLepFinalPhi = stable->phi();
            genNegChargedLepFinalE = stable->energy();
            genNegChargedLepPDG = cand->pdgId();
         }
      }
   }

   // Generated leptons in orig. event
   origGenPosChargedLepPt = origGenPosChargedLepEta = origGenPosChargedLepPhi = origGenPosChargedLepE = -1.0f;
   origGenPosChargedLepFinalPt = origGenPosChargedLepFinalEta = origGenPosChargedLepFinalPhi = origGenPosChargedLepFinalE = -1.0f;
   origGenPosChargedLepPDG = 0;
   origGenNegChargedLepPt = origGenNegChargedLepEta = origGenNegChargedLepPhi = origGenNegChargedLepE = -1.0f;
   origGenNegChargedLepFinalPt = origGenNegChargedLepFinalEta = origGenNegChargedLepFinalPhi = origGenNegChargedLepFinalE = -1.0f;
   origGenNegChargedLepPDG = 0;
   if(isEmbedded && !isSS && !isData)
   {
      edm::Handle<std::vector<reco::GenParticle> > origGenParticlesPtr;
      iEvent.getByLabel(origGenParticlesSource, origGenParticlesPtr);

      // Find exactly two opposite sign ME muons
      for(unsigned int i = 0; i < origGenParticlesPtr->size(); ++i)
      {
         const reco::Candidate* cand = &(*origGenParticlesPtr)[i];
         if(cand->numberOfMothers() != 1 || cand->mother(0)->pdgId() != 23) continue;
	 if(abs(cand->pdgId()) != 11 && abs(cand->pdgId()) != 13 && abs(cand->pdgId()) != 15)
            continue;

         const reco::Candidate* stable = getStableParticle(cand);

         if(cand->pdgId() < 0)
	 {
            origGenPosChargedLepPt = cand->pt();
            origGenPosChargedLepEta = cand->eta();
            origGenPosChargedLepPhi = cand->phi();
            origGenPosChargedLepE = cand->energy();
            origGenPosChargedLepFinalPt = stable->pt();
            origGenPosChargedLepFinalEta = stable->eta();
            origGenPosChargedLepFinalPhi = stable->phi();
            origGenPosChargedLepFinalE = stable->energy();
            origGenPosChargedLepPDG = cand->pdgId();
	 }
         else
	 {
            origGenNegChargedLepPt = cand->pt();
            origGenNegChargedLepEta = cand->eta();
            origGenNegChargedLepPhi = cand->phi();
            origGenNegChargedLepE = cand->energy();
            origGenNegChargedLepFinalPt = stable->pt();
            origGenNegChargedLepFinalEta = stable->eta();
            origGenNegChargedLepFinalPhi = stable->phi();
            origGenNegChargedLepFinalE = stable->energy();
            origGenNegChargedLepPDG = cand->pdgId();
	 }
      }
   }

   // Counting jets in madgraph DY+nJets
   if(!isData)
   {
      edm::Handle<std::vector<reco::GenParticle> > genParticlesPtr;
      if(isEmbedded)
         iEvent.getByLabel(origGenParticlesSource, genParticlesPtr);
      else
         iEvent.getByLabel(genParticlesSource, genParticlesPtr);

      nGenMEJets = 0;
      bool countJets = false;
      for(unsigned int i = 0; i < genParticlesPtr->size(); ++i)
      {
         const reco::Candidate* cand = &(*genParticlesPtr)[i];
	 if(cand->status() != 3) continue;

         const int pdg = abs(cand->pdgId());
	 if(pdg == 23) countJets = true;
	 else if(countJets)
	 {
	   if(pdg == 1 || pdg == 2 || pdg == 3 || pdg == 4 || pdg == 5 || pdg == 6 || pdg == 21)
             ++nGenMEJets;
	   if(pdg == 11 || pdg == 13 || pdg == 15)
	     countJets = false;
	 }
      }
   }

   // Get weight for embedded events
   minVisPtFilterWeight = 1.0f;
   zmumuEvtSelEffCorrWeight = 1.0f;
   muonRadiationFilter = true;
   muonRadiationFilter2Sel1 = true;
   muonRadiationFilter2Sel2 = true;
   muonRadiationFilter2Sel3 = true;

   if(isEmbedded)
   {
      edm::Handle<GenFilterInfo> hGenFilterInfo;
      iEvent.getByLabel(edm::InputTag("generator", "minVisPtFilter", "EmbeddedRECO"), hGenFilterInfo);
      minVisPtFilterWeight = hGenFilterInfo->filterEfficiency();

      if(isRHEmbedded)
      {
         edm::Handle<double> hZmumuEvtSelEffCorrWeight;
         iEvent.getByLabel(edm::InputTag("ZmumuEvtSelEffCorrWeightProducer", "weight", "EmbeddedRECO"), hZmumuEvtSelEffCorrWeight);

         // Some events don't have this... I don't yet understand the reason why
	 if(hZmumuEvtSelEffCorrWeight.isValid())
           zmumuEvtSelEffCorrWeight = *hZmumuEvtSelEffCorrWeight;

         edm::Handle<bool> hMuonRadiationFilter;
         iEvent.getByLabel(edm::InputTag("muonRadiationFilter", "", "EmbeddedRECO"), hMuonRadiationFilter);

         // Some events don't have this... I don't yet understand the reason why
	 if(hMuonRadiationFilter.isValid())
           muonRadiationFilter = *hMuonRadiationFilter;

         edm::Handle<bool> hMuonRadiationFilter2Sel1, hMuonRadiationFilter2Sel2, hMuonRadiationFilter2Sel3;
	 iEvent.getByLabel(edm::InputTag("muonRadiationFilter2ResultStripSel1"), hMuonRadiationFilter2Sel1);
	 iEvent.getByLabel(edm::InputTag("muonRadiationFilter2ResultStripSel2"), hMuonRadiationFilter2Sel2);
	 iEvent.getByLabel(edm::InputTag("muonRadiationFilter2ResultStripSel3"), hMuonRadiationFilter2Sel3);
	 muonRadiationFilter2Sel1 = *hMuonRadiationFilter2Sel1;
	 muonRadiationFilter2Sel2 = *hMuonRadiationFilter2Sel2;
	 muonRadiationFilter2Sel3 = *hMuonRadiationFilter2Sel3;
      }
   }

   // Fill primary vertex histograms for PileUp reweighting
   h_nPV->Fill(nPV);
   if(!isData) h_nTrueInteractions->Fill(nTrueInteractions);
   if(isEmbedded && !isSS) h_nOrigPV->Fill(nOrigPV);

   // Only fill tree for two good muons either in event or orig. event,
   // or if generator muons are within acceptance
   // if writeAllEvents is set, write all di-muon events
   bool fillTree = (writeAllEvents && abs(genPosChargedLepPDG) == 13 && abs(genNegChargedLepPDG) == 13);
   if(posMuonPt >= 0.0f || (isEmbedded && !isSS && origPosMuonPt >= 0.0f))
      fillTree = true;
   if(abs(genPosChargedLepPDG) == 13 && abs(genNegChargedLepPDG) == 13 && std::max(genPosChargedLepPt, genNegChargedLepPt) >= 17. && std::min(genPosChargedLepPt, genNegChargedLepPt) >= 8. && fabs(genPosChargedLepEta) <= 2.5 && fabs(genNegChargedLepEta) <= 2.5)
     fillTree = true;

   if(fillTree)
      tree->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
MuMuNTupleProducer::beginJob()
{
  edm::Service<TFileService> fs;

  tree = fs->make<TTree>("MuonTree", "MuonTree");

  tree->Branch("Run", &run, "Run/i");
  tree->Branch("Lumi", &lumi, "Lumi/i");
  tree->Branch("Event", &event, "Event/i");

  tree->Branch("MinVisPtFilterWeight", &minVisPtFilterWeight, "MinVisPtFilterWeight/F");
  tree->Branch("ZmumuEvtSelEffCorrWeight", &zmumuEvtSelEffCorrWeight, "ZmumuEvtSelEffCorrWeight/F");
  tree->Branch("MuonRadiationFilter", &muonRadiationFilter, "MuonRadiationFilter/O");
  tree->Branch("MuonRadiationFilter2Sel1", &muonRadiationFilter2Sel1, "MuonRadiationFilter2Sel1/O");
  tree->Branch("MuonRadiationFilter2Sel2", &muonRadiationFilter2Sel2, "MuonRadiationFilter2Sel2/O");
  tree->Branch("MuonRadiationFilter2Sel3", &muonRadiationFilter2Sel3, "MuonRadiationFilter2Sel3/O");

  tree->Branch("NPV", &nPV, "NPV/i");
  if(!isData) tree->Branch("NTrueInteractions", &nTrueInteractions, "NTrueInteractions/F");
  if(!isData) tree->Branch("NGenMEJets", &nGenMEJets, "NGenMEJets/I");

  tree->Branch("HltIsoMu17", &hltIsoMu17, "HltIsoMu17/O");
  tree->Branch("HltIsoMu24", &hltIsoMu24, "HltIsoMu24/O");
  tree->Branch("HltIsoMu24Eta2p1", &hltIsoMu24Eta2p1, "HltIsoMu24Eta2p1/O");
  tree->Branch("HltDoubleMu7", &hltDoubleMu7, "HltDoubleMu7/O");
  tree->Branch("HltMu13Mu8", &hltMu13Mu8, "HltMu13Mu8/O");
  tree->Branch("HltMu17Mu8", &hltMu17Mu8, "HltMu17Mu8/O");

  tree->Branch("Rho", &rho, "Rho/F");

  tree->Branch("CaloMet", &caloMet, "CaloMet/F");
  tree->Branch("CaloMetPhi", &caloMetPhi, "CaloMetPhi/F");
  tree->Branch("CaloMetSigXX", &caloMetSigXX, "CaloMetSigXX/F");
  tree->Branch("CaloMetSigXY", &caloMetSigXY, "CaloMetSigXY/F");
  tree->Branch("CaloMetSigYX", &caloMetSigYX, "CaloMetSigYX/F");
  tree->Branch("CaloMetSigYY", &caloMetSigYY, "CaloMetSigYY/F");

  tree->Branch("PfMet", &pfMet, "PfMet/F");
  tree->Branch("PfMetPhi", &pfMetPhi, "PfMetPhi/F");
  tree->Branch("PfMetSigXX", &pfMetSigXX, "PfMetSigXX/F");
  tree->Branch("PfMetSigXY", &pfMetSigXY, "PfMetSigXY/F");
  tree->Branch("PfMetSigYX", &pfMetSigYX, "PfMetSigYX/F");
  tree->Branch("PfMetSigYY", &pfMetSigYY, "PfMetSigYY/F");

  tree->Branch("PfType1CorrectedMet", &pfType1CorrectedMet, "PfType1CorrectedMet/F");
  tree->Branch("PfType1CorrectedMetPhi", &pfType1CorrectedMetPhi, "PfType1CorrectedMetPhi/F");
  tree->Branch("PfType1CorrectedMetSigXX", &pfType1CorrectedMetSigXX, "PfType1CorrectedMetSigXX/F");
  tree->Branch("PfType1CorrectedMetSigXY", &pfType1CorrectedMetSigXY, "PfType1CorrectedMetSigXY/F");
  tree->Branch("PfType1CorrectedMetSigYX", &pfType1CorrectedMetSigYX, "PfType1CorrectedMetSigYX/F");
  tree->Branch("PfType1CorrectedMetSigYY", &pfType1CorrectedMetSigYY, "PfType1CorrectedMetSigYY/F");

  if(!mvaPfMetSource.label().empty())
  {
    tree->Branch("MvaPfMet", &mvaPfMet, "MvaPfMet/F");
    tree->Branch("MvaPfMetPhi", &mvaPfMetPhi, "MvaPfMetPhi/F");
    tree->Branch("MvaPfMetSigXX", &mvaPfMetSigXX, "MvaPfMetSigXX/F");
    tree->Branch("MvaPfMetSigXY", &mvaPfMetSigXY, "MvaPfMetSigXY/F");
    tree->Branch("MvaPfMetSigYX", &mvaPfMetSigYX, "MvaPfMetSigYX/F");
    tree->Branch("MvaPfMetSigYY", &mvaPfMetSigYY, "MvaPfMetSigYY/F");
  }

  tree->Branch("PosMuonPt", &posMuonPt, "PosMuonPt/F");
  tree->Branch("PosMuonEta", &posMuonEta, "PosMuonEta/F");
  tree->Branch("PosMuonPhi", &posMuonPhi, "PosMuonPhi/F");
  tree->Branch("PosMuonE", &posMuonE, "PosMuonE/F");
  tree->Branch("PosMuonTrackerIso03", &posMuonTrackerIso03, "PosMuonTrackerIso03/F");
  tree->Branch("PosMuonEcalIso03", &posMuonEcalIso03, "PosMuonEcalIso03/F");
  tree->Branch("PosMuonHcalIso03", &posMuonHcalIso03, "PosMuonHcalIso03/F");
  tree->Branch("PosMuonChargedHadronPfIso04", &posMuonChargedHadronPfIso04, "PosMuonChargedHadronPfIso04/F");
  tree->Branch("PosMuonChargedParticlePtPfIso04", &posMuonChargedParticlePtPfIso04, "PosMuonChargedParticlePtPfIso04/F");
  tree->Branch("PosMuonNeutralHadronEtPfIso04", &posMuonNeutralHadronEtPfIso04, "PosMuonNeutralHadronEtPfIso04/F");
  tree->Branch("PosMuonPhotonEtPfIso04", &posMuonPhotonEtPfIso04, "PosMuonPhotonEtPfIso04/F");
  tree->Branch("PosMuonPUPtPfIso04", &posMuonPUPtPfIso04, "PosMuonPUPtPfIso04/F");
  tree->Branch("PosMuonQuality", &posMuonQuality, "PosMuonQuality/O");
  tree->Branch("PosMuonHltIsoMu17", &posMuonHltIsoMu17, "PosMuonHltIsoMu17/O");
  tree->Branch("PosMuonHltIsoMu24", &posMuonHltIsoMu24, "PosMuonHltIsoMu24/O");
  tree->Branch("PosMuonHltIsoMu24Eta2p1", &posMuonHltIsoMu24Eta2p1, "PosMuonHltIsoMu24Eta2p1/O");
  tree->Branch("PosMuonHltDoubleMu7", &posMuonHltDoubleMu7, "PosMuonHltDoubleMu7/O");
  tree->Branch("PosMuonHltMu13Mu8Leg8", &posMuonHltMu13Mu8Leg8, "PosMuonHltMu13Mu8Leg8/O");
  tree->Branch("PosMuonHltMu13Mu8Leg13", &posMuonHltMu13Mu8Leg13, "PosMuonHltMu13Mu8Leg13/O");
  tree->Branch("PosMuonHltMu17Mu8Leg8", &posMuonHltMu17Mu8Leg8, "PosMuonHltMu17Mu8Leg8/O");
  tree->Branch("PosMuonHltMu17Mu8Leg17", &posMuonHltMu17Mu8Leg17, "PosMuonHltMu17Mu8Leg17/O");

  tree->Branch("NegMuonPt", &negMuonPt, "NegMuonPt/F");
  tree->Branch("NegMuonEta", &negMuonEta, "NegMuonEta/F");
  tree->Branch("NegMuonPhi", &negMuonPhi, "NegMuonPhi/F");
  tree->Branch("NegMuonE", &negMuonE, "NegMuonE/F");
  tree->Branch("NegMuonTrackerIso03", &negMuonTrackerIso03, "NegMuonTrackerIso03/F");
  tree->Branch("NegMuonEcalIso03", &negMuonEcalIso03, "NegMuonEcalIso03/F");
  tree->Branch("NegMuonHcalIso03", &negMuonHcalIso03, "NegMuonHcalIso03/F");
  tree->Branch("NegMuonChargedHadronPfIso04", &negMuonChargedHadronPfIso04, "NegMuonChargedHadronPfIso04/F");
  tree->Branch("NegMuonChargedParticlePtPfIso04", &negMuonChargedParticlePtPfIso04, "NegMuonChargedParticlePtPfIso04/F");
  tree->Branch("NegMuonNeutralHadronEtPfIso04", &negMuonNeutralHadronEtPfIso04, "NegMuonNeutralHadronEtPfIso04/F");
  tree->Branch("NegMuonPhotonEtPfIso04", &negMuonPhotonEtPfIso04, "NegMuonPhotonEtPfIso04/F");
  tree->Branch("NegMuonPUPtPfIso04", &negMuonPUPtPfIso04, "NegMuonPUPtPfIso04/F");
  tree->Branch("NegMuonQuality", &negMuonQuality, "NegMuonQuality/O");
  tree->Branch("NegMuonHltIsoMu17", &negMuonHltIsoMu17, "NegMuonHltIsoMu17/O");
  tree->Branch("NegMuonHltIsoMu24", &negMuonHltIsoMu24, "NegMuonHltIsoMu24/O");
  tree->Branch("NegMuonHltIsoMu24Eta2p1", &negMuonHltIsoMu24Eta2p1, "NegMuonHltIsoMu24Eta2p1/O");
  tree->Branch("NegMuonHltDoubleMu7", &negMuonHltDoubleMu7, "NegMuonHltDoubleMu7/O");
  tree->Branch("NegMuonHltMu13Mu8Leg8", &negMuonHltMu13Mu8Leg8, "NegMuonHltMu13Mu8Leg8/O");
  tree->Branch("NegMuonHltMu13Mu8Leg13", &negMuonHltMu13Mu8Leg13, "NegMuonHltMu13Mu8Leg13/O");
  tree->Branch("NegMuonHltMu17Mu8Leg8", &negMuonHltMu17Mu8Leg8, "NegMuonHltMu17Mu8Leg8/O");
  tree->Branch("NegMuonHltMu17Mu8Leg17", &negMuonHltMu17Mu8Leg17, "NegMuonHltMu17Mu8Leg17/O");

  if(isEmbedded || !isData)
  {
    tree->Branch("GenPosChargedLepPt", &genPosChargedLepPt, "GenPosChargedLepPt/F");
    tree->Branch("GenPosChargedLepEta", &genPosChargedLepEta, "GenPosChargedLepEta/F");
    tree->Branch("GenPosChargedLepPhi", &genPosChargedLepPhi, "GenPosChargedLepPhi/F");
    tree->Branch("GenPosChargedLepE", &genPosChargedLepE, "GenPosChargedLepE/F");
    tree->Branch("GenPosChargedLepFinalPt", &genPosChargedLepFinalPt, "GenPosChargedLepFinalPt/F");
    tree->Branch("GenPosChargedLepFinalEta", &genPosChargedLepFinalEta, "GenPosChargedLepFinalEta/F");
    tree->Branch("GenPosChargedLepFinalPhi", &genPosChargedLepFinalPhi, "GenPosChargedLepFinalPhi/F");
    tree->Branch("GenPosChargedLepFinalE", &genPosChargedLepFinalE, "GenPosChargedLepFinalE/F");
    tree->Branch("GenPosChargedLepPDG", &genPosChargedLepPDG, "GenPosChargedLepPDG/I");

    tree->Branch("GenNegChargedLepPt", &genNegChargedLepPt, "GenNegChargedLepPt/F");
    tree->Branch("GenNegChargedLepEta", &genNegChargedLepEta, "GenNegChargedLepEta/F");
    tree->Branch("GenNegChargedLepPhi", &genNegChargedLepPhi, "GenNegChargedLepPhi/F");
    tree->Branch("GenNegChargedLepE", &genNegChargedLepE, "GenNegChargedLepE/F");
    tree->Branch("GenNegChargedLepFinalPt", &genNegChargedLepFinalPt, "GenNegChargedLepFinalPt/F");
    tree->Branch("GenNegChargedLepFinalEta", &genNegChargedLepFinalEta, "GenNegChargedLepFinalEta/F");
    tree->Branch("GenNegChargedLepFinalPhi", &genNegChargedLepFinalPhi, "GenNegChargedLepFinalPhi/F");
    tree->Branch("GenNegChargedLepFinalE", &genNegChargedLepFinalE, "GenNegChargedLepFinalE/F");
    tree->Branch("GenNegChargedLepPDG", &genNegChargedLepPDG, "GenNegChargedLepPDG/I");
  }

  tree->Branch("nJets", &nJets, "nJets/i");
  tree->Branch("jetPt", jetPt, "jetPt[nJets]/F");
  tree->Branch("jetEta", jetEta, "jetEta[nJets]/F");
  tree->Branch("jetPhi", jetPhi, "jetPhi[nJets]/F");
  tree->Branch("jetEnergy", jetEnergy, "jetEnergy[nJets]/F");
  tree->Branch("jetBTag", jetBTag, "jetBTag[nJets]/F");

  tree->Branch("nTracks5", &nTracks5, "nTracks5/i");
  tree->Branch("nTracks10", &nTracks10, "nTracks10/i");
  tree->Branch("nTracks20", &nTracks20, "nTracks20/i");
  tree->Branch("nTracks30", &nTracks30, "nTracks30/i");
  tree->Branch("nTracks40", &nTracks40, "nTracks40/i");

  tree->Branch("nGlobalMuons", &nGlobalMuons, "nGlobalMuons/i");
  tree->Branch("nStandaloneMuons", &nStandaloneMuons, "nStandaloneMuons/i");
  tree->Branch("nPFMuons", &nPFMuons, "nPFMuons/i");

  if(isEmbedded && !isGenEmbedded && !isSS)
  {
    tree->Branch("NOrigPV", &nOrigPV, "NOrigPV/i");

    tree->Branch("OrigHltIsoMu17", &origHltIsoMu17, "OrigHltIsoMu17/O");
    tree->Branch("OrigHltDoubleMu7", &origHltDoubleMu7, "OrigHltDoubleMu7/O");
    tree->Branch("OrigHltMu13Mu8", &origHltMu13Mu8, "OrigHltMu13Mu8/O");
    tree->Branch("OrigHltMu17Mu8", &origHltMu17Mu8, "OrigHltMu17Mu8/O");

    tree->Branch("OrigRho", &origRho, "OrigRho/F");

    tree->Branch("OrigCaloMet", &origCaloMet, "OrigCaloMet/F");
    tree->Branch("OrigCaloMetPhi", &origCaloMetPhi, "OrigCaloMetPhi/F");
    tree->Branch("OrigPfMet", &origPfMet, "OrigPfet/F");
    tree->Branch("OrigPfMetPhi", &origPfMetPhi, "OrigPfMetPhi/F");

    tree->Branch("OrigPosMuonPt", &origPosMuonPt, "OrigPosMuonPt/F");
    tree->Branch("OrigPosMuonEta", &origPosMuonEta, "OrigPosMuonEta/F");
    tree->Branch("OrigPosMuonPhi", &origPosMuonPhi, "OrigPosMuonPhi/F");
    tree->Branch("OrigPosMuonE", &origPosMuonE, "OrigPosMuonE/F");
    tree->Branch("OrigPosMuonTrackerIso03", &origPosMuonTrackerIso03, "OrigPosMuonTrackerIso03/F");
    tree->Branch("OrigPosMuonEcalIso03", &origPosMuonEcalIso03, "OrigPosMuonEcalIso03/F");
    tree->Branch("OrigPosMuonHcalIso03", &origPosMuonHcalIso03, "OrigPosMuonHcalIso03/F");
    tree->Branch("OrigPosMuonChargedHadronPfIso04", &origPosMuonChargedHadronPfIso04, "OrigPosMuonChargedHadronPfIso04/F");
    tree->Branch("OrigPosMuonChargedParticlePtPfIso04", &origPosMuonChargedParticlePtPfIso04, "OrigPosMuonChargedParticlePtPfIso04/F");
    tree->Branch("OrigPosMuonNeutralHadronEtPfIso04", &origPosMuonNeutralHadronEtPfIso04, "OrigPosMuonNeutralHadronEtPfIso04/F");
    tree->Branch("OrigPosMuonPhotonEtPfIso04", &origPosMuonPhotonEtPfIso04, "OrigPosMuonPhotonEtPfIso04/F");
    tree->Branch("OrigPosMuonPUPtPfIso04", &origPosMuonPUPtPfIso04, "OrigPosMuonPUPtPfIso04/F");
    tree->Branch("OrigPosMuonQuality", &origPosMuonQuality, "OrigPosMuonQuality/O");

    tree->Branch("OrigNegMuonPt", &origNegMuonPt, "OrigNegMuonPt/F");
    tree->Branch("OrigNegMuonEta", &origNegMuonEta, "OrigNegMuonEta/F");
    tree->Branch("OrigNegMuonPhi", &origNegMuonPhi, "OrigNegMuonPhi/F");
    tree->Branch("OrigNegMuonE", &origNegMuonE, "OrigNegMuonE/F");
    tree->Branch("OrigNegMuonTrackerIso03", &origNegMuonTrackerIso03, "OrigNegMuonTrackerIso03/F");
    tree->Branch("OrigNegMuonEcalIso03", &origNegMuonEcalIso03, "OrigNegMuonEcalIso03/F");
    tree->Branch("OrigNegMuonHcalIso03", &origNegMuonHcalIso03, "OrigNegMuonHcalIso03/F");
    tree->Branch("OrigNegMuonChargedHadronPfIso04", &origNegMuonChargedHadronPfIso04, "OrigNegMuonChargedHadronPfIso04/F");
    tree->Branch("OrigNegMuonChargedParticlePtPfIso04", &origNegMuonChargedParticlePtPfIso04, "OrigNegMuonChargedParticlePtPfIso04/F");
    tree->Branch("OrigNegMuonNeutralHadronEtPfIso04", &origNegMuonNeutralHadronEtPfIso04, "OrigNegMuonNeutralHadronEtPfIso04/F");
    tree->Branch("OrigNegMuonPhotonEtPfIso04", &origNegMuonPhotonEtPfIso04, "OrigNegMuonPhotonEtPfIso04/F");
    tree->Branch("OrigNegMuonPUPtPfIso04", &origNegMuonPUPtPfIso04, "OrigNegMuonPUPtPfIso04/F");
    tree->Branch("OrigNegMuonQuality", &origNegMuonQuality, "OrigNegMuonQuality/O");
  }

  if(isEmbedded && !isSS && !isData)
  {
    tree->Branch("OrigGenPosChargedLepPt", &origGenPosChargedLepPt, "OrigGenPosChargedLepPt/F");
    tree->Branch("OrigGenPosChargedLepEta", &origGenPosChargedLepEta, "OrigGenPosChargedLepEta/F");
    tree->Branch("OrigGenPosChargedLepPhi", &origGenPosChargedLepPhi, "OrigGenPosChargedLepPhi/F");
    tree->Branch("OrigGenPosChargedLepE", &origGenPosChargedLepE, "OrigGenPosChargedLepE/F");
    tree->Branch("OrigGenPosChargedLepFinalPt", &origGenPosChargedLepFinalPt, "OrigGenPosChargedLepFinalPt/F");
    tree->Branch("OrigGenPosChargedLepFinalEta", &origGenPosChargedLepFinalEta, "OrigGenPosChargedLepFinalEta/F");
    tree->Branch("OrigGenPosChargedLepFinalPhi", &origGenPosChargedLepFinalPhi, "OrigGenPosChargedLepFinalPhi/F");
    tree->Branch("OrigGenPosChargedLepFinalE", &origGenPosChargedLepFinalE, "OrigGenPosChargedLepFinalE/F");
    tree->Branch("OrigGenPosChargedLepPDG", &origGenPosChargedLepPDG, "OrigGenPosChargedLepPDG/I");

    tree->Branch("OrigGenNegChargedLepPt", &origGenNegChargedLepPt, "OrigGenNegChargedLepPt/F");
    tree->Branch("OrigGenNegChargedLepEta", &origGenNegChargedLepEta, "OrigGenNegChargedLepEta/F");
    tree->Branch("OrigGenNegChargedLepPhi", &origGenNegChargedLepPhi, "OrigGenNegChargedLepPhi/F");
    tree->Branch("OrigGenNegChargedLepE", &origGenNegChargedLepE, "OrigGenNegChargedLepE/F");
    tree->Branch("OrigGenNegChargedLepFinalPt", &origGenNegChargedLepFinalPt, "OrigGenNegChargedLepFinalPt/F");
    tree->Branch("OrigGenNegChargedLepFinalEta", &origGenNegChargedLepFinalEta, "OrigGenNegChargedLepFinalEta/F");
    tree->Branch("OrigGenNegChargedLepFinalPhi", &origGenNegChargedLepFinalPhi, "OrigGenNegChargedLepFinalPhi/F");
    tree->Branch("OrigGenNegChargedLepFinalE", &origGenNegChargedLepFinalE, "OrigGenNegChargedLepFinalE/F");
    tree->Branch("OrigGenNegChargedLepPDG", &origGenNegChargedLepPDG, "OrigGenNegChargedLepPDG/I");
  }

  h_nPV = fs->make<TH1F>("h_nPV", "Number of primary vertices", 50, -0.5, 49.5);
  if(isEmbedded && !isSS) h_nOrigPV = fs->make<TH1F>("h_nOrigPV", "Number of primary vertices in original event", 50, -0.5, 49.5);
  if(!isData) h_nTrueInteractions = fs->make<TH1F>("h_nTrueInteractions", "Number of true interactions", 100, -0.5, 99.5);
  h_dilepmass = fs->make<TH1F>("h_dilepmass", "Invariant di-muon mass", 150, 0, 150);
  h_origdilepmass = fs->make<TH1F>("h_origdilepmass", "Original Invariant di-muon mass", 150, 0, 150);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuMuNTupleProducer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
MuMuNTupleProducer::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
  bool changed = false;
  hltConfiguration.init(iRun, iSetup, triggerEventSource.process(), changed);

  hltIsoMu17Index = getTrigger(hltConfiguration, "HLT_IsoMu17_v");
  hltIsoMu24Index = getTrigger(hltConfiguration, "HLT_IsoMu24_v");
  hltIsoMu24Eta2p1Index = getTrigger(hltConfiguration, "HLT_IsoMu24_eta2p1_v");
  hltDoubleMu7Index = getTrigger(hltConfiguration, "HLT_DoubleMu7_v");
  hltMu13Mu8Index = getTrigger(hltConfiguration, "HLT_Mu13_Mu8_v");
  hltMu17Mu8Index = getTrigger(hltConfiguration, "HLT_Mu17_Mu8_v");

  const std::vector<std::string> HLT_ISOMU17_FILTERS = { "hltSingleMuIsoL3IsoFiltered17" };
  const std::vector<std::string> HLT_ISOMU24_FILTERS = { "hltSingleMuIsoL3IsoFiltered24", "hltSingleMuL2QualIsoL3IsoFiltered24", "hltL3crIsoL1sMu16L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15" };
  const std::vector<std::string> HLT_ISOMU24ETA2P1_FILTERS = { "hltL3IsoL1sMu14Eta2p1L1f0L2f14QL2IsoL3f24L3IsoFiltered", "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoFiltered10", "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15" };
  const std::vector<std::string> HLT_DOUBLEMU7_FILTERS = { "hltDiMuonL3PreFiltered7" };
  const std::vector<std::string> HLT_MU13MU8_FILTERS8 = { "hltDiMuonL3PreFiltered8", "hltDiMuonL3p5PreFiltered8", "hltL3pfL1DoubleMu10MuOpenL1f0L2pf0L3PreFiltered8", "hltL3pfL1DoubleMu10MuOpenOR3p5L1f0L2pf0L3PreFiltered8" };
  const std::vector<std::string> HLT_MU13MU8_FILTERS13 = { "hltSingleMu13L3Filtered13", "hltL3fL1DoubleMu10MuOpenL1f0L2f10L3Filtered13", "hltL3fL1DoubleMu10MuOpenOR3p5L1f0L2f10L3Filtered13" };
  const std::vector<std::string> HLT_MU17MU8_FILTERS8 = { "hltDiMuonL3PreFiltered8", "hltDiMuonL3p5PreFiltered8", "hltL3pfL1DoubleMu10MuOpenL1f0L2pf0L3PreFiltered8", "hltL3pfL1DoubleMu10MuOpenOR3p5L1f0L2pf0L3PreFiltered8" };
  const std::vector<std::string> HLT_MU17MU8_FILTERS17 = { "hltSingleMu13L3Filtered17", "hltL3fL1DoubleMu10MuOpenL1f0L2f10L3Filtered17", "hltL3fL1DoubleMu10MuOpenOR3p5L1f0L2f10L3Filtered17" };

  if(hltIsoMu17Index != NO_TRIGGER) hltIsoMu17Filter = getTriggerFilter(hltConfiguration, hltIsoMu17Index, HLT_ISOMU17_FILTERS);
  if(hltIsoMu24Index != NO_TRIGGER) hltIsoMu24Filter = getTriggerFilter(hltConfiguration, hltIsoMu24Index, HLT_ISOMU24_FILTERS);
  if(hltIsoMu24Eta2p1Index != NO_TRIGGER) hltIsoMu24Eta2p1Filter = getTriggerFilter(hltConfiguration, hltIsoMu24Eta2p1Index, HLT_ISOMU24ETA2P1_FILTERS);
  if(hltDoubleMu7Index != NO_TRIGGER) hltDoubleMu7Filter = getTriggerFilter(hltConfiguration, hltDoubleMu7Index, HLT_DOUBLEMU7_FILTERS);
  if(hltMu13Mu8Index != NO_TRIGGER) hltMu13Mu8Filter8 = getTriggerFilter(hltConfiguration, hltMu13Mu8Index, HLT_MU13MU8_FILTERS8);
  if(hltMu13Mu8Index != NO_TRIGGER) hltMu13Mu8Filter13 = getTriggerFilter(hltConfiguration, hltMu13Mu8Index, HLT_MU13MU8_FILTERS13);
  if(hltMu17Mu8Index != NO_TRIGGER) hltMu17Mu8Filter8 = getTriggerFilter(hltConfiguration, hltMu17Mu8Index, HLT_MU17MU8_FILTERS8);
  if(hltMu17Mu8Index != NO_TRIGGER) hltMu17Mu8Filter17 = getTriggerFilter(hltConfiguration, hltMu17Mu8Index, HLT_MU17MU8_FILTERS17);
}

// ------------ method called when ending the processing of a run  ------------
void 
MuMuNTupleProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
MuMuNTupleProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
MuMuNTupleProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuMuNTupleProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuMuNTupleProducer);
