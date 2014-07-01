// -*- C++ -*-
//
// Package:    SinglePionNTupleProducer
// Class:      SinglePionNTupleProducer
// 
/**\class SinglePionNTupleProducer SinglePionNTupleProducer.cc DesyHTauTau/SinglePionNTupleProducer/src/SinglePionNTupleProducer.cc

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

namespace
{

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

const reco::Candidate::LorentzVector getNoRadVector(const reco::Candidate* part)
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
    {
      if(orig->numberOfDaughters() == 0)
        return orig->p4();

      // Sum up all daugthers of orig that are not photons
      reco::Candidate::LorentzVector vec(0,0,0,0);
      for(unsigned int i = 0; i < orig->numberOfDaughters(); ++i)
      {
        const reco::Candidate* cand = orig->daughter(i);
	if(abs(cand->pdgId()) != 22) vec += cand->p4();
      }

      return vec;
    }
  }

  assert(false);
  return reco::Candidate::LorentzVector(0,0,0,0);
}

const reco::Candidate* getPionOrKaonDecay(const reco::Candidate* part)
{
  const reco::Candidate* cand1 = NULL;
  const reco::Candidate* cand2 = NULL;
  for(unsigned int i = 0; i < part->numberOfDaughters(); ++i)
  {
    const reco::Candidate* cand = part->daughter(i);

    if(abs(cand->pdgId()) != 22)
    {
      if(!cand1) cand1 = cand;
      else if(!cand2) cand2 = cand;
      else return NULL;
    }
  }

  assert(abs(cand1->pdgId()) == 16 || abs(cand2->pdgId()) == 16);

  const reco::Candidate* vis = cand1;
  if(abs(vis->pdgId()) == 16) vis = cand2;

  const int pion = 211;
  const int kaon = 321;
  if(abs(vis->pdgId()) != pion && abs(vis->pdgId()) != kaon)
    return NULL;

  return vis;
}

} // anonymous namespace

class SinglePionNTupleProducer : public edm::EDAnalyzer
{
public:
  explicit SinglePionNTupleProducer(const edm::ParameterSet&);
  ~SinglePionNTupleProducer();

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
  bool isEmbedded;

  edm::InputTag pileupSummaryInfoSource;
  edm::InputTag genParticlesSource;
  edm::InputTag origGenParticlesSource;

  TTree* tree;
  TH1F* h_nTrueInteractions;

  unsigned int run;
  unsigned int lumi;
  unsigned int event;

  float minVisPtFilterWeight;
  float zmumuEvtSelEffCorrWeight;
  float tauSpinWeight;

  bool muonRadiationFilter;
  bool muonRadiationFilter2Sel1;
  bool muonRadiationFilter2Sel2;
  bool muonRadiationFilter2Sel3;

  float nTrueInteractions;

  float genPosTauPt;
  float genPosTauEta;
  float genPosTauPhi;
  float genPosTauE;
  float genPosTauFinalPt;
  float genPosTauFinalEta;
  float genPosTauFinalPhi;
  float genPosTauFinalE;
  float genPosTauVisPt;
  float genPosTauVisEta;
  float genPosTauVisPhi;
  float genPosTauVisE;
  int genPosTauDecayPDG;

  float genNegTauPt;
  float genNegTauEta;
  float genNegTauPhi;
  float genNegTauE;
  float genNegTauFinalPt;
  float genNegTauFinalEta;
  float genNegTauFinalPhi;
  float genNegTauFinalE;
  float genNegTauVisPt;
  float genNegTauVisEta;
  float genNegTauVisPhi;
  float genNegTauVisE;
  int genNegTauDecayPDG;

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

SinglePionNTupleProducer::SinglePionNTupleProducer(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
   isEmbedded = iConfig.getParameter<bool>("isEmbedded");

   pileupSummaryInfoSource = iConfig.getParameter<edm::InputTag>("PileupSummaryInfoSource");
   genParticlesSource = iConfig.getParameter<edm::InputTag>("GenParticlesSource");
   if(isEmbedded) origGenParticlesSource = iConfig.getParameter<edm::InputTag>("OrigGenParticlesSource");;
}


SinglePionNTupleProducer::~SinglePionNTupleProducer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}

void
SinglePionNTupleProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   run = iEvent.id().run();
   lumi = iEvent.getLuminosityBlock().luminosityBlock();
   event = iEvent.id().event();

   edm::Handle<std::vector<PileupSummaryInfo> > puInfoPtr;
   iEvent.getByLabel(pileupSummaryInfoSource, puInfoPtr);

   nTrueInteractions = -1.f;
   if(puInfoPtr.isValid())
   {
     for(std::vector<PileupSummaryInfo>::const_iterator iter = puInfoPtr->begin(); iter != puInfoPtr->end(); ++iter)
       if(iter->getBunchCrossing() == 0)
         nTrueInteractions = iter->getTrueNumInteractions();

     if(nTrueInteractions < 0.f)
       throw cms::Exception("SinglePionNTupleProducer") << "No Pileup information for BX = 0 present";
   }

   h_nTrueInteractions->Fill(nTrueInteractions);

   // Generated leptons
   edm::Handle<std::vector<reco::GenParticle> > genParticlesPtr;
   iEvent.getByLabel(genParticlesSource, genParticlesPtr);

   // Find exactly two opposite sign taus from Z or H
   const reco::Candidate* posTau = NULL;
   const reco::Candidate* negTau = NULL;
   for(unsigned int i = 0; i < genParticlesPtr->size(); ++i)
   {
      const reco::Candidate* cand = &(*genParticlesPtr)[i];
      if(cand->numberOfMothers() != 1 || (cand->mother(0)->pdgId() != 23 && cand->mother(0)->pdgId() != 25)) continue;
      if(abs(cand->pdgId()) != 15)
            continue;
      assert(isEmbedded || cand->status() == 3);
      assert(cand->mother(0)->numberOfDaughters() == 3); // ???
      assert(cand->mother(0)->mother(0)->pdgId() != 23);

      if(cand->pdgId() < 0)
      {
         assert(posTau == NULL);
	 posTau = cand;
      }
      else
      {
         assert(negTau == NULL);
	 negTau = cand;
      }
   }

   // No tau decays
   if(posTau == NULL || negTau == NULL) return;

   // Check the tau decays
   const reco::Candidate* posStableTau = getStableParticle(posTau);
   const reco::Candidate* negStableTau = getStableParticle(negTau);

   const reco::Candidate::LorentzVector posNoRadTau = getNoRadVector(posStableTau);
   const reco::Candidate::LorentzVector negNoRadTau = getNoRadVector(negStableTau);

   const reco::Candidate* posPiKaTau = getPionOrKaonDecay(posStableTau);
   const reco::Candidate* negPiKaTau = getPionOrKaonDecay(negStableTau);

   if(posPiKaTau == NULL && negPiKaTau == NULL) return;

   genPosTauPt = posTau->pt();
   genPosTauEta = posTau->eta();
   genPosTauPhi = posTau->phi();
   genPosTauE = posTau->energy();
   genPosTauFinalPt = posNoRadTau.pt();
   genPosTauFinalEta = posNoRadTau.eta();
   genPosTauFinalPhi = posNoRadTau.phi();
   genPosTauFinalE = posNoRadTau.energy();
   if(posPiKaTau != NULL)
   {
     genPosTauVisPt = posPiKaTau->pt();
     genPosTauVisEta = posPiKaTau->eta();
     genPosTauVisPhi = posPiKaTau->phi();
     genPosTauVisE = posPiKaTau->energy();
     genPosTauDecayPDG = posPiKaTau->pdgId();
   }
   else
   {
     genPosTauVisPt = genPosTauVisEta = genPosTauVisPhi = genPosTauVisE = -1.f;
     genPosTauDecayPDG = 0;
   }

   genNegTauPt = negTau->pt();
   genNegTauEta = negTau->eta();
   genNegTauPhi = negTau->phi();
   genNegTauE = negTau->energy();
   genNegTauFinalPt = negNoRadTau.pt();
   genNegTauFinalEta = negNoRadTau.eta();
   genNegTauFinalPhi = negNoRadTau.phi();
   genNegTauFinalE = negNoRadTau.energy();
   if(negPiKaTau != NULL)
   {
     genNegTauVisPt = negPiKaTau->pt();
     genNegTauVisEta = negPiKaTau->eta();
     genNegTauVisPhi = negPiKaTau->phi();
     genNegTauVisE = negPiKaTau->energy();
     genNegTauDecayPDG = negPiKaTau->pdgId();
   }
   else
   {
     genNegTauVisPt = genNegTauVisEta = genNegTauVisPhi = genNegTauVisE = -1.f;
     genNegTauDecayPDG = 0.;
   }

   // Get weight for embedded events
   minVisPtFilterWeight = 1.0f;
   zmumuEvtSelEffCorrWeight = 1.0f;
   tauSpinWeight = 1.0f;

   muonRadiationFilter = true;
   muonRadiationFilter2Sel1 = true;
   muonRadiationFilter2Sel2 = true;
   muonRadiationFilter2Sel3 = true;

   if(isEmbedded)
   {
      edm::Handle<double> hTauSpinnerWT;
      iEvent.getByLabel(edm::InputTag("TauSpinnerReco", "TauSpinnerWT", "EmbeddedSPIN"), hTauSpinnerWT);
      tauSpinWeight = *hTauSpinnerWT;

      edm::Handle<GenFilterInfo> hGenFilterInfo;
      iEvent.getByLabel(edm::InputTag("generator", "minVisPtFilter", "EmbeddedRECO"), hGenFilterInfo);
      minVisPtFilterWeight = hGenFilterInfo->filterEfficiency();

      edm::Handle<double> hZmumuEvtSelEffCorrWeight;
      iEvent.getByLabel(edm::InputTag("ZmumuEvtSelEffCorrWeightProducer", "weight", "EmbeddedRECO"), hZmumuEvtSelEffCorrWeight);

      // Some events don't have this... I don't yet understand the reason why
      if(!hZmumuEvtSelEffCorrWeight.isValid()) return;
      zmumuEvtSelEffCorrWeight = *hZmumuEvtSelEffCorrWeight;

      edm::Handle<bool> hMuonRadiationFilter;
      iEvent.getByLabel(edm::InputTag("muonRadiationFilter", "", "EmbeddedRECO"), hMuonRadiationFilter);
      if(!hMuonRadiationFilter.isValid()) return; // skip event(?)
      muonRadiationFilter = *hMuonRadiationFilter;

      edm::Handle<bool> hMuonRadiationFilter2Sel1, hMuonRadiationFilter2Sel2, hMuonRadiationFilter2Sel3;
      iEvent.getByLabel(edm::InputTag("muonRadiationFilter2ResultStripSel1"), hMuonRadiationFilter2Sel1);
      iEvent.getByLabel(edm::InputTag("muonRadiationFilter2ResultStripSel2"), hMuonRadiationFilter2Sel2);
      iEvent.getByLabel(edm::InputTag("muonRadiationFilter2ResultStripSel3"), hMuonRadiationFilter2Sel3);
      muonRadiationFilter2Sel1 = *hMuonRadiationFilter2Sel1;
      muonRadiationFilter2Sel2 = *hMuonRadiationFilter2Sel2;
      muonRadiationFilter2Sel3 = *hMuonRadiationFilter2Sel3;
   }

   // Generated leptons in orig. event
   origGenPosChargedLepPt = origGenPosChargedLepEta = origGenPosChargedLepPhi = origGenPosChargedLepE = -1.0f;
   origGenPosChargedLepFinalPt = origGenPosChargedLepFinalEta = origGenPosChargedLepFinalPhi = origGenPosChargedLepFinalE = -1.0f;
   origGenPosChargedLepPDG = 0;
   origGenNegChargedLepPt = origGenNegChargedLepEta = origGenNegChargedLepPhi = origGenNegChargedLepE = -1.0f;
   origGenNegChargedLepFinalPt = origGenNegChargedLepFinalEta = origGenNegChargedLepFinalPhi = origGenNegChargedLepFinalE = -1.0f;
   origGenNegChargedLepPDG = 0;
   if(isEmbedded)
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
	 const reco::Candidate::LorentzVector noRad = getNoRadVector(stable);

         if(cand->pdgId() < 0)
	 {
            origGenPosChargedLepPt = cand->pt();
            origGenPosChargedLepEta = cand->eta();
            origGenPosChargedLepPhi = cand->phi();
            origGenPosChargedLepE = cand->energy();
            origGenPosChargedLepFinalPt = noRad.pt();
            origGenPosChargedLepFinalEta = noRad.eta();
            origGenPosChargedLepFinalPhi = noRad.phi();
            origGenPosChargedLepFinalE = noRad.energy();
            origGenPosChargedLepPDG = cand->pdgId();
	 }
         else
	 {
            origGenNegChargedLepPt = cand->pt();
            origGenNegChargedLepEta = cand->eta();
            origGenNegChargedLepPhi = cand->phi();
            origGenNegChargedLepE = cand->energy();
            origGenNegChargedLepFinalPt = noRad.pt();
            origGenNegChargedLepFinalEta = noRad.eta();
            origGenNegChargedLepFinalPhi = noRad.phi();
            origGenNegChargedLepFinalE = noRad.energy();
            origGenNegChargedLepPDG = cand->pdgId();
	 }
      }
   }

   tree->Fill();
}

void 
SinglePionNTupleProducer::beginJob()
{
  edm::Service<TFileService> fs;

  tree = fs->make<TTree>("TauTree", "TauTree");

  tree->Branch("Run", &run, "Run/i");
  tree->Branch("Lumi", &lumi, "Lumi/i");
  tree->Branch("Event", &event, "Event/i");

  tree->Branch("MinVisPtFilterWeight", &minVisPtFilterWeight, "MinVisPtFilterWeight/F");
  tree->Branch("ZmumuEvtSelEffCorrWeight", &zmumuEvtSelEffCorrWeight, "ZmumuEvtSelEffCorrWeight/F");
  tree->Branch("TauSpinWeight", &tauSpinWeight, "TauSpinWeight/F");

  tree->Branch("MuonRadiationFilter", &muonRadiationFilter, "MuonRadiationFilter/O");
  tree->Branch("MuonRadiationFilter2Sel1", &muonRadiationFilter2Sel1, "MuonRadiationFilter2Sel1/O");
  tree->Branch("MuonRadiationFilter2Sel2", &muonRadiationFilter2Sel2, "MuonRadiationFilter2Sel2/O");
  tree->Branch("MuonRadiationFilter2Sel3", &muonRadiationFilter2Sel3, "MuonRadiationFilter2Sel3/O");

  tree->Branch("NTrueInteractions", &nTrueInteractions, "NTrueInteractions/F");

  tree->Branch("GenPosTauPt", &genPosTauPt, "GenPosTauPt/F");
  tree->Branch("GenPosTauEta", &genPosTauEta, "GenPosTauEta/F");
  tree->Branch("GenPosTauPhi", &genPosTauPhi, "GenPosTauPhi/F");
  tree->Branch("GenPosTauE", &genPosTauE, "GenPosTauE/F");
  tree->Branch("GenPosTauFinalPt", &genPosTauFinalPt, "GenPosTauFinalPt/F");
  tree->Branch("GenPosTauFinalEta", &genPosTauFinalEta, "GenPosTauFinalEta/F");
  tree->Branch("GenPosTauFinalPhi", &genPosTauFinalPhi, "GenPosTauFinalPhi/F");
  tree->Branch("GenPosTauFinalE", &genPosTauFinalE, "GenPosTauFinalE/F");
  tree->Branch("GenPosTauVisPt", &genPosTauVisPt, "GenPosTauVisPt/F");
  tree->Branch("GenPosTauVisEta", &genPosTauVisEta, "GenPosTauVisEta/F");
  tree->Branch("GenPosTauVisPhi", &genPosTauVisPhi, "GenPosTauVisPhi/F");
  tree->Branch("GenPosTauVisE", &genPosTauVisE, "GenPosTauVisE/F");
  tree->Branch("GenPosTauDecayPDG", &genPosTauDecayPDG, "GenPosTauDecayPDG/I");

  tree->Branch("GenNegTauPt", &genNegTauPt, "GenNegTauPt/F");
  tree->Branch("GenNegTauEta", &genNegTauEta, "GenNegTauEta/F");
  tree->Branch("GenNegTauPhi", &genNegTauPhi, "GenNegTauPhi/F");
  tree->Branch("GenNegTauE", &genNegTauE, "GenNegTauE/F");
  tree->Branch("GenNegTauFinalPt", &genNegTauFinalPt, "GenNegTauFinalPt/F");
  tree->Branch("GenNegTauFinalEta", &genNegTauFinalEta, "GenNegTauFinalEta/F");
  tree->Branch("GenNegTauFinalPhi", &genNegTauFinalPhi, "GenNegTauFinalPhi/F");
  tree->Branch("GenNegTauFinalE", &genNegTauFinalE, "GenNegTauFinalE/F");
  tree->Branch("GenNegTauVisPt", &genNegTauVisPt, "GenNegTauVisPt/F");
  tree->Branch("GenNegTauVisEta", &genNegTauVisEta, "GenNegTauVisEta/F");
  tree->Branch("GenNegTauVisPhi", &genNegTauVisPhi, "GenNegTauVisPhi/F");
  tree->Branch("GenNegTauVisE", &genNegTauVisE, "GenNegTauVisE/F");
  tree->Branch("GenNegTauDecayPDG", &genNegTauDecayPDG, "GenNegTauDecayPDG/I");

  if(isEmbedded)
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

  h_nTrueInteractions = fs->make<TH1F>("h_nTrueInteractions", "Number of true interactions", 100, -0.5, 99.5);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SinglePionNTupleProducer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
SinglePionNTupleProducer::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
SinglePionNTupleProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
SinglePionNTupleProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
SinglePionNTupleProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SinglePionNTupleProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SinglePionNTupleProducer);
