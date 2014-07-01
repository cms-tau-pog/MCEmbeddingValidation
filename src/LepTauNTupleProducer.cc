// -*- C++ -*-
//
// Package:    LepTauNTupleProducer
// Class:      LepTauNTupleProducer
// 
/**\class LepTauNTupleProducer LepTauNTupleProducer.cc DesyHTauTau/LepTauNTupleProducer/src/LepTauNTupleProducer.cc

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

#include "DataFormats/RecoCandidate/interface/IsoDepositVetos.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "SimDataFormats/GeneratorProducts/interface/GenFilterInfo.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "TMath.h"
#include "TTree.h"
#include "TH1F.h"
#include "Math/VectorUtil.h"
#include "TVector3.h"

enum Channel {
  MUTAU,
  ETAU
};

class Electron: public reco::GsfElectron
{
public:
  Electron(const reco::GsfElectronRef& ref,
           const edm::Handle<reco::ConversionCollection>& hConversions,
           const edm::Handle<reco::BeamSpot>& hBeamSpot,
           const edm::ValueMap<float>& id,
           const edm::ValueMap<reco::IsoDeposit> isoDepsChargedParticles,
           const edm::ValueMap<reco::IsoDeposit> isoDepsChargedHadrons,
           const edm::ValueMap<reco::IsoDeposit> isoDepsNeutralHadrons,
           const edm::ValueMap<reco::IsoDeposit> isoDepsPhotons,
           const edm::ValueMap<reco::IsoDeposit> isoDepsPU):
    reco::GsfElectron(*ref),
    hasConversion_(ConversionTools::hasMatchedConversion(*ref, hConversions, hBeamSpot->position(), true, 2.0, 1e-6, 0)),
    mvaId_(id[ref])
  {
    reco::isodeposit::Direction dir = reco::isodeposit::Direction(ref->eta(), ref->phi());
    reco::isodeposit::ConeVeto pf_cone_veto_charged(dir, ref->isEB() ? 0.01 : 0.015);
    reco::isodeposit::ConeVeto pf_cone_veto_photons(dir, ref->isEB() ? 0.08 : 0.08);
    reco::isodeposit::ConeVeto pf_cone_veto_other(dir, 0.0);
    reco::isodeposit::ThresholdVeto pf_threshold_veto(0.0);

    std::vector<reco::isodeposit::AbsVeto*> vetosPFCharged, vetosPFPhotons, vetosPFOther;
    vetosPFCharged.push_back(&pf_cone_veto_charged);
    vetosPFCharged.push_back(&pf_threshold_veto);
    vetosPFPhotons.push_back(&pf_cone_veto_photons);
    vetosPFPhotons.push_back(&pf_threshold_veto);
    vetosPFOther.push_back(&pf_threshold_veto);

    chargedParticleIso_ = isoDepsChargedParticles[ref].depositWithin(0.4, vetosPFCharged);
    chargedHadronIso_ = isoDepsChargedHadrons[ref].depositWithin(0.4, vetosPFCharged);
    neutralHadronIso_ = isoDepsNeutralHadrons[ref].depositWithin(0.4, vetosPFOther);
    photonIso_ = isoDepsPhotons[ref].depositWithin(0.4, vetosPFPhotons);
    puIso_ = isoDepsPU[ref].depositWithin(0.4, vetosPFOther);
  }

  bool hasConversion() const { return hasConversion_; }
  float mvaId() const { return mvaId_; }

  float chargedParticleIso() const { return chargedParticleIso_; }
  float chargedHadronIso() const { return chargedHadronIso_; }
  float neutralHadronIso() const { return neutralHadronIso_; }
  float photonIso() const { return photonIso_; }
  float puIso() const { return puIso_; }

private:
  bool hasConversion_;
  float mvaId_;

  float chargedParticleIso_;
  float chargedHadronIso_;
  float neutralHadronIso_;
  float photonIso_;
  float puIso_;
};

class Muon: public reco::Muon
{
public:
  Muon(const reco::MuonRef& ref,
       const edm::ValueMap<reco::IsoDeposit> isoDepsChargedParticles,
       const edm::ValueMap<reco::IsoDeposit> isoDepsChargedHadrons,
       const edm::ValueMap<reco::IsoDeposit> isoDepsNeutralHadrons,
       const edm::ValueMap<reco::IsoDeposit> isoDepsPhotons,
       const edm::ValueMap<reco::IsoDeposit> isoDepsPU):
    reco::Muon(*ref)
  {
    reco::isodeposit::Direction dir = reco::isodeposit::Direction(ref->eta(), ref->phi());
    reco::isodeposit::ConeVeto pf_cone_veto_charged(dir, 0.0001);
    reco::isodeposit::ThresholdVeto pf_threshold_veto_charged(0.0);
    reco::isodeposit::ConeVeto pf_cone_veto(dir, 0.01);
    reco::isodeposit::ThresholdVeto pf_threshold_veto(0.5);

    std::vector<reco::isodeposit::AbsVeto*> vetosPFCharged;
    vetosPFCharged.push_back(&pf_cone_veto_charged);
    vetosPFCharged.push_back(&pf_threshold_veto_charged);
    std::vector<reco::isodeposit::AbsVeto*> vetosPF;
    vetosPF.push_back(&pf_cone_veto);
    vetosPF.push_back(&pf_threshold_veto);

    chargedParticleIso_ = isoDepsChargedParticles[ref].depositWithin(0.4, vetosPFCharged);
    chargedHadronIso_ = isoDepsChargedHadrons[ref].depositWithin(0.4, vetosPFCharged);
    neutralHadronIso_ = isoDepsNeutralHadrons[ref].depositWithin(0.4, vetosPF);
    photonIso_ = isoDepsPhotons[ref].depositWithin(0.4, vetosPF);
    puIso_ = isoDepsPU[ref].depositWithin(0.4, vetosPF);
  }

  float chargedParticleIso() const { return chargedParticleIso_; }
  float chargedHadronIso() const { return chargedHadronIso_; }
  float neutralHadronIso() const { return neutralHadronIso_; }
  float photonIso() const { return photonIso_; }
  float puIso() const { return puIso_; }
private:
  float chargedParticleIso_;
  float chargedHadronIso_;
  float neutralHadronIso_;
  float photonIso_;
  float puIso_;
};

class Tau: public reco::PFTau
{
public:
  Tau(const reco::PFTauRef& ref,
      const reco::PFTauDiscriminator& byDecayModeFinding,
      const reco::PFTauDiscriminator& byCombinedIsolationDBSumPtCorr3HitsRaw,
      const reco::PFTauDiscriminator& byLooseCombinedIsolationDBSumPtCorr3Hits,
      const reco::PFTauDiscriminator& byMediumCombinedIsolationDBSumPtCorr3Hits,
      const reco::PFTauDiscriminator& byTightCombinedIsolationDBSumPtCorr3Hits,
      const reco::PFTauDiscriminator& byLooseElectronRejection,
      const reco::PFTauDiscriminator& byLooseElectronRejectionMVA3,
      const reco::PFTauDiscriminator& byMediumElectronRejectionMVA3,
      const reco::PFTauDiscriminator& byTightElectronRejectionMVA3,
      const reco::PFTauDiscriminator& byLooseMuonRejection,
      const reco::PFTauDiscriminator& byTightMuonRejection):
    reco::PFTau(*ref),
    decayModeFinding_(byDecayModeFinding[ref] > 0.5),
    isolation3Hits_(byCombinedIsolationDBSumPtCorr3HitsRaw[ref]),
    looseIsolation3Hits_(byLooseCombinedIsolationDBSumPtCorr3Hits[ref] > 0.5),
    mediumIsolation3Hits_(byMediumCombinedIsolationDBSumPtCorr3Hits[ref] > 0.5),
    tightIsolation3Hits_(byTightCombinedIsolationDBSumPtCorr3Hits[ref] > 0.5),
    againstElectronLoose_(byLooseElectronRejection[ref] > 0.5),
    againstElectronLooseMVA3_(byLooseElectronRejectionMVA3[ref] > 0.5),
    againstElectronMediumMVA3_(byMediumElectronRejectionMVA3[ref] > 0.5),
    againstElectronTightMVA3_(byTightElectronRejectionMVA3[ref] > 0.5),
    againstMuonLoose_(byLooseMuonRejection[ref] > 0.5),
    againstMuonTight_(byTightMuonRejection[ref] > 0.5)
  {}

  bool decayModeFinding() const { return decayModeFinding_; }
  float isolation3Hits() const { return isolation3Hits_; }
  bool looseIsolation3Hits() const { return looseIsolation3Hits_; }
  bool mediumIsolation3Hits() const { return mediumIsolation3Hits_; }
  bool tightIsolation3Hits() const { return tightIsolation3Hits_; }
  bool againstElectronLoose() const { return againstElectronLoose_; }
  bool againstElectronLooseMVA3() const { return againstElectronLooseMVA3_; }
  bool againstElectronMediumMVA3() const { return againstElectronMediumMVA3_; }
  bool againstElectronTightMVA3() const { return againstElectronTightMVA3_; }
  bool againstMuonLoose() const { return againstMuonLoose_; }
  bool againstMuonTight() const { return againstMuonTight_; }

private:
  bool decayModeFinding_;
  float isolation3Hits_;
  bool looseIsolation3Hits_;
  bool mediumIsolation3Hits_;
  bool tightIsolation3Hits_;
  bool againstElectronLoose_;
  bool againstElectronLooseMVA3_;
  bool againstElectronMediumMVA3_;
  bool againstElectronTightMVA3_;
  bool againstMuonLoose_;
  bool againstMuonTight_;
};

namespace
{

static bool isOrigQualityMuon(const reco::Muon& muon, const edm::Handle<std::vector<reco::Vertex> >& vertexPtr)
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

bool isOrigSuperior(const reco::Muon& candidate, const reco::Muon& to, const edm::Handle<std::vector<reco::Vertex> >& vertexPtr)
{
  const bool candidateQuality = isOrigQualityMuon(candidate, vertexPtr);
  const bool toQuality = isOrigQualityMuon(to, vertexPtr);

  if(candidateQuality != toQuality)
    return candidateQuality == true;
  else
    return candidate.pt() > to.pt();
}

std::pair<const reco::Muon*, const reco::Muon*> getMuons(const std::vector<reco::Muon>& muons, const edm::Handle<std::vector<reco::Vertex> >& vertexPtr, bool isSS)
{
   const reco::Muon* posMuon = NULL;
   const reco::Muon* negMuon = NULL;

   if(!isSS)
   {
      for(std::vector<reco::Muon>::const_iterator iter = muons.begin(); iter != muons.end(); ++iter)
      {
         if(iter->charge() > 0 && (!posMuon || isOrigSuperior(*iter, *posMuon, vertexPtr))) posMuon = &*iter;
         if(iter->charge() < 0 && (!negMuon || isOrigSuperior(*iter, *negMuon, vertexPtr))) negMuon = &*iter;
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
         if(iter->charge() > 0 && (!posMuon1 || isOrigSuperior(*iter, *posMuon1, vertexPtr))) { posMuon2 = posMuon1; posMuon1 = &*iter; }
         else if(iter->charge() > 0 && (!posMuon2 || isOrigSuperior(*iter, *posMuon2, vertexPtr))) { posMuon2 = &*iter; }

         if(iter->charge() < 0 && (!negMuon1 || isOrigSuperior(*iter, *negMuon1, vertexPtr))) { negMuon2 = negMuon1; negMuon1 = &*iter; }
         else if(iter->charge() < 0 && (!negMuon2 || isOrigSuperior(*iter, *negMuon2, vertexPtr))) { negMuon2 = &*iter; }
      }

      // Decide whether to use the positive SS pair or the negative one
      if(!negMuon2) { posMuon = posMuon1; negMuon = posMuon2; }
      else if(!posMuon2) { posMuon = negMuon1; negMuon = negMuon2; }
      else if(isOrigSuperior(*posMuon2, *negMuon2, vertexPtr)) { posMuon = posMuon1; negMuon = posMuon2; }
      else { posMuon = negMuon1; negMuon = negMuon2; }
   }

   return std::make_pair(posMuon, negMuon);
}

static bool looseElectronId(const Electron& elec)
{
  const float mvaId = elec.mvaId();
  const float eta = fabs(elec.superCluster()->eta());

  if(elec.pt() < 20 && eta < 0.8) return mvaId > 0.925;
  else if(elec.pt() < 20 && eta >= 0.8 && eta < 1.479) return mvaId > 0.915;
  else if(elec.pt() < 20 && eta >= 1.479) return mvaId > 0.965;
  else if(elec.pt() >= 20 && eta < 0.8) return mvaId > 0.905;
  else if(elec.pt() >= 20 && eta >= 0.8 && eta < 1.479) return mvaId > 0.955;
  else if(elec.pt() >= 20 && eta >= 1.479) return mvaId > 0.975;
  else throw cms::Exception("looseElectronId") << "Invalid pt, eta: pt=" << elec.pt() << ", eta=" << eta;
}

static float leptonIso(const Electron& elec)
{
  const float charged = elec.chargedParticleIso();
  const float neutral = elec.neutralHadronIso();
  const float photon = elec.photonIso();
  const float pu = elec.puIso();

  return (charged + std::max(0., neutral + photon - 0.5 * pu)) / elec.pt();
}

static float leptonIso(const Muon& muon)
{
  const float charged = muon.chargedParticleIso();
  const float neutral = muon.neutralHadronIso();
  const float photon = muon.photonIso();
  const float pu = muon.puIso();

  return (charged + std::max(0., neutral + photon - 0.5 * pu)) / muon.pt();
}

static bool isQuality(const Electron& elec, const edm::Handle<std::vector<reco::Vertex> >& vertexPtr)
{
  if(vertexPtr->empty()) return false;
  if(fabs(elec.gsfTrack()->hitPattern().numberOfLostPixelHits()) > 0) return false;
  if(fabs(elec.gsfTrack()->dxy((*vertexPtr)[0].position())) >= 0.045) return false;
  if(fabs(elec.gsfTrack()->dz((*vertexPtr)[0].position())) >= 0.2) return false;
  if(elec.hasConversion()) return false;
  if(!looseElectronId(elec)) return false;
  return true;
}

static bool isQuality(const reco::Muon& muon, const edm::Handle<std::vector<reco::Vertex> >& vertexPtr)
{
  if(vertexPtr->empty()) return false;
  if(!muon.isGlobalMuon()) return false;
  if(!muon.isPFMuon()) return false;
  if(muon.globalTrack()->normalizedChi2() >= 10.0) return false;
  if(muon.globalTrack()->hitPattern().numberOfValidMuonHits() == 0) return false;
  if(muon.numberOfMatchedStations() <= 1) return false;
  if(fabs(muon.muonBestTrack()->dxy((*vertexPtr)[0].position())) >= 0.02) return false;
  if(fabs(muon.muonBestTrack()->dz((*vertexPtr)[0].position())) >= 0.2) return false;
  if(muon.innerTrack()->hitPattern().numberOfValidPixelHits() == 0) return false;
  if(muon.track()->hitPattern().trackerLayersWithMeasurement() <= 5) return false;
  return true;
}

static bool isQuality(const Tau& tau, const edm::Handle<std::vector<reco::Vertex> >& vertexPtr, Channel channel)
{
  if(vertexPtr->empty()) return false;
  if(fabs(tau.vertex().z() - (*vertexPtr)[0].z()) >= 0.2) return false;
  if(!tau.decayModeFinding()) return false;

  if(channel == ETAU)
  {
    if(!tau.againstElectronLooseMVA3()) return false;
    if(!tau.againstMuonLoose()) return false;
  }
  else
  {
    if(!tau.againstElectronLoose()) return false;
    if(!tau.againstMuonTight()) return false;
  }

  return true;
}

bool isSuperior(const Electron& candidate, const Electron& to, const edm::Handle<std::vector<reco::Vertex> >& vertexPtr)
{
  const bool candidateQuality = isQuality(candidate, vertexPtr);
  const bool toQuality = isQuality(to, vertexPtr);

  if(candidateQuality != toQuality)
    return candidateQuality == true;
  else if((candidate.pt() > 15.) != (to.pt() > 15.))
    return candidate.pt() > to.pt();
  else
    return leptonIso(candidate) < leptonIso(to);
}

bool isSuperior(const Muon& candidate, const Muon& to, const edm::Handle<std::vector<reco::Vertex> >& vertexPtr)
{
  const bool candidateQuality = isQuality(candidate, vertexPtr);
  const bool toQuality = isQuality(to, vertexPtr);

  if(candidateQuality != toQuality)
    return candidateQuality == true;
  else if((candidate.pt() > 15.) != (to.pt() > 15.))
    return candidate.pt() > to.pt();
  else
    return leptonIso(candidate) < leptonIso(to);
}

bool isSuperior(const Tau& candidate, const Tau& to, const edm::Handle<std::vector<reco::Vertex> >& vertexPtr, Channel channel)
{
  const bool candidateQuality = isQuality(candidate, vertexPtr, channel);
  const bool toQuality = isQuality(to, vertexPtr, channel);

  if(candidateQuality != toQuality)
    return candidateQuality == true;
  else if((candidate.pt() > 20.) != (to.pt() > 20.))
    return candidate.pt() > to.pt();
  else
    return candidate.isolation3Hits() < to.isolation3Hits();
}

template<typename TLepton>
std::pair<const TLepton*, const Tau*> getBetterPair(const std::pair<const TLepton*, const Tau*> first, const std::pair<const TLepton*, const Tau*> second, const edm::Handle<std::vector<reco::Vertex> >& vertexPtr, Channel channel)
{
   const bool qualityLep1 = first.first != NULL && isQuality(*first.first, vertexPtr);
   const bool qualityLep2 = second.first != NULL && isQuality(*second.first, vertexPtr);
   const bool qualityTau1 = first.second != NULL && isQuality(*first.second, vertexPtr, channel) && first.second->looseIsolation3Hits();
   const bool qualityTau2 = second.second != NULL && isQuality(*second.second, vertexPtr, channel) && second.second->looseIsolation3Hits();

   // First, decide by lepton quality
   if(qualityLep1 && !qualityLep2)
   {
      return first;
   }
   else if(!qualityLep1 && qualityLep2)
   {
      return second;
   }
   else if(!qualityLep1)
   {
      assert(!qualityLep2);
      if(!qualityTau1 && !qualityTau2) return std::pair<const TLepton*, const Tau*>(NULL, NULL);
      else if(qualityTau1 && !qualityTau2) return first;
      else if(!qualityTau1 && qualityTau2) return second;
      else if(first.second->isolation3Hits() < second.second->isolation3Hits()) return first;
      else return second;
   }
   else
   {
      assert(qualityLep2);
      if(!qualityTau1 && !qualityTau2)
        if(leptonIso(*first.first) < leptonIso(*second.first))
          return first;
        else
          return second;
      else if(qualityTau1 && !qualityTau2) return first;
      else if(!qualityTau1 && qualityTau2) return second;
      else if(first.second->isolation3Hits() < second.second->isolation3Hits()) return first;
      else return second;
   }
}

template<typename TLepton>
std::pair<const TLepton*, const Tau*> getLepTau(const std::vector<TLepton>& leptons, const std::vector<Tau>& taus, const edm::Handle<std::vector<reco::Vertex> >& vertexPtr, Channel channel, bool isSS)
{
   const TLepton* posLep = NULL;
   const TLepton* negLep = NULL;
   const Tau* posTau = NULL;
   const Tau* negTau = NULL;

   for(typename std::vector<TLepton>::const_iterator iter = leptons.begin(); iter != leptons.end(); ++iter)
   {
      if(iter->charge() > 0 && (!posLep || isSuperior(*iter, *posLep, vertexPtr)))
         posLep = &*iter;
      if(iter->charge() < 0 && (!negLep || isSuperior(*iter, *negLep, vertexPtr)))
         negLep = &*iter;
   }

   for(std::vector<Tau>::const_iterator iter = taus.begin(); iter != taus.end(); ++iter)
   {
      if(iter->charge() > 0 && (!posTau || isSuperior(*iter, *posTau, vertexPtr, channel)))
         posTau = &*iter;
      if(iter->charge() < 0 && (!negTau || isSuperior(*iter, *negTau, vertexPtr, channel)))
         negTau = &*iter;
   }

   if(!isSS)
   {
      return getBetterPair(std::make_pair(posLep, negTau), std::make_pair(negLep, posTau), vertexPtr, channel);
   }
   else
   {
      return getBetterPair(std::make_pair(posLep, posTau), std::make_pair(negLep, negTau), vertexPtr, channel);
   }
}

const reco::Candidate* getVisTauDecay(const reco::Candidate& genTau)
{
   const reco::Candidate* curTau = &genTau;
   while(curTau != NULL && abs(curTau->pdgId()) == 15)
   {
      const reco::Candidate* nextTau = NULL;
      for(unsigned int i = 0; i < curTau->numberOfDaughters(); ++i)
      {
         const reco::Candidate* daughter = curTau->daughter(i);
         if(abs(daughter->pdgId()) == 12 || abs(daughter->pdgId()) == 14 || abs(daughter->pdgId()) == 16) continue;

         if(abs(daughter->pdgId()) == 22) continue; // skip gamma radiation

         assert(nextTau == NULL);
         nextTau = daughter;
      }

      curTau = nextTau;
   }

   return curTau;
}

const reco::Candidate::LorentzVector getStableParticle(const reco::Candidate* part)
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

} // anonymous namespace

class LepTauNTupleProducer : public edm::EDAnalyzer
{
public:
  explicit LepTauNTupleProducer(const edm::ParameterSet&);
  ~LepTauNTupleProducer();

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
  const bool isData;
  const bool isEmbedded;
  const bool isRHEmbedded;
  const Channel channel;

  edm::InputTag triggerSource;
  edm::InputTag triggerEventSource;
  edm::InputTag origTriggerSource;

  edm::InputTag jetSource;
  edm::InputTag btagSource;
  edm::InputTag trackSource;
  edm::InputTag pfNoPileupSource;
  edm::InputTag pfPileupSource;

  edm::InputTag electronSource;
  edm::InputTag muonSource;
  edm::InputTag tauSource;

  edm::InputTag origMuonSource;

  edm::InputTag conversionsSource;
  edm::InputTag electronIDSource;
  edm::InputTag electronIsoDepsChargedParticlesSource;
  edm::InputTag electronIsoDepsChargedHadronsSource;
  edm::InputTag electronIsoDepsNeutralHadronsSource;
  edm::InputTag electronIsoDepsPhotonsSource;
  edm::InputTag electronIsoDepsPUSource;

  edm::InputTag muonIsoDepsChargedParticlesSource;
  edm::InputTag muonIsoDepsChargedHadronsSource;
  edm::InputTag muonIsoDepsNeutralHadronsSource;
  edm::InputTag muonIsoDepsPhotonsSource;
  edm::InputTag muonIsoDepsPUSource;

  edm::InputTag tauByDecayModeFindingSource;
  edm::InputTag tauByCombinedIsolationDBSumPtCorr3HitsRawSource;
  edm::InputTag tauByLooseCombinedIsolationDBSumPtCorr3HitsSource;
  edm::InputTag tauByMediumCombinedIsolationDBSumPtCorr3HitsSource;
  edm::InputTag tauByTightCombinedIsolationDBSumPtCorr3HitsSource;
  edm::InputTag tauByLooseElectronRejectionSource;
  edm::InputTag tauByLooseElectronRejectionMVA3Source;
  edm::InputTag tauByMediumElectronRejectionMVA3Source;
  edm::InputTag tauByTightElectronRejectionMVA3Source;
  edm::InputTag tauByLooseMuonRejectionSource;
  edm::InputTag tauByTightMuonRejectionSource;

  edm::InputTag caloMetSource;
  edm::InputTag pfMetSource;
  edm::InputTag pfType1CorrectedMetSource;
  edm::InputTag mvaPfMetSource;

  edm::InputTag beamSpotSource;
  edm::InputTag vertexSource;
  edm::InputTag pileupSummaryInfoSource;
  edm::InputTag origVertexSource;

  edm::InputTag genParticlesSource;
  edm::InputTag origGenParticlesSource;

  HLTConfigProvider hltConfiguration;

  TTree* tree;
  TH1F* h_nPV;
  TH1F* h_nTrueInteractions;

  unsigned int run;
  unsigned int lumi;
  unsigned int event;

  bool origHltIsoMu17;
  bool origHltDoubleMu7;
  bool origHltMu13Mu8;
  bool origHltMu17Mu8;

  float minVisPtFilterWeight;
  float tauSpinWeight;
  float zmumuEvtSelEffCorrWeight;
  float muonRadiationWeight;
  bool muonRadiationFilter;
  bool muonRadiationFilter2Sel1;
  bool muonRadiationFilter2Sel2;
  bool muonRadiationFilter2Sel3;

  unsigned int nPV;
  float nTrueInteractions;

  float origPosMuonPt;
  float origPosMuonEta;
  float origPosMuonPhi;
  float origPosMuonE;
  float origNegMuonPt;
  float origNegMuonEta;
  float origNegMuonPhi;
  float origNegMuonE;

  float beamSpotX;
  float beamSpotY;
  float beamSpotZ;

  float vertexX;
  float vertexY;
  float vertexZ;
  float vertexChi2;
  float vertexNdof;
  bool vertexFake;
  bool vertexValid;

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

  float electronPt;
  float electronEta;
  float electronPhi;
  float electronE;
  int electronCharge;
  float electronChargedParticlePfIso04;
  float electronChargedHadronPfIso04;
  float electronNeutralHadronEtPfIso04;
  float electronPhotonEtPfIso04;
  float electronPUPtPfIso04;
  bool electronQuality;
  bool electronHasConversion;
  float electronMva;

  float muonPt;
  float muonEta;
  float muonPhi;
  float muonE;
  int muonCharge;
  float muonChargedParticlePfIso04;
  float muonChargedHadronPfIso04;
  float muonNeutralHadronEtPfIso04;
  float muonPhotonEtPfIso04;
  float muonPUPtPfIso04;
  bool muonQuality;

  float tauPt;
  float tauEta;
  float tauPhi;
  float tauE;
  int tauCharge;
  bool tauQuality;
  bool tauDecayModeFinding;
  float tauIsolation3Hits;
  bool tauLooseIsolation3Hits;
  bool tauMediumIsolation3Hits;
  bool tauTightIsolation3Hits;
  bool tauAgainstMuonLoose;
  bool tauAgainstMuonTight;
  bool tauAgainstElectronLoose;
  bool tauAgainstElectronLooseMVA3;
  bool tauAgainstElectronMediumMVA3;
  bool tauAgainstElectronTightMVA3;

  float electronGenPt;
  float electronGenEta;
  float electronGenPhi;
  float electronGenE;
  float electronGenFinalPt;
  float electronGenFinalEta;
  float electronGenFinalPhi;
  float electronGenFinalE;
  float electronGenPtVis;
  float electronGenEtaVis;
  float electronGenPhiVis;
  float electronGenEVis;
  int electronGenCharge;
  float muonGenPt;
  float muonGenEta;
  float muonGenPhi;
  float muonGenE;
  float muonGenFinalPt;
  float muonGenFinalEta;
  float muonGenFinalPhi;
  float muonGenFinalE;
  float muonGenPtVis;
  float muonGenEtaVis;
  float muonGenPhiVis;
  float muonGenEVis;
  int muonGenCharge;
  float tauGenPt;
  float tauGenEta;
  float tauGenPhi;
  float tauGenE;
  float tauGenFinalPt;
  float tauGenFinalEta;
  float tauGenFinalPhi;
  float tauGenFinalE;
  float tauGenPtVis;
  float tauGenEtaVis;
  float tauGenPhiVis;
  float tauGenEVis;
  int tauGenCharge;
  int tauGenDecayPDG;

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

  unsigned int nJets;
  float jetPt[20];
  float jetEta[20];
  float jetPhi[20];
  float jetEnergy[20];
  float jetBTag[20];
};

static Channel get_channel_from_string(const std::string& str)
{
  if(str == "mutau") return MUTAU;
  else if(str == "etau") return ETAU;
  else throw cms::Exception("LepTauNTupleProducer") << "Invalid channel: " << str;
}

LepTauNTupleProducer::LepTauNTupleProducer(const edm::ParameterSet& iConfig):
   isData(iConfig.getParameter<bool>("isData")),
   isEmbedded(iConfig.getParameter<bool>("isEmbedded")),
   isRHEmbedded(isEmbedded && iConfig.getParameter<bool>("isRHEmbedded")),
   channel(get_channel_from_string(iConfig.getParameter<std::string>("channel")))
{
   //now do what ever initialization is needed
   triggerSource = iConfig.getParameter<edm::InputTag>("TriggerSource");
   triggerEventSource = iConfig.getParameter<edm::InputTag>("TriggerEventSource");
   origTriggerSource = iConfig.getParameter<edm::InputTag>("OrigTriggerSource");

   jetSource = iConfig.getParameter<edm::InputTag>("JetSource");
   btagSource = iConfig.getParameter<edm::InputTag>("BTagSource");
   trackSource = iConfig.getParameter<edm::InputTag>("TrackSource");
   pfNoPileupSource = iConfig.getParameter<edm::InputTag>("PfNoPileupSource");
   pfPileupSource = iConfig.getParameter<edm::InputTag>("PfPileupSource");

   electronSource = iConfig.getParameter<edm::InputTag>("ElectronSource");
   muonSource = iConfig.getParameter<edm::InputTag>("MuonSource");
   tauSource = iConfig.getParameter<edm::InputTag>("TauSource");

   if(isEmbedded) origMuonSource = iConfig.getParameter<edm::InputTag>("OrigMuonSource");

   conversionsSource = iConfig.getParameter<edm::InputTag>("ConversionsSource");
   electronIDSource = iConfig.getParameter<edm::InputTag>("ElectronIDSource");
   electronIsoDepsChargedParticlesSource = iConfig.getParameter<edm::InputTag>("ElectronIsoDepsChargedParticlesSource");
   electronIsoDepsChargedHadronsSource = iConfig.getParameter<edm::InputTag>("ElectronIsoDepsChargedHadronsSource");
   electronIsoDepsNeutralHadronsSource = iConfig.getParameter<edm::InputTag>("ElectronIsoDepsNeutralHadronsSource");
   electronIsoDepsPhotonsSource = iConfig.getParameter<edm::InputTag>("ElectronIsoDepsPhotonsSource");
   electronIsoDepsPUSource = iConfig.getParameter<edm::InputTag>("ElectronIsoDepsPUSource");

   muonIsoDepsChargedParticlesSource = iConfig.getParameter<edm::InputTag>("MuonIsoDepsChargedParticlesSource");
   muonIsoDepsChargedHadronsSource = iConfig.getParameter<edm::InputTag>("MuonIsoDepsChargedHadronsSource");
   muonIsoDepsNeutralHadronsSource = iConfig.getParameter<edm::InputTag>("MuonIsoDepsNeutralHadronsSource");
   muonIsoDepsPhotonsSource = iConfig.getParameter<edm::InputTag>("MuonIsoDepsPhotonsSource");
   muonIsoDepsPUSource = iConfig.getParameter<edm::InputTag>("MuonIsoDepsPUSource");

   tauByDecayModeFindingSource = iConfig.getParameter<edm::InputTag>("TauByDecayModeFindingSource");
   tauByCombinedIsolationDBSumPtCorr3HitsRawSource = iConfig.getParameter<edm::InputTag>("TauByCombinedIsolationDBSumPtCorr3HitsRawSource");
   tauByLooseCombinedIsolationDBSumPtCorr3HitsSource = iConfig.getParameter<edm::InputTag>("TauByLooseCombinedIsolationDBSumPtCorr3HitsSource");
   tauByMediumCombinedIsolationDBSumPtCorr3HitsSource = iConfig.getParameter<edm::InputTag>("TauByMediumCombinedIsolationDBSumPtCorr3HitsSource");
   tauByTightCombinedIsolationDBSumPtCorr3HitsSource = iConfig.getParameter<edm::InputTag>("TauByTightCombinedIsolationDBSumPtCorr3HitsSource");
   tauByLooseElectronRejectionSource = iConfig.getParameter<edm::InputTag>("TauByLooseElectronRejectionSource");
   tauByLooseElectronRejectionMVA3Source = iConfig.getParameter<edm::InputTag>("TauByLooseElectronRejectionMVA3Source");
   tauByMediumElectronRejectionMVA3Source = iConfig.getParameter<edm::InputTag>("TauByMediumElectronRejectionMVA3Source");
   tauByTightElectronRejectionMVA3Source = iConfig.getParameter<edm::InputTag>("TauByTightElectronRejectionMVA3Source");
   tauByLooseMuonRejectionSource = iConfig.getParameter<edm::InputTag>("TauByLooseMuonRejectionSource");
   tauByTightMuonRejectionSource = iConfig.getParameter<edm::InputTag>("TauByTightMuonRejectionSource");

   caloMetSource = iConfig.getParameter<edm::InputTag>("CaloMetSource");
   pfMetSource = iConfig.getParameter<edm::InputTag>("PfMetSource");
   pfType1CorrectedMetSource = iConfig.getParameter<edm::InputTag>("PfType1CorrectedMetSource");
   mvaPfMetSource = iConfig.getParameter<edm::InputTag>("MvaPfMetSource");

   beamSpotSource = iConfig.getParameter<edm::InputTag>("BeamSpotSource");
   vertexSource = iConfig.getParameter<edm::InputTag>("VertexSource");
   pileupSummaryInfoSource = iConfig.getParameter<edm::InputTag>("PileupSummaryInfoSource");
   if(isEmbedded) origVertexSource = iConfig.getParameter<edm::InputTag>("OrigVertexSource");

   genParticlesSource = iConfig.getParameter<edm::InputTag>("GenParticlesSource");;
   if(isEmbedded) origGenParticlesSource = iConfig.getParameter<edm::InputTag>("OrigGenParticlesSource");;
}

LepTauNTupleProducer::~LepTauNTupleProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

// ------------ method called for each event  ------------
void
LepTauNTupleProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   run = iEvent.id().run();
   lumi = iEvent.getLuminosityBlock().luminosityBlock();
   event = iEvent.id().event();

   // Trigger
   edm::Handle<edm::TriggerResults> triggerResultsPtr, origTriggerResultsPtr;
   iEvent.getByLabel(triggerSource, triggerResultsPtr);
   if(isEmbedded) iEvent.getByLabel(origTriggerSource, origTriggerResultsPtr);

   if(isEmbedded)
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

   // BeamSpot
   edm::Handle<reco::BeamSpot> beamSpotPtr;
   iEvent.getByLabel(beamSpotSource, beamSpotPtr);

   // Primary vertices
   edm::Handle<std::vector<reco::Vertex> > vertexPtr;
   iEvent.getByLabel(vertexSource, vertexPtr);
   nPV = vertexPtr->size();

   beamSpotX = beamSpotPtr->position().x();
   beamSpotY = beamSpotPtr->position().y();
   beamSpotZ = beamSpotPtr->position().z();

   if(!vertexPtr->empty())
   {
      vertexX = (*vertexPtr)[0].x();
      vertexY = (*vertexPtr)[0].y();
      vertexZ = (*vertexPtr)[0].z();
      vertexChi2 = (*vertexPtr)[0].chi2();
      vertexNdof = (*vertexPtr)[0].ndof();
      vertexFake = (*vertexPtr)[0].isFake();
      vertexValid = (*vertexPtr)[0].isValid();
   }
   else
   {
      vertexX = beamSpotPtr->position().x();
      vertexY = beamSpotPtr->position().y();
      vertexZ = beamSpotPtr->position().z();
      vertexChi2 = 0.0f;
      vertexNdof = 0.0f;
      vertexFake = true;
      vertexValid = false;
   }

   if(!isData)
   {
      edm::Handle<std::vector<PileupSummaryInfo> > puInfoPtr;
      iEvent.getByLabel(pileupSummaryInfoSource, puInfoPtr);

      nTrueInteractions = -1.f;
      for(std::vector<PileupSummaryInfo>::const_iterator iter = puInfoPtr->begin(); iter != puInfoPtr->end(); ++iter)
        if(iter->getBunchCrossing() == 0)
          nTrueInteractions = iter->getTrueNumInteractions();

      if(nTrueInteractions < 0.f)
        throw cms::Exception("LepTauNTupleProducer") << "No Pileup information for BX = 0 present";
   }
   else
   {
      nTrueInteractions = -1.f;
   }

   // Orig muons
   if(isEmbedded)
   {
      edm::Handle<std::vector<reco::CompositeCandidate> > origMuonPtr;
      iEvent.getByLabel(origMuonSource, origMuonPtr);

      if(origMuonPtr.isValid())
      {
         const reco::CompositeCandidate& cand = (*origMuonPtr)[0];
         assert(cand.numberOfDaughters() == 2);

         const reco::Candidate& posMuon = cand.daughter(0)->charge() > 0 ? *cand.daughter(0) : *cand.daughter(1);
         const reco::Candidate& negMuon = cand.daughter(0)->charge() < 0 ? *cand.daughter(0) : *cand.daughter(1);

         origPosMuonPt = posMuon.pt();
         origPosMuonEta = posMuon.eta();
         origPosMuonPhi = posMuon.phi();
         origPosMuonE = posMuon.energy();

         origNegMuonPt = negMuon.pt();
         origNegMuonEta = negMuon.eta();
         origNegMuonPhi = negMuon.phi();
         origNegMuonE = negMuon.energy();
      }
      else
      {
         edm::Handle<std::vector<reco::Muon> > origMuonPtr;
         iEvent.getByLabel(edm::InputTag("muons", "", "RECO"), origMuonPtr);
         edm::Handle<std::vector<reco::Vertex> > origVertexPtr;
         iEvent.getByLabel(origVertexSource, origVertexPtr);

         const std::pair<const reco::Muon*, const reco::Muon*> muons = getMuons(*origMuonPtr, origVertexPtr, false);
         const reco::Muon* posMuon = muons.first;
         const reco::Muon* negMuon = muons.second;

         // origMuons must always be available for embedded samples
         assert(posMuon != NULL && negMuon != NULL);

         origPosMuonPt = posMuon->pt();
         origPosMuonEta = posMuon->eta();
         origPosMuonPhi = posMuon->phi();
         origPosMuonE = posMuon->energy();

         origNegMuonPt = negMuon->pt();
         origNegMuonEta = negMuon->eta();
         origNegMuonPhi = negMuon->phi();
         origNegMuonE = negMuon->energy();
      }
   }

   // Reconstructed Electrons
   edm::Handle<std::vector<reco::GsfElectron> > electronPtr;
   iEvent.getByLabel(electronSource, electronPtr);

   edm::Handle<reco::ConversionCollection> conversionsPtr;
   iEvent.getByLabel(conversionsSource, conversionsPtr);

   edm::Handle<edm::ValueMap<float> > electronIDPtr;
   iEvent.getByLabel(electronIDSource, electronIDPtr);

   edm::Handle<edm::ValueMap<reco::IsoDeposit> > electronIsoDepsChargedParticlesPtr, electronIsoDepsChargedHadronsPtr, electronIsoDepsNeutralHadronsPtr, electronIsoDepsPhotonsPtr, electronIsoDepsPUPtr;
   iEvent.getByLabel(electronIsoDepsChargedParticlesSource, electronIsoDepsChargedParticlesPtr);
   iEvent.getByLabel(electronIsoDepsChargedHadronsSource, electronIsoDepsChargedHadronsPtr);
   iEvent.getByLabel(electronIsoDepsNeutralHadronsSource, electronIsoDepsNeutralHadronsPtr);
   iEvent.getByLabel(electronIsoDepsPhotonsSource, electronIsoDepsPhotonsPtr);
   iEvent.getByLabel(electronIsoDepsPUSource, electronIsoDepsPUPtr);

   std::vector<Electron> electrons;
   electrons.reserve(electronPtr->size());
   for(unsigned int i = 0; i < electronPtr->size(); ++i)
   {
      reco::GsfElectronRef ref(electronPtr, i);
      Electron elec(ref, conversionsPtr, beamSpotPtr, *electronIDPtr, *electronIsoDepsChargedParticlesPtr, *electronIsoDepsChargedHadronsPtr, *electronIsoDepsNeutralHadronsPtr, *electronIsoDepsPhotonsPtr, *electronIsoDepsPUPtr);
      electrons.push_back(elec);
   }

   // Reconstructed Muons
   edm::Handle<std::vector<reco::Muon> > muonPtr;
   iEvent.getByLabel(muonSource, muonPtr);

   edm::Handle<edm::ValueMap<reco::IsoDeposit> > muonIsoDepsChargedParticlesPtr, muonIsoDepsChargedHadronsPtr, muonIsoDepsNeutralHadronsPtr, muonIsoDepsPhotonsPtr, muonIsoDepsPUPtr;
   iEvent.getByLabel(muonIsoDepsChargedParticlesSource, muonIsoDepsChargedParticlesPtr);
   iEvent.getByLabel(muonIsoDepsChargedHadronsSource, muonIsoDepsChargedHadronsPtr);
   iEvent.getByLabel(muonIsoDepsNeutralHadronsSource, muonIsoDepsNeutralHadronsPtr);
   iEvent.getByLabel(muonIsoDepsPhotonsSource, muonIsoDepsPhotonsPtr);
   iEvent.getByLabel(muonIsoDepsPUSource, muonIsoDepsPUPtr);

   std::vector<Muon> muons;
   muons.reserve(muonPtr->size());
   for(unsigned int i = 0; i < muonPtr->size(); ++i)
   {
      reco::MuonRef ref(muonPtr, i);
      Muon muon(ref, *muonIsoDepsChargedParticlesPtr, *muonIsoDepsChargedHadronsPtr, *muonIsoDepsNeutralHadronsPtr, *muonIsoDepsPhotonsPtr, *muonIsoDepsPUPtr);
      muons.push_back(muon);
   }

   // Reconstructed Taus
   edm::Handle<std::vector<reco::PFTau> > tauPtr;
   iEvent.getByLabel(tauSource, tauPtr);

   edm::Handle<reco::PFTauDiscriminator> byDecayModeFinding;
   edm::Handle<reco::PFTauDiscriminator> byCombinedIsolationDBSumPtCorr3HitsRaw;
   edm::Handle<reco::PFTauDiscriminator> byLooseCombinedIsolationDBSumPtCorr3Hits;
   edm::Handle<reco::PFTauDiscriminator> byMediumCombinedIsolationDBSumPtCorr3Hits;
   edm::Handle<reco::PFTauDiscriminator> byTightCombinedIsolationDBSumPtCorr3Hits;
   edm::Handle<reco::PFTauDiscriminator> byLooseElectronRejection;
   edm::Handle<reco::PFTauDiscriminator> byLooseElectronRejectionMVA3;
   edm::Handle<reco::PFTauDiscriminator> byMediumElectronRejectionMVA3;
   edm::Handle<reco::PFTauDiscriminator> byTightElectronRejectionMVA3;
   edm::Handle<reco::PFTauDiscriminator> byLooseMuonRejection;
   edm::Handle<reco::PFTauDiscriminator> byTightMuonRejection;

   iEvent.getByLabel(tauByDecayModeFindingSource, byDecayModeFinding);
   iEvent.getByLabel(tauByCombinedIsolationDBSumPtCorr3HitsRawSource, byCombinedIsolationDBSumPtCorr3HitsRaw);
   iEvent.getByLabel(tauByLooseCombinedIsolationDBSumPtCorr3HitsSource, byLooseCombinedIsolationDBSumPtCorr3Hits);
   iEvent.getByLabel(tauByMediumCombinedIsolationDBSumPtCorr3HitsSource, byMediumCombinedIsolationDBSumPtCorr3Hits);
   iEvent.getByLabel(tauByTightCombinedIsolationDBSumPtCorr3HitsSource, byTightCombinedIsolationDBSumPtCorr3Hits);
   iEvent.getByLabel(tauByLooseElectronRejectionSource, byLooseElectronRejection);
   iEvent.getByLabel(tauByLooseElectronRejectionMVA3Source, byLooseElectronRejectionMVA3);
   iEvent.getByLabel(tauByMediumElectronRejectionMVA3Source, byMediumElectronRejectionMVA3);
   iEvent.getByLabel(tauByTightElectronRejectionMVA3Source, byTightElectronRejectionMVA3);
   iEvent.getByLabel(tauByLooseMuonRejectionSource, byLooseMuonRejection);
   iEvent.getByLabel(tauByTightMuonRejectionSource, byTightMuonRejection);

   std::vector<Tau> taus;
   taus.reserve(tauPtr->size());
   for(unsigned i = 0; i < tauPtr->size(); ++i)
   {
      reco::PFTauRef tauRef(tauPtr, i);

      Tau tau(tauRef, *byDecayModeFinding, *byCombinedIsolationDBSumPtCorr3HitsRaw, *byLooseCombinedIsolationDBSumPtCorr3Hits, *byMediumCombinedIsolationDBSumPtCorr3Hits, *byTightCombinedIsolationDBSumPtCorr3Hits, *byLooseElectronRejection, *byLooseElectronRejectionMVA3, *byMediumElectronRejectionMVA3, *byTightElectronRejectionMVA3, *byLooseMuonRejection, *byTightMuonRejection);

      taus.push_back(tau);
   }

   const Electron* electron = NULL;
   const Muon* muon = NULL;
   const Tau* tau = NULL;

   if(channel == ETAU)
   {
      const std::pair<const Electron*, const Tau*> pair = getLepTau(electrons, taus, vertexPtr, channel, false);
      electron = pair.first;
      tau = pair.second;
   }
   else
   {
      const std::pair<const Muon*, const Tau*> pair = getLepTau(muons, taus, vertexPtr, channel, false);
      muon = pair.first;
      tau = pair.second;
   }

   if(electron != NULL)
   {
      electronPt = electron->pt();
      electronEta = electron->eta();
      electronPhi = electron->phi();
      electronE = electron->energy();
      electronCharge = electron->charge();
      electronChargedParticlePfIso04 = electron->chargedParticleIso();
      electronChargedHadronPfIso04 = electron->chargedHadronIso();
      electronNeutralHadronEtPfIso04 = electron->neutralHadronIso();
      electronPhotonEtPfIso04 = electron->photonIso();
      electronPUPtPfIso04 = electron->puIso();
      electronQuality = isQuality(*electron, vertexPtr);
      electronHasConversion = electron->hasConversion();
      electronMva = electron->mvaId();
   }
   else
   {
      electronPt = electronEta = electronPhi = electronE = -1.0f;
      electronChargedParticlePfIso04 = electronChargedHadronPfIso04 = electronNeutralHadronEtPfIso04 = electronPhotonEtPfIso04 = electronPUPtPfIso04 = -1.0f;
      electronCharge = 0; electronQuality = false; electronHasConversion = true; electronMva = -1.0f;
   }

   if(muon != NULL)
   {
      muonPt = muon->pt();
      muonEta = muon->eta();
      muonPhi = muon->phi();
      muonE = muon->energy();
      muonCharge = muon->charge();
      muonChargedParticlePfIso04 = muon->chargedParticleIso();
      muonChargedHadronPfIso04 = muon->chargedHadronIso();
      muonNeutralHadronEtPfIso04 = muon->neutralHadronIso();
      muonPhotonEtPfIso04 = muon->photonIso();
      muonPUPtPfIso04 = muon->puIso();
      muonQuality = isQuality(*muon, vertexPtr);
   }
   else
   {
      // no muon reconstructed. Set all variables to -1 to indicate this.
      muonPt = muonEta = muonPhi = muonE = -1.0f;
      muonChargedParticlePfIso04 = muonChargedHadronPfIso04 = muonNeutralHadronEtPfIso04 = muonPhotonEtPfIso04 = muonPUPtPfIso04 = -1.0f;
      muonCharge = 0; muonQuality = false;
   }

   if(tau != NULL)
   {
      tauPt = tau->pt();
      tauEta = tau->eta();
      tauPhi = tau->phi();
      tauE = tau->energy();
      tauCharge = tau->charge();
      tauQuality = isQuality(*tau, vertexPtr, channel);
      tauDecayModeFinding = tau->decayModeFinding();
      tauIsolation3Hits = tau->isolation3Hits();
      tauLooseIsolation3Hits = tau->looseIsolation3Hits();
      tauMediumIsolation3Hits = tau->mediumIsolation3Hits();
      tauTightIsolation3Hits = tau->tightIsolation3Hits();
      tauAgainstMuonLoose = tau->againstMuonLoose();
      tauAgainstMuonTight = tau->againstMuonTight();
      tauAgainstElectronLoose = tau->againstElectronLoose();
      tauAgainstElectronLooseMVA3 = tau->againstElectronLooseMVA3();
      tauAgainstElectronMediumMVA3 = tau->againstElectronMediumMVA3();
      tauAgainstElectronTightMVA3 = tau->againstElectronTightMVA3();
   }
   else
   {
      // no tau reconstructed. Set all variables to -1 to indicate this.
      tauPt = tauEta = tauPhi = tauE = -1.0f;
      tauIsolation3Hits = -1.0f;
      tauCharge = 0;
      tauQuality = false;
      tauDecayModeFinding = tauLooseIsolation3Hits = tauMediumIsolation3Hits = tauTightIsolation3Hits = tauAgainstMuonLoose = tauAgainstMuonTight = tauAgainstElectronLoose = tauAgainstElectronLooseMVA3 = tauAgainstElectronMediumMVA3 = tauAgainstElectronTightMVA3 = false;
   }

   // Jets
   edm::Handle<std::vector<reco::PFJet> > jetPtr;
   edm::Handle<reco::JetTagCollection> bTagPtr;
   iEvent.getByLabel(jetSource, jetPtr);
   iEvent.getByLabel(btagSource, bTagPtr);

   nJets = 0;
   for(unsigned int i = 0; i < jetPtr->size() && nJets < 20; ++i)
   {
      const reco::Jet* jet = &(*jetPtr)[i];
      if((muon != NULL && ROOT::Math::VectorUtil::DeltaR(jet->p4(), muon->p4()) < 0.3) || (tau != NULL && ROOT::Math::VectorUtil::DeltaR(jet->p4(), tau->p4()) < 0.3))
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
      if(bTagIndex == bTagPtr->size()) throw cms::Exception("LepTauNTupleProducer") << "No BTag found for jet!";

      ++nJets;
   }

   // Generated taus
   const reco::Candidate* genLepVis = NULL;
   const reco::Candidate* genTauVis = NULL;

   const reco::Candidate* genTauLep = NULL;
   const reco::Candidate* genTauHad = NULL;

   muonGenPt = muonGenEta = muonGenPhi = muonGenE = -1.0f;
   muonGenFinalPt = muonGenFinalEta = muonGenFinalPhi = muonGenFinalE = -1.0f;
   muonGenPtVis = muonGenEtaVis = muonGenPhiVis = muonGenEVis = -1.0f;
   muonGenCharge = 0;
   electronGenPt = electronGenEta = electronGenPhi = electronGenE = -1.0f;
   electronGenFinalPt = electronGenFinalEta = electronGenFinalPhi = electronGenFinalE = -1.0f;
   electronGenPtVis = electronGenEtaVis = electronGenPhiVis = electronGenEVis = -1.0f;
   electronGenCharge = 0;
   tauGenPt = tauGenEta = tauGenPhi = tauGenE = -1.0f;
   tauGenFinalPt = tauGenFinalEta = tauGenFinalPhi = tauGenFinalE = -1.0f;
   tauGenPtVis = tauGenEtaVis = tauGenPhiVis = tauGenEVis = -1.0f;
   tauGenCharge = tauGenDecayPDG = 0;

   if(!isData)
   {
      edm::Handle<std::vector<reco::GenParticle> > genParticlesPtr;
      iEvent.getByLabel(genParticlesSource, genParticlesPtr);

      // Find exactly two opposite sign ME taus
      bool havePosTau = false;
      bool haveNegTau = false;
      const reco::Candidate* posGenTau = NULL;
      const reco::Candidate* negGenTau = NULL;
      for(unsigned int i = 0; i < genParticlesPtr->size(); ++i)
      {
         const reco::Candidate* cand = &(*genParticlesPtr)[i];

	 if(abs(cand->pdgId()) != 15)
            continue;
         if(cand->mother(0)->pdgId() != 23 && cand->mother(0)->pdgId() != 25) continue;

         if(cand->pdgId() < 0)
         {
            if(havePosTau)
               posGenTau = NULL; /* more than one */
            else
               posGenTau = cand;
            havePosTau = true;
         }
         else
         {
            if(haveNegTau)
               negGenTau = NULL; /* more than one */
            else
               negGenTau = cand;
            haveNegTau = true;
         }
      }

      if(posGenTau != NULL && negGenTau != NULL)
      {
         const reco::Candidate::LorentzVector posGenTauStable = getStableParticle(posGenTau);
         const reco::Candidate::LorentzVector negGenTauStable = getStableParticle(negGenTau);
         const reco::Candidate* posGenTauVis = getVisTauDecay(*posGenTau);
         const reco::Candidate* negGenTauVis = getVisTauDecay(*negGenTau);

         assert(posGenTauVis != NULL); // must decay somehow
         assert(negGenTauVis != NULL); // must decay somehow

         const int lepPDG = ((channel == ETAU) ? 11 : 13);

         const reco::Candidate::LorentzVector* genTauLepFinal = NULL;
         const reco::Candidate::LorentzVector* genTauHadFinal = NULL;

         if(abs(posGenTauVis->pdgId()) == lepPDG) { genLepVis = posGenTauVis; genTauLep = posGenTau; genTauLepFinal = &posGenTauStable; }
         else if(abs(posGenTauVis->pdgId()) != 11 && abs(posGenTauVis->pdgId()) != 13) { genTauVis = posGenTauVis; genTauHad = posGenTau; genTauHadFinal = &posGenTauStable; }

         if(abs(negGenTauVis->pdgId()) == lepPDG) { genLepVis = negGenTauVis; genTauLep = negGenTau; genTauLepFinal = &negGenTauStable; }
         else if(abs(negGenTauVis->pdgId()) != 11 && abs(negGenTauVis->pdgId()) != 13) { genTauVis = negGenTauVis; genTauHad = negGenTau; genTauHadFinal = &negGenTauStable; }

         if(genLepVis != NULL && genTauVis != NULL)
         {
            if(channel == ETAU)
            {
               electronGenPt = genTauLep->pt();
               electronGenEta = genTauLep->eta();
               electronGenPhi = genTauLep->phi();
               electronGenE = genTauLep->energy();
               electronGenFinalPt = genTauLepFinal->pt();
               electronGenFinalEta = genTauLepFinal->eta();
               electronGenFinalPhi = genTauLepFinal->phi();
               electronGenFinalE = genTauLepFinal->energy();
               electronGenPtVis = genLepVis->pt();
               electronGenEtaVis = genLepVis->eta();
               electronGenPhiVis = genLepVis->phi();
               electronGenEVis = genLepVis->energy();
               electronGenCharge = genTauLep->pdgId() > 0 ? -1 : +1;
            }
            else
            {
               muonGenPt = genTauLep->pt();
               muonGenEta = genTauLep->eta();
               muonGenPhi = genTauLep->phi();
               muonGenE = genTauLep->energy();
               muonGenFinalPt = genTauLepFinal->pt();
               muonGenFinalEta = genTauLepFinal->eta();
               muonGenFinalPhi = genTauLepFinal->phi();
               muonGenFinalE = genTauLepFinal->energy();
               muonGenPtVis = genLepVis->pt();
               muonGenEtaVis = genLepVis->eta();
               muonGenPhiVis = genLepVis->phi();
               muonGenEVis = genLepVis->energy();
               muonGenCharge = genTauLep->pdgId() > 0 ? -1 : +1;
            }

            tauGenPt = genTauHad->pt();
            tauGenEta = genTauHad->eta();
            tauGenPhi = genTauHad->phi();
            tauGenE = genTauHad->energy();
            tauGenFinalPt = genTauHadFinal->pt();
            tauGenFinalEta = genTauHadFinal->eta();
            tauGenFinalPhi = genTauHadFinal->phi();
            tauGenFinalE = genTauHadFinal->energy();
            tauGenPtVis = genTauVis->pt();
            tauGenEtaVis = genTauVis->eta();
            tauGenPhiVis = genTauVis->phi();
            tauGenEVis = genTauVis->energy();
            tauGenCharge = genTauHad->pdgId() > 0 ? -1 : +1;
            tauGenDecayPDG = genTauVis->pdgId();
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
   if(isEmbedded && !isData)
   {
      edm::Handle<std::vector<reco::GenParticle> > origGenParticlesPtr;
      iEvent.getByLabel(origGenParticlesSource, origGenParticlesPtr);

      // Find exactly two opposite sign ME muons
      for(unsigned int i = 0; i < origGenParticlesPtr->size(); ++i)
      {
         const reco::Candidate* cand = &(*origGenParticlesPtr)[i];
         if(cand->numberOfMothers() != 1 || cand->mother(0)->pdgId() != 23) continue;
	 if(fabs(cand->pdgId()) != 11 && fabs(cand->pdgId()) != 13 && fabs(cand->pdgId()) != 15)
            continue;

         const reco::Candidate::LorentzVector stable = getStableParticle(cand);

         if(cand->pdgId() < 0)
	 {
            origGenPosChargedLepPt = cand->pt();
            origGenPosChargedLepEta = cand->eta();
            origGenPosChargedLepPhi = cand->phi();
            origGenPosChargedLepE = cand->energy();
            origGenPosChargedLepFinalPt = stable.pt();
            origGenPosChargedLepFinalEta = stable.eta();
            origGenPosChargedLepFinalPhi = stable.phi();
            origGenPosChargedLepFinalE = stable.energy();
            origGenPosChargedLepPDG = cand->pdgId();
	 }
         else
	 {
            origGenNegChargedLepPt = cand->pt();
            origGenNegChargedLepEta = cand->eta();
            origGenNegChargedLepPhi = cand->phi();
            origGenNegChargedLepE = cand->energy();
            origGenNegChargedLepFinalPt = stable.pt();
            origGenNegChargedLepFinalEta = stable.eta();
            origGenNegChargedLepFinalPhi = stable.phi();
            origGenNegChargedLepFinalE = stable.energy();
            origGenNegChargedLepPDG = cand->pdgId();
	 }
      }
   }

   // Get weight for embedded events
   minVisPtFilterWeight = 1.0f;
   tauSpinWeight = 1.0f;
   zmumuEvtSelEffCorrWeight = 1.0f;
   muonRadiationWeight = 1.0f;
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
         edm::Handle<double> hTauSpinnerWT;
         iEvent.getByLabel(edm::InputTag("TauSpinnerReco", "TauSpinnerWT", "EmbeddedSPIN"), hTauSpinnerWT);
         tauSpinWeight = *hTauSpinnerWT;

         edm::Handle<double> hZmumuEvtSelEffCorrWeight;
         iEvent.getByLabel(edm::InputTag("ZmumuEvtSelEffCorrWeightProducer", "weight", "EmbeddedRECO"), hZmumuEvtSelEffCorrWeight);
	 if(!hZmumuEvtSelEffCorrWeight.isValid()) return; // skip event(?)
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
   }

   // Fill primary vertex histograms for PileUp reweighting
   const float totalWeight = minVisPtFilterWeight * tauSpinWeight * zmumuEvtSelEffCorrWeight * muonRadiationWeight;
   h_nPV->Fill(nPV, totalWeight);
   if(!isData) h_nTrueInteractions->Fill(nTrueInteractions, totalWeight);

   //  Fill event if either:
   //   * genlevel muon+tau are found within the orig di-muon embedding acceptance
   //   * genlevel muon+tau are found with visible pT > 10 GeV
   //   * quality muon > 15 GeV AND quality and isolated tau > 20 GeV are found
   bool keepEvent = false;
   if(genTauHad != NULL && genTauLep != NULL && std::max(genTauHad->pt(), genTauLep->pt()) > 17. && std::min(genTauHad->pt(), genTauLep->pt()) > 8. && fabs(genTauHad->eta()) < 2.5 && fabs(genTauLep->eta()) < 2.5)
      keepEvent = true;
   if(genLepVis != NULL && genTauVis != NULL && genLepVis->pt() > 10. && genTauVis->pt() > 10.)
      keepEvent = true;
   else if(muon != NULL && tau != NULL && isQuality(*muon, vertexPtr) && isQuality(*tau, vertexPtr, channel) && tau->looseIsolation3Hits() && muon->pt() > 15. && tau->pt() > 20.)
      keepEvent = true;
   else if(electron != NULL && tau != NULL && isQuality(*electron, vertexPtr) && isQuality(*tau, vertexPtr, channel) && tau->looseIsolation3Hits() && electron->pt() > 15. && tau->pt() > 20.)
      keepEvent = true;

   if(keepEvent)
      tree->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void 
LepTauNTupleProducer::beginJob()
{
  edm::Service<TFileService> fs;

  tree = fs->make<TTree>("MuonTree", "MuonTree");

  tree->Branch("Run", &run, "Run/i");
  tree->Branch("Lumi", &lumi, "Lumi/i");
  tree->Branch("Event", &event, "Event/i");

  tree->Branch("MinVisPtFilterWeight", &minVisPtFilterWeight, "MinVisPtFilterWeight/F");
  tree->Branch("TauSpinWeight", &tauSpinWeight, "TauSpinWeight/F");
  tree->Branch("ZmumuEvtSelEffCorrWeight", &zmumuEvtSelEffCorrWeight, "ZmumuEvtSelEffCorrWeight/F");
  tree->Branch("MuonRadiationWeight", &muonRadiationWeight, "MuonRadiationWeight/F");
  tree->Branch("MuonRadiationFilter", &muonRadiationFilter, "MuonRadiationFilter/O");
  tree->Branch("MuonRadiationFilter2Sel1", &muonRadiationFilter2Sel1, "MuonRadiationFilter2Sel1/O");
  tree->Branch("MuonRadiationFilter2Sel2", &muonRadiationFilter2Sel2, "MuonRadiationFilter2Sel2/O");
  tree->Branch("MuonRadiationFilter2Sel3", &muonRadiationFilter2Sel3, "MuonRadiationFilter2Sel3/O");

  tree->Branch("NPV", &nPV, "NPV/i");
  if(!isData) tree->Branch("NTrueInteractions", &nTrueInteractions, "NTrueInteractions/F");

  if(isEmbedded)
  {
    tree->Branch("OrigHltIsoMu17", &origHltIsoMu17, "OrigHltIsoMu17/O");
    tree->Branch("OrigHltDoubleMu7", &origHltDoubleMu7, "OrigHltDoubleMu7/O");
    tree->Branch("OrigHltMu13Mu8", &origHltMu13Mu8, "OrigHltMu13Mu8/O");
    tree->Branch("OrigHltMu17Mu8", &origHltMu17Mu8, "OrigHltMu17Mu8/O");

    tree->Branch("OrigPosMuonPt", &origPosMuonPt, "OrigPosMuonPt/F");
    tree->Branch("OrigPosMuonEta", &origPosMuonEta, "OrigPosMuonEta/F");
    tree->Branch("OrigPosMuonPhi", &origPosMuonPhi, "OrigPosMuonPhi/F");
    tree->Branch("OrigPosMuonE", &origPosMuonE, "OrigPosMuonE/F");
    tree->Branch("OrigNegMuonPt", &origNegMuonPt, "OrigNegMuonPt/F");
    tree->Branch("OrigNegMuonEta", &origNegMuonEta, "OrigNegMuonEta/F");
    tree->Branch("OrigNegMuonPhi", &origNegMuonPhi, "OrigNegMuonPhi/F");
    tree->Branch("OrigNegMuonE", &origNegMuonE, "OrigNegMuonE/F");
  }

  if(isEmbedded && !isData)
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

  tree->Branch("BeamSpotX", &beamSpotX, "BeamSpotX/F");
  tree->Branch("BeamSpotY", &beamSpotY, "BeamSpotY/F");
  tree->Branch("BeamSpotZ", &beamSpotZ, "BeamSpotZ/F");

  tree->Branch("VertexX", &vertexX, "VertexX/F");
  tree->Branch("VertexY", &vertexY, "VertexY/F");
  tree->Branch("VertexZ", &vertexZ, "VertexZ/F");
  tree->Branch("VertexChi2", &vertexChi2, "VertexChi2/F");
  tree->Branch("VertexNdof", &vertexNdof, "VertexNdof/F");
  tree->Branch("VertexFake", &vertexFake, "VertexFake/O");
  tree->Branch("VertexValid", &vertexValid, "VertexValid/O");

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

  if(channel == ETAU)
  {
    tree->Branch("ElectronPt", &electronPt, "ElectronPt/F");
    tree->Branch("ElectronEta", &electronEta, "ElectronEta/F");
    tree->Branch("ElectronPhi", &electronPhi, "ElectronPhi/F");
    tree->Branch("ElectronE", &electronE, "ElectronE/F");
    tree->Branch("ElectronCharge", &electronCharge, "ElectronCharge/I");
    tree->Branch("ElectronChargedParticlePtPfIso04", &electronChargedParticlePfIso04, "ElectronChargedParticlePtPfIso04/F");
    tree->Branch("ElectronChargedHadronPtPfIso04", &electronChargedHadronPfIso04, "ElectronChargedHadronPtPfIso04/F");
    tree->Branch("ElectronNeutralHadronEtPfIso04", &electronNeutralHadronEtPfIso04, "ElectronNeutralHadronEtPfIso04/F");
    tree->Branch("ElectronPhotonEtPfIso04", &electronPhotonEtPfIso04, "ElectronPhotonEtPfIso04/F");
    tree->Branch("ElectronPUPtPfIso04", &electronPUPtPfIso04, "ElectronPUPtPfIso04/F");
    tree->Branch("ElectronQuality", &electronQuality, "ElectronQuality/O");
    tree->Branch("ElectronHasConversion", &electronHasConversion, "ElectronHasConversion/O");
    tree->Branch("ElectronMVA", &electronMva, "ElectronMVA/F");
  }
  else
  {
    tree->Branch("MuonPt", &muonPt, "MuonPt/F");
    tree->Branch("MuonEta", &muonEta, "MuonEta/F");
    tree->Branch("MuonPhi", &muonPhi, "MuonPhi/F");
    tree->Branch("MuonE", &muonE, "MuonE/F");
    tree->Branch("MuonCharge", &muonCharge, "MuonCharge/I");
    tree->Branch("MuonChargedParticlePtPfIso04", &muonChargedParticlePfIso04, "MuonChargedParticlePtPfIso04/F");
    tree->Branch("MuonChargedHadronPtPfIso04", &muonChargedHadronPfIso04, "MuonChargedHadronPtPfIso04/F");
    tree->Branch("MuonNeutralHadronEtPfIso04", &muonNeutralHadronEtPfIso04, "MuonNeutralHadronEtPfIso04/F");
    tree->Branch("MuonPhotonEtPfIso04", &muonPhotonEtPfIso04, "MuonPhotonEtPfIso04/F");
    tree->Branch("MuonPUPtPfIso04", &muonPUPtPfIso04, "MuonPUPtPfIso04/F");
    tree->Branch("MuonQuality", &muonQuality, "MuonQuality/O");
  }

  tree->Branch("TauPt", &tauPt, "TauPt/F");
  tree->Branch("TauEta", &tauEta, "TauEta/F");
  tree->Branch("TauPhi", &tauPhi, "TauPhi/F");
  tree->Branch("TauE", &tauE, "TauE/F");
  tree->Branch("TauCharge", &tauCharge, "TauCharge/I");
  tree->Branch("TauQuality", &tauQuality, "TauQuality/O");
  tree->Branch("TauDecayModeFinding", &tauDecayModeFinding, "TauDecayModeFinding/O");
  tree->Branch("TauIsolation3Hits", &tauIsolation3Hits, "TauIsolation3Hits/F");
  tree->Branch("TauLooseIsolation3Hits", &tauLooseIsolation3Hits, "TauLooseIsolation3Hits/O");
  tree->Branch("TauMediumIsolation3Hits", &tauMediumIsolation3Hits, "TauMediumIsolation3Hits/O");
  tree->Branch("TauTightIsolation3Hits", &tauTightIsolation3Hits, "TauTightIsolation3Hits/O");
  tree->Branch("TauAgainstMuonLoose", &tauAgainstMuonLoose, "TauAgainstMuonLoose/O");
  tree->Branch("TauAgainstMuonTight", &tauAgainstMuonTight, "TauAgainstMuonTight/O");
  tree->Branch("TauAgainstElectronLoose", &tauAgainstElectronLoose, "TauAgainstElectronLoose/O");
  tree->Branch("TauAgainstElectronLooseMVA3", &tauAgainstElectronLooseMVA3, "TauAgainstElectronLooseMVA3/O");
  tree->Branch("TauAgainstElectronMediumMVA3", &tauAgainstElectronMediumMVA3, "TauAgainstElectronMediumMVA3/O");
  tree->Branch("TauAgainstElectronTightMVA3", &tauAgainstElectronTightMVA3, "TauAgainstElectronTightMVA3/O");

  tree->Branch("NJets", &nJets, "NJets/i");
  tree->Branch("JetPt", jetPt, "JetPt[NJets]/F");
  tree->Branch("JetEta", jetEta, "JetEta[NJets]/F");
  tree->Branch("JetPhi", jetPhi, "JetPhi[NJets]/F");
  tree->Branch("JetEnergy", jetEnergy, "JetEnergy[NJets]/F");
  tree->Branch("JetBTag", jetBTag, "JetBTag[NJets]/F");

  if(!isData)
  {
    if(channel == ETAU)
    {
      tree->Branch("ElectronGenPt", &electronGenPt, "ElectronGenPt/F");
      tree->Branch("ElectronGenEta", &electronGenEta, "ElectronGenEta/F");
      tree->Branch("ElectronGenPhi", &electronGenPhi, "ElectronGenPhi/F");
      tree->Branch("ElectronGenE", &electronGenE, "ElectronGenE/F");
      tree->Branch("ElectronGenFinalPt", &electronGenFinalPt, "ElectronGenFinalPt/F");
      tree->Branch("ElectronGenFinalEta", &electronGenFinalEta, "ElectronGenFinalEta/F");
      tree->Branch("ElectronGenFinalPhi", &electronGenFinalPhi, "ElectronGenFinalPhi/F");
      tree->Branch("ElectronGenFinalE", &electronGenFinalE, "ElectronGenFinalE/F");
      tree->Branch("ElectronGenPtVis", &electronGenPtVis, "ElectronGenPtVis/F");
      tree->Branch("ElectronGenEtaVis", &electronGenEtaVis, "ElectronGenEtaVis/F");
      tree->Branch("ElectronGenPhiVis", &electronGenPhiVis, "ElectronGenPhiVis/F");
      tree->Branch("ElectronGenEVis", &electronGenEVis, "ElectronGenEVis/F");
      tree->Branch("ElectronGenCharge", &electronGenCharge, "ElectronGenCharge/I");
    }
    else
    {
      tree->Branch("MuonGenPt", &muonGenPt, "MuonGenPt/F");
      tree->Branch("MuonGenEta", &muonGenEta, "MuonGenEta/F");
      tree->Branch("MuonGenPhi", &muonGenPhi, "MuonGenPhi/F");
      tree->Branch("MuonGenE", &muonGenE, "MuonGenE/F");
      tree->Branch("MuonGenFinalPt", &muonGenFinalPt, "MuonGenFinalPt/F");
      tree->Branch("MuonGenFinalEta", &muonGenFinalEta, "MuonGenFinalEta/F");
      tree->Branch("MuonGenFinalPhi", &muonGenFinalPhi, "MuonGenFinalPhi/F");
      tree->Branch("MuonGenFinalE", &muonGenFinalE, "MuonGenFinalE/F");
      tree->Branch("MuonGenPtVis", &muonGenPtVis, "MuonGenPtVis/F");
      tree->Branch("MuonGenEtaVis", &muonGenEtaVis, "MuonGenEtaVis/F");
      tree->Branch("MuonGenPhiVis", &muonGenPhiVis, "MuonGenPhiVis/F");
      tree->Branch("MuonGenEVis", &muonGenEVis, "MuonGenEVis/F");
      tree->Branch("MuonGenCharge", &muonGenCharge, "MuonGenCharge/I");
    }

    tree->Branch("TauGenPt", &tauGenPt, "TauGenPt/F");
    tree->Branch("TauGenEta", &tauGenEta, "TauGenEta/F");
    tree->Branch("TauGenPhi", &tauGenPhi, "TauGenPhi/F");
    tree->Branch("TauGenE", &tauGenE, "TauGenE/F");
    tree->Branch("TauGenFinalPt", &tauGenFinalPt, "TauGenFinalPt/F");
    tree->Branch("TauGenFinalEta", &tauGenFinalEta, "TauGenFinalEta/F");
    tree->Branch("TauGenFinalPhi", &tauGenFinalPhi, "TauGenFinalPhi/F");
    tree->Branch("TauGenFinalE", &tauGenFinalE, "TauGenFinalE/F");
    tree->Branch("TauGenPtVis", &tauGenPtVis, "TauGenPtVis/F");
    tree->Branch("TauGenEtaVis", &tauGenEtaVis, "TauGenEtaVis/F");
    tree->Branch("TauGenPhiVis", &tauGenPhiVis, "TauGenPhiVis/F");
    tree->Branch("TauGenEVis", &tauGenEVis, "TauGenEVis/F");
    tree->Branch("TauGenCharge", &tauGenCharge, "TauGenCharge/I");
    tree->Branch("TauGenDecayPDG", &tauGenDecayPDG, "TauGenDecayPDG/I");
  }

  h_nPV = fs->make<TH1F>("h_nPV", "Number of primary vertices", 100, -0.5, 99.5);
  if(!isData) h_nTrueInteractions = fs->make<TH1F>("h_nTrueInteractions", "Number of true interactions", 100, -0.5, 99.5);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
LepTauNTupleProducer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
LepTauNTupleProducer::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
  bool changed = false;
  hltConfiguration.init(iRun, iSetup, triggerEventSource.process(), changed);
}

// ------------ method called when ending the processing of a run  ------------
void 
LepTauNTupleProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
LepTauNTupleProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
LepTauNTupleProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
LepTauNTupleProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(LepTauNTupleProducer);
