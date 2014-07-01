# -*- coding: utf-8 -*-

import FWCore.ParameterSet.Config as cms
import os

muonRadiationFilter2 = cms.EDFilter("MuonRadiationFilter2",
    srcSelectedMuons = cms.InputTag(''), # CV: replaced in embeddingCustomizeAll.py
    srcPFCandidates = cms.InputTag('particleFlow'),
    # parameters for reconstruction of eta x phi strips                                
    seedPtECAL = cms.double(3.),
    seedDeltaR = cms.double(0.50), # wrt. muon                          
    stripCandidatesParticleIds = cms.vint32(2, 3, 4), # CV: 2 = e, 3 = mu, 4 = gamma, cf. DataFormats/ParticleFlowCandidate/interface/PFCandidate.h                                    
    stripEtaAssociationDistance = cms.double(0.05),
    stripPhiAssociationDistance = cms.double(0.2),
    # selection of muon FSR events                                 
    stripSelection = cms.VPSet(
        cms.PSet(
            minPt = cms.double(10.),
            maxDeltaR = cms.double(0.30), # wrt. muon
            maxHoverE = cms.double(0.25),
            applyMassWindowSelection = cms.bool(False)
        ),
        cms.PSet(
            minPt = cms.double(5.),
            maxDeltaR = cms.double(0.50), # wrt. muon
            maxHoverE = cms.double(0.25),
            applyMassWindowSelection = cms.bool(True)
        ),
        cms.PSet(
            minPt = cms.double(5.),
            maxDeltaR = cms.double(0.10), # wrt. muon
            maxHoverE = cms.double(0.50),
            applyMassWindowSelection = cms.bool(False)
        )
    ),
    # track isolation parameters
    dRvetoCone = cms.double(1.e-3),
    dRisoCone = cms.double(0.4),
    minTrackPt = cms.double(1.0),
    maxTransverseImpactParameter = cms.double(0.03), # wrt. muon vertex
    maxDeltaZ = cms.double(0.2),                     # wrt. muon vertex
    maximumSumPtCut = cms.double(3.0),
    # global flags
    invert = cms.bool(False),                                     
    filter = cms.bool(True),
    verbosity = cms.int32(0)
)

muonRadiationFilterSequence2 = cms.Sequence(
    muonRadiationFilter2
)
