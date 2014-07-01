import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        # New RH:
	#'/store/user/aburgmei/embedding/20140303-RHEmbed-Prod-MC-v4-MuTau-NoRot-WCaloNoise/rhembTauTau_Summer12_DR53X_PU_S10_START53_V7A_v2_ReplaceRecMuons_0_spinned.root'
	# PF Emb:
	'file:/nfs/dust/cms/user/aburgm/embedding/files/pfemb/A8240750-9BE8-E211-B46E-0023AEFDE9E4.root'
	# Normal MC:
	#'/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0002/FEF4E41A-46D4-E111-9594-0025B3E06424.root'
	#'/store/mc/Summer12_DR53X/GluGluToHToTauTau_M-125_8TeV-powheg-pythia6/AODSIM/PU_S10_START53_V7A-v1/0000/00E903E2-9FE9-E111-8B1E-003048FF86CA.root'
    ),
    skipEvents = cms.untracked.uint32(0)
)
process.TFileService = cms.Service("TFileService", fileName = cms.string('output.root'))

useMvaMet = False
isEmbedded = True
isRHEmbedded = False
isData = False
channel = 'mutau'

# For grid-control override:
if not "@IS_EMBEDDED@".startswith('@'):
	isEmbedded = ("@IS_EMBEDDED@".lower() == 'true')
if not "@IS_RH_EMBEDDED@".startswith('@'):
	isRHEmbedded = ("@IS_RH_EMBEDDED@".lower() == 'true')
if not "@IS_DATA@".startswith('@'):
	isData = ("@IS_DATA@".lower() == 'true')
if not "@CHANNEL@".startswith('@'):
	channel = "@CHANNEL@"

## Geometry and Detector Conditions
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

if isData:
	process.GlobalTag.globaltag = cms.string('FT_53_V21_AN4::All')
else:
	process.GlobalTag.globaltag = cms.string('START53_V23::All')

############################################
# Jet energy corrections

process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load('JetMETCorrections.Configuration.JetCorrectionProducers_cff')

if isData:
	jetSource = cms.InputTag("ak5PFJetsL1FastL2L3Residual")
else:
	jetSource = cms.InputTag("ak5PFJetsL1FastL2L3")

############################################
# BTagging
process.load('RecoBTag/Configuration/RecoBTag_cff')
process.load('RecoJets/JetAssociationProducers/ak5JTA_cff')

process.ak5JetTracksAssociatorAtVertex.jets = cms.InputTag("ak5PFJets") # TODO: Should this be corrected jets?
if isEmbedded and not isRHEmbedded:
	process.ak5JetTracksAssociatorAtVertex.tracks = cms.InputTag("tmfTracks")

if isEmbedded and not isRHEmbedded:
	trackSource = cms.InputTag("tmfTracks")
else:
	trackSource = cms.InputTag("generalTracks")

############################################
# MET
process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")

if useMvaMet:
	process.load("JetMETCorrections.METPUSubtraction.mvaPFMET_leptons_cff")
	if isData:
		process.calibratedAK5PFJetsForPFMEtMVA.correctors = cms.vstring("ak5PFL1FastL2L3Residual")
	else:
		process.calibratedAK5PFJetsForPFMEtMVA.correctors = cms.vstring("ak5PFL1FastL2L3")
	process.pfMEtMVA.verbosity = cms.int32(0)
	mvaPfMetSource = cms.InputTag("pfMEtMVA")
else:
	mvaPfMetSource = cms.InputTag("")

############################################
# Tau
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")

############################################
# Electron ID
process.load('EgammaAnalysis.ElectronTools.electronIdMVAProducer_cfi')

############################################
# Lepton Iso
process.load('CommonTools.ParticleFlow.pfNoPileUp_cff')
process.pfPileUp.PFCandidates = cms.InputTag('particleFlow')
process.pfNoPileUp.bottomCollection = cms.InputTag('particleFlow')

from CommonTools.ParticleFlow.ParticleSelectors.pfAllChargedHadrons_cfi import pfAllChargedHadrons
from CommonTools.ParticleFlow.ParticleSelectors.pfAllNeutralHadrons_cfi import pfAllNeutralHadrons
from CommonTools.ParticleFlow.ParticleSelectors.pfAllPhotons_cfi import pfAllPhotons
from CommonTools.ParticleFlow.ParticleSelectors.pfAllMuons_cfi import pfAllMuons
from CommonTools.ParticleFlow.ParticleSelectors.pfAllElectrons_cfi import pfAllElectrons
from CommonTools.ParticleFlow.Isolation.tools_cfi import *

pfAllChargedParticlesPdgIds = []
pfAllChargedParticlesPdgIds.extend(pfAllChargedHadrons.pdgId.value())
pfAllChargedParticlesPdgIds.extend(pfAllMuons.pdgId.value())
pfAllChargedParticlesPdgIds.extend(pfAllElectrons.pdgId.value())

process.pfChargedParticles = cms.EDFilter("PdgIdPFCandidateSelector", src = cms.InputTag('pfNoPileUp'), pdgId = cms.vint32(pfAllChargedParticlesPdgIds))
process.pfChargedHadrons = cms.EDFilter("PdgIdPFCandidateSelector", src = cms.InputTag('pfNoPileUp'), pdgId = cms.vint32(pfAllChargedHadrons.pdgId.value()))
process.pfNeutralHadrons = cms.EDFilter("PdgIdPFCandidateSelector", src = cms.InputTag('particleFlow'), pdgId = cms.vint32(pfAllNeutralHadrons.pdgId.value()))
process.pfPhotons = cms.EDFilter("PdgIdPFCandidateSelector", src = cms.InputTag('particleFlow'), pdgId = cms.vint32(pfAllPhotons.pdgId.value()))
process.pfPU = cms.EDFilter("PdgIdPFCandidateSelector", src = cms.InputTag('pfPileUp'), pdgId = cms.vint32(pfAllChargedParticlesPdgIds))

process.pfElIsoDepositChargedParticles = isoDepositReplace("gsfElectrons", "pfChargedParticles")
process.pfElIsoDepositChargedHadrons = isoDepositReplace("gsfElectrons", "pfChargedHadrons")
process.pfElIsoDepositNeutralHadrons = isoDepositReplace("gsfElectrons", "pfNeutralHadrons")
process.pfElIsoDepositPhotons = isoDepositReplace("gsfElectrons", "pfPhotons")
process.pfElIsoDepositPU = isoDepositReplace("gsfElectrons", "pfPU")

process.pfMuIsoDepositChargedParticles = isoDepositReplace("muons", "pfChargedParticles")
process.pfMuIsoDepositChargedHadrons = isoDepositReplace("muons", "pfChargedHadrons")
process.pfMuIsoDepositNeutralHadrons = isoDepositReplace("muons", "pfNeutralHadrons")
process.pfMuIsoDepositPhotons = isoDepositReplace("muons", "pfPhotons")
process.pfMuIsoDepositPU = isoDepositReplace("muons", "pfPU")

process.isoSequence = cms.Sequence(process.pfNoPileUpSequence * process.pfChargedParticles * process.pfChargedHadrons * process.pfNeutralHadrons * process.pfPhotons * process.pfPU * process.pfElIsoDepositChargedParticles * process.pfElIsoDepositChargedHadrons * process.pfElIsoDepositNeutralHadrons * process.pfElIsoDepositPhotons * process.pfElIsoDepositPU * process.pfMuIsoDepositChargedParticles * process.pfMuIsoDepositChargedHadrons * process.pfMuIsoDepositNeutralHadrons * process.pfMuIsoDepositPhotons * process.pfMuIsoDepositPU)

############################################
# Radiation filter
process.load("Validation.MCEmbedding.muonRadiationFilter2_cfi")

process.muonRadiationFilter2.srcSelectedMuons = cms.InputTag('goldenZmumuCandidatesGe2IsoMuons', '', 'EmbeddedRECO')
process.muonRadiationFilter2.srcPFCandidates = cms.InputTag('particleFlow', '', 'RECO')
process.muonRadiationFilter2.filter = cms.bool(False)

process.muonRadiationFilter2ResultStripSel1 = process.muonRadiationFilter2.clone(
	stripSelection = cms.VPSet(process.muonRadiationFilter2.stripSelection[0]))
process.muonRadiationFilter2ResultStripSel2 = process.muonRadiationFilter2.clone(
	stripSelection = cms.VPSet(process.muonRadiationFilter2.stripSelection[1]))
process.muonRadiationFilter2ResultStripSel3 = process.muonRadiationFilter2.clone(
	stripSelection = cms.VPSet(process.muonRadiationFilter2.stripSelection[2]))

############################################
# EDAnalyzer

process.TauAna = cms.EDAnalyzer('LepTauNTupleProducer',
   isEmbedded = cms.bool(isEmbedded),
   isRHEmbedded = cms.bool(isRHEmbedded),
   isData = cms.bool(isData),
   channel = cms.string(channel),

   TriggerSource = cms.InputTag("TriggerResults"),
   TriggerEventSource = cms.InputTag("hltTriggerSummaryAOD"),

   OrigTriggerSource = cms.InputTag("TriggerResults", "", "HLT"),

   ElectronSource = cms.InputTag("gsfElectrons"),
   MuonSource = cms.InputTag("muons"),
   TauSource = cms.InputTag("hpsPFTauProducer"),

   OrigMuonSource = cms.InputTag("goldenZmumuCandidatesGe2IsoMuons"),

   ConversionsSource = cms.InputTag("allConversions"),
   ElectronIDSource = cms.InputTag("mvaNonTrigV0"),
   ElectronIsoDepsChargedParticlesSource = cms.InputTag("pfElIsoDepositChargedParticles"),
   ElectronIsoDepsChargedHadronsSource = cms.InputTag("pfElIsoDepositChargedHadrons"),
   ElectronIsoDepsNeutralHadronsSource = cms.InputTag("pfElIsoDepositNeutralHadrons"),
   ElectronIsoDepsPhotonsSource = cms.InputTag("pfElIsoDepositPhotons"),
   ElectronIsoDepsPUSource = cms.InputTag("pfElIsoDepositPU"),

   MuonIsoDepsChargedParticlesSource = cms.InputTag("pfMuIsoDepositChargedParticles"),
   MuonIsoDepsChargedHadronsSource = cms.InputTag("pfMuIsoDepositChargedHadrons"),
   MuonIsoDepsNeutralHadronsSource = cms.InputTag("pfMuIsoDepositNeutralHadrons"),
   MuonIsoDepsPhotonsSource = cms.InputTag("pfMuIsoDepositPhotons"),
   MuonIsoDepsPUSource = cms.InputTag("pfMuIsoDepositPU"),

   TauByDecayModeFindingSource = cms.InputTag("hpsPFTauDiscriminationByDecayModeFinding"),
   TauByCombinedIsolationDBSumPtCorr3HitsRawSource = cms.InputTag("hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits"),
   TauByLooseCombinedIsolationDBSumPtCorr3HitsSource = cms.InputTag("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits"),
   TauByMediumCombinedIsolationDBSumPtCorr3HitsSource = cms.InputTag("hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits"),
   TauByTightCombinedIsolationDBSumPtCorr3HitsSource = cms.InputTag("hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits"),
   TauByLooseElectronRejectionSource = cms.InputTag("hpsPFTauDiscriminationByLooseElectronRejection"),
   TauByLooseElectronRejectionMVA3Source = cms.InputTag("hpsPFTauDiscriminationByMVA3LooseElectronRejection"),
   TauByMediumElectronRejectionMVA3Source = cms.InputTag("hpsPFTauDiscriminationByMVA3MediumElectronRejection"),
   TauByTightElectronRejectionMVA3Source = cms.InputTag("hpsPFTauDiscriminationByMVA3TightElectronRejection"),
   TauByLooseMuonRejectionSource = cms.InputTag("hpsPFTauDiscriminationByLooseMuonRejection"),
   TauByTightMuonRejectionSource = cms.InputTag("hpsPFTauDiscriminationByTightMuonRejection"),

   JetSource = jetSource,
   BTagSource = cms.InputTag("combinedSecondaryVertexBJetTags"),
   TrackSource = trackSource,
   PfNoPileupSource = cms.InputTag("pfNoPileUp"),
   PfPileupSource = cms.InputTag("pfPileUp"),

   CaloMetSource = cms.InputTag("metNoHF"),
   PfMetSource = cms.InputTag("pfMet"),
   PfType1CorrectedMetSource = cms.InputTag("pfType1CorrectedMet"),
   MvaPfMetSource = mvaPfMetSource,

   BeamSpotSource = cms.InputTag("offlineBeamSpot"),
   VertexSource = cms.InputTag("offlinePrimaryVertices"),
   PileupSummaryInfoSource = cms.InputTag("addPileupInfo", "" ,"HLT"),
   OrigVertexSource = cms.InputTag("offlinePrimaryVertices", "", "RECO"),

   GenParticlesSource = cms.InputTag("genParticles"),
   OrigGenParticlesSource = cms.InputTag("genParticles", "", "SIM"),
)

process.p = cms.Path(process.ak5PFJetsL1FastL2L3 * process.ak5JetTracksAssociatorAtVertex * process.btagging * process.recoTauClassicHPSSequence * process.mvaNonTrigV0 * process.isoSequence * process.producePFMETCorrections * process.TauAna)
if isData:
	process.p.replace(process.ak5PFJetsL1FastL2L3, process.ak5PFJetsL1FastL2L3Residual)
if useMvaMet:
	process.p.replace(process.producePFMETCorrections, process.producePFMETCorrections * process.pfMEtMVAsequence)
if isEmbedded and isRHEmbedded:
	process.p.replace(process.TauAna, process.muonRadiationFilter2ResultStripSel1 * process.muonRadiationFilter2ResultStripSel2 * process.muonRadiationFilter2ResultStripSel3 * process.TauAna)
