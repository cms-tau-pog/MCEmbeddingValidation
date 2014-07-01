import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    	# Normal MC:
        #'/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/00037C53-AAD1-E111-B1BE-003048D45F38.root'
	# PF embedded:
	#'/store/user/aburgmei/embedding/20140228-PFEmbed-MuMu/pfembMuMu_mc_Summer12_DR53X_PU_S10_START53_V7A_v2_ReplaceRecMuons_tau0_pt_5_embedded.root'
	# RH embedded v5:
	'/store/user/aburgmei/embedding/20140305-RHEmbed-Prod-MC-v5-MuMu/rhembTauTau_Summer12_DR53X_PU_S10_START53_V7A_v2_ReplaceRecMuons_500_embed_AOD.root'
	# RH GEN embedded v5:
	#'/store/user/aburgmei/embedding/20140321-RHEmbed-Prod-MC-v5-MuMu-Gen/rhembTauTau_Summer12_DR53X_PU_S10_START53_V7A_v2_ReplaceGenMuons_0_spinned.root'
    )
)
process.TFileService = cms.Service("TFileService", fileName = cms.string('output.root'))

useMvaMet = False
writeAllEvents = False
isEmbedded = True
isGenEmbedded = False
isRHEmbedded = True
isData = False

# For grid-control override:
if not "@WRITE_ALL_EVENTS@".startswith('@'):
	writeAllEvents = ("@WRITE_ALL_EVENTS@".lower() == 'true')
if not "@IS_EMBEDDED@".startswith('@'):
	isEmbedded = ("@IS_EMBEDDED@".lower() == 'true')
if not "@IS_RH_EMBEDDED@".startswith('@'):
	isRHEmbedded = ("@IS_RH_EMBEDDED@".lower() == 'true')
if not "@IS_GEN_EMBEDDED@".startswith('@'):
	isGenEmbedded = ("@IS_GEN_EMBEDDED@".lower() == 'true')
if not "@IS_DATA@".startswith('@'):
	isData = ("@IS_DATA@".lower() == 'true')

## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

if isData:
	process.GlobalTag.globaltag = cms.string('FT_53_V21_AN4::All')
else:
	process.GlobalTag.globaltag = cms.string('START53_V23::All')

############################################
# RhoNeutral calculation

from CommonTools.ParticleFlow.ParticleSelectors.pfAllChargedHadrons_cfi import pfAllChargedHadrons
from CommonTools.ParticleFlow.ParticleSelectors.pfAllNeutralHadrons_cfi import pfAllNeutralHadrons
from CommonTools.ParticleFlow.ParticleSelectors.pfAllPhotons_cfi import pfAllPhotons
from CommonTools.ParticleFlow.ParticleSelectors.pfAllMuons_cfi import pfAllMuons
from CommonTools.ParticleFlow.ParticleSelectors.pfAllElectrons_cfi import pfAllElectrons

pfNeutralCandPdgIds = []
pfNeutralCandPdgIds.extend(pfAllNeutralHadrons.pdgId.value())
pfNeutralCandPdgIds.extend(pfAllPhotons.pdgId.value())

process.pfNeutralCands = cms.EDFilter("PdgIdPFCandidateSelector", src = cms.InputTag('particleFlow'), pdgId = cms.vint32(pfNeutralCandPdgIds))

process.load('RecoJets.JetProducers.kt4PFJets_cfi')
process.kt6PFJetsNeutral = process.kt4PFJets.clone()
process.kt6PFJetsNeutral.rParam = cms.double(0.6)
process.kt6PFJetsNeutral.doAreaFastjet = cms.bool( True )
process.kt6PFJetsNeutral.doRhoFastjet = cms.bool( True )
process.kt6PFJetsNeutral.Rho_EtaMax = cms.double(2.5)
process.kt6PFJetsNeutral.src = cms.InputTag('pfNeutralCands')
process.rhoSequence = cms.Sequence(process.pfNeutralCands * process.kt6PFJetsNeutral)

############################################
# Isolation
from CommonTools.ParticleFlow.Isolation.tools_cfi import *

process.load('CommonTools.ParticleFlow.pfNoPileUp_cff')
process.pfPileUp.PFCandidates = cms.InputTag('particleFlow')
process.pfNoPileUp.bottomCollection = cms.InputTag('particleFlow')

pfAllChargedParticlesPdgIds = []
pfAllChargedParticlesPdgIds.extend(pfAllChargedHadrons.pdgId.value())
pfAllChargedParticlesPdgIds.extend(pfAllMuons.pdgId.value())
pfAllChargedParticlesPdgIds.extend(pfAllElectrons.pdgId.value())

process.pfChargedParticles = cms.EDFilter("PdgIdPFCandidateSelector", src = cms.InputTag('pfNoPileUp'), pdgId = cms.vint32(pfAllChargedParticlesPdgIds))
process.pfChargedHadrons = cms.EDFilter("PdgIdPFCandidateSelector", src = cms.InputTag('pfNoPileUp'), pdgId = cms.vint32(pfAllChargedHadrons.pdgId.value()))
process.pfNeutralHadrons = cms.EDFilter("PdgIdPFCandidateSelector", src = cms.InputTag('particleFlow'), pdgId = cms.vint32(pfAllNeutralHadrons.pdgId.value()))
process.pfPhotons = cms.EDFilter("PdgIdPFCandidateSelector", src = cms.InputTag('particleFlow'), pdgId = cms.vint32(pfAllPhotons.pdgId.value()))
process.pfPU = cms.EDFilter("PdgIdPFCandidateSelector", src = cms.InputTag('pfPileUp'), pdgId = cms.vint32(pfAllChargedParticlesPdgIds))

process.pfMuIsoDepositChargedParticles = isoDepositReplace("muons", "pfChargedParticles")
process.pfMuIsoDepositChargedHadrons = isoDepositReplace("muons", "pfChargedHadrons")
process.pfMuIsoDepositNeutralHadrons = isoDepositReplace("muons", "pfNeutralHadrons")
process.pfMuIsoDepositPhotons = isoDepositReplace("muons", "pfPhotons")
process.pfMuIsoDepositPU = isoDepositReplace("muons", "pfPU")

process.isoSequence = cms.Sequence(process.pfNoPileUpSequence * process.pfChargedParticles * process.pfChargedHadrons * process.pfNeutralHadrons * process.pfPhotons * process.pfPU * process.pfMuIsoDepositChargedParticles * process.pfMuIsoDepositChargedHadrons * process.pfMuIsoDepositNeutralHadrons * process.pfMuIsoDepositPhotons * process.pfMuIsoDepositPU)

############################################
# Embed
if isEmbedded:
	process.origPfNeutralCands = cms.EDFilter("PdgIdPFCandidateSelector", src = cms.InputTag('particleFlow',"","RECO"), pdgId = cms.vint32(pfNeutralCandPdgIds))

	process.origKt6PFJetsNeutral = process.kt4PFJets.clone()
	process.origKt6PFJetsNeutral.rParam = cms.double(0.6)
	process.origKt6PFJetsNeutral.doAreaFastjet = cms.bool( True )
	process.origKt6PFJetsNeutral.doRhoFastjet = cms.bool( True )
	process.origKt6PFJetsNeutral.Rho_EtaMax = cms.double(2.5)
	process.origKt6PFJetsNeutral.src = cms.InputTag('origPfNeutralCands')

	process.rhoSequence *= process.origPfNeutralCands
	process.rhoSequence *= process.origKt6PFJetsNeutral

	triggerSource = cms.InputTag("TriggerResults", "", "EmbeddedRECO")
	triggerEventSource = cms.InputTag("hltTriggerSummaryAOD", "", "EmbeddedRECO")
	origTriggerSource = cms.InputTag("TriggerResults", "", "HLT")
        origVertexSource = cms.InputTag("offlinePrimaryVertices", "", "RECO")

	origCaloMetSource = cms.InputTag("metNoHF", "" , "RECO")
	origPfMetSource = cms.InputTag("pfMet", "", "RECO")

	origMuonSource = cms.InputTag("muons","","RECO")
	origRhoNeutralSource = cms.InputTag("origKt6PFJetsNeutral","rho","")

	origMuonChargedHadronIsolation = cms.InputTag("origPfmuIsoDepositChargedHadronCandidates")
	origMuonNeutralHadronIsolation = cms.InputTag("origPfmuIsoDepositNeutralHadronCandidates")
	origMuonPhotonIsolation = cms.InputTag("origPfmuIsoDepositPhotonCandidates")

	genParticlesSource = cms.InputTag("genParticles", "", "EmbeddedRECO")
	if not isData: origGenParticlesSource = cms.InputTag("genParticles", "", "SIM") # was: "HLT"
	else: origGenParticlesSource = cms.InputTag("")
else:
	triggerSource = cms.InputTag("TriggerResults", "", "HLT")
	triggerEventSource = cms.InputTag("hltTriggerSummaryAOD", "", "HLT")
	origTriggerSource = cms.InputTag("")
        origVertexSource = cms.InputTag("")

	origCaloMetSource = cms.InputTag("")
	origPfMetSource = cms.InputTag("")

	origMuonSource = cms.InputTag("")
	origRhoNeutralSource = cms.InputTag("")

	origMuonChargedHadronIsolation = cms.InputTag("")
	origMuonNeutralHadronIsolation = cms.InputTag("")
	origMuonPhotonIsolation = cms.InputTag("")

	if not isData: genParticlesSource = cms.InputTag("genParticles", "", "SIM") # was: "HLT"
	else: genParticlesSource = cms.InputTag("")
	origGenParticlesSource = cms.InputTag("")

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

############################################
# Rochester muon momentum corrections
#process.correctedMuons = cms.EDProducer('RochesterCorrMuonProducer',
#	src = cms.InputTag('muons'),
#	isMC = cms.bool(not isData))

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

process.diMuons = cms.EDAnalyzer('MuMuNTupleProducer',
   writeAllEvents = cms.bool(writeAllEvents),
   isEmbedded = cms.bool(isEmbedded),
   isRHEmbedded = cms.bool(isRHEmbedded),
   isGenEmbedded = cms.bool(isGenEmbedded),
   isData = cms.bool(isData),
   isSS = cms.bool(False),

   TriggerSource = triggerSource,
   TriggerEventSource = triggerEventSource,
   OrigTriggerSource = origTriggerSource,

   MuonIsoDepsChargedParticlesSource = cms.InputTag("pfMuIsoDepositChargedParticles"),
   MuonIsoDepsChargedHadronsSource = cms.InputTag("pfMuIsoDepositChargedHadrons"),
   MuonIsoDepsNeutralHadronsSource = cms.InputTag("pfMuIsoDepositNeutralHadrons"),
   MuonIsoDepsPhotonsSource = cms.InputTag("pfMuIsoDepositPhotons"),
   MuonIsoDepsPUSource = cms.InputTag("pfMuIsoDepositPU"),

   # Use this for study of mumu input events:
   #MuonSource = cms.InputTag("correctedMuons"),
   # Use this for study of embedded mumu->mumu samples:
   MuonSource = cms.InputTag("muons"),
   JetSource = jetSource,
   BTagSource = cms.InputTag("combinedSecondaryVertexBJetTags"),
   TrackSource = cms.InputTag("generalTracks"),
   OrigMuonSource = origMuonSource,
   RhoNeutralSource = cms.InputTag("kt6PFJetsNeutral","rho",""),
   OrigRhoNeutralSource = origRhoNeutralSource,
   CaloMetSource = cms.InputTag("metNoHF"),
   PfMetSource = cms.InputTag("pfMet"),
   PfType1CorrectedMetSource = cms.InputTag("pfType1CorrectedMet"),
   MvaPfMetSource = mvaPfMetSource,
   OrigCaloMetSource = origCaloMetSource,
   OrigPfMetSource = origPfMetSource,
   VertexSource = cms.InputTag("offlinePrimaryVertices"),
   PileupSummaryInfoSource = cms.InputTag("addPileupInfo", "", "HLT"),
   OrigVertexSource = origVertexSource,

   GenParticlesSource = genParticlesSource,
   OrigGenParticlesSource = origGenParticlesSource
)

process.diMuonsSS = process.diMuons.clone(
   isSS = cms.bool(True)
)

process.p = cms.Path(process.ak5PFJetsL1FastL2L3 * process.rhoSequence * process.ak5JetTracksAssociatorAtVertex * process.btagging * process.isoSequence * process.producePFMETCorrections * process.diMuons * process.diMuonsSS)
if isData:
	process.p.replace(process.ak5PFJetsL1FastL2L3, process.ak5PFJetsL1FastL2L3Residual)
if useMvaMet:
	process.p.replace(process.producePFMETCorrections, process.producePFMETCorrections * process.pfMEtMVAsequence)
if isEmbedded and isRHEmbedded:
	process.p.replace(process.diMuons, process.muonRadiationFilter2ResultStripSel1 * process.muonRadiationFilter2ResultStripSel2 * process.muonRadiationFilter2ResultStripSel3 * process.diMuons)
