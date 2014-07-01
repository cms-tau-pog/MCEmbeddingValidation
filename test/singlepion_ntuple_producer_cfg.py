import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
	'/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/00037C53-AAD1-E111-B1BE-003048D45F38.root'
	#'/store/mc/Summer12_DR53X/GluGluToHToTauTau_M-125_8TeV-powheg-pythia6/AODSIM/PU_S10_START53_V7A-v1/0000/00E903E2-9FE9-E111-8B1E-003048FF86CA.root'
    )
)
process.TFileService = cms.Service("TFileService", fileName = cms.string('output.root'))

isData = False
isEmbedded = False
isGenEmbedded = False
if not "@IS_EMBEDDED@".startswith('@'):
	isEmbedded = ("@IS_EMBEDDED@".lower() == 'true')
if not "@IS_GEN_EMBEDDED@".startswith('@'):
	isGenEmbedded = ("@IS_GEN_EMBEDDED@".lower() == 'true')

## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

if isData:
	process.GlobalTag.globaltag = cms.string('FT_53_V21_AN4::All')
else:
	process.GlobalTag.globaltag = cms.string('START53_V23::All')

if isEmbedded:
	genParticlesSource = cms.InputTag("genParticles", "", "EmbeddedRECO")
	origGenParticlesSource = cms.InputTag("genParticles", "", "SIM")
else:
	genParticlesSource = cms.InputTag("genParticles", "", "SIM") # was: "HLT"
	origGenParticlesSource = cms.InputTag("")

############################################
# Radiation filter
process.load("Validation.MCEmbedding.muonRadiationFilter2_cfi")

if isGenEmbedded:
	process.muonRadiationFilter2.srcSelectedMuons = cms.InputTag('genMuonsFromZs', '', 'EmbeddedRECO')
elif isEmbedded:
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

process.spinTest = cms.EDAnalyzer('SinglePionNTupleProducer',
   isEmbedded = cms.bool(isEmbedded),
   isData = cms.bool(isData),

   PileupSummaryInfoSource = cms.InputTag("addPileupInfo", "", "HLT"),
   GenParticlesSource = genParticlesSource,
   OrigGenParticlesSource = origGenParticlesSource
)

if isEmbedded:
	process.p = cms.Path(process.muonRadiationFilter2ResultStripSel1 * process.muonRadiationFilter2ResultStripSel2 * process.muonRadiationFilter2ResultStripSel3 * process.spinTest)
else:
	process.p = cms.Path(process.spinTest)
