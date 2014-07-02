#!/usr/bin/env python

# This script can be used to regenerate the ZmumuEvtSelEffCorrWeights.root file which might be necessary
# when the muon ID or pile-up conditions change. It should be run on a file created by
# AnalyzeEmbedMuMuNTuples run on a normal (non-embedded) Z->mumu Monte Carlo with PDGfilter set to 13
# (to skip tautau->mumu decays).

# There must be two iterations. In the first iteration, the pt vs. pt histogram is generated. The created
# file should then be used to run AnalyzeEmbedMuMuNTuples again, which then applies the pt-dependent
# corrections. Afterwards, this script should be run again, and it computes the eta-dependent corrections
# on top of the pt ones. The output of this second iteration of the script is the final and can be used
# for to correct embedded samples.

import sys
sys.argv.append('-b')

import ROOT

sourceFile = ROOT.TFile('out-embedtest_Summer12_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-v1-filter13.root')

def get(file, name):
	histo = file.Get('%s' % name)
	if not histo: return None
	return histo

def divide(numerator, denominator):
	histo = numerator.Clone('%s_ratio' % numerator.GetName())
	histo.Divide(denominator)
	return histo

# Level 1:
ptGen = get(sourceFile, 'zmumuGen/genPt')
ptRec = get(sourceFile, 'zmumuRec/recPt')

ratioPt = divide(ptRec, ptGen)

# Level 2 (might not be available when doing only the first iteration).
etaGen = get(sourceFile, 'zmumuGen/genEta')
etaRec = get(sourceFile, 'zmumuRec_PtCorr/recEta')

if etaRec: ratioEta = divide(etaRec, etaGen)

# The same for the pt 20/10 and eta 2.4 case:
ptGenEmb = get(sourceFile, 'zmumuGenEmb/genPt')
ptRecEmb = get(sourceFile, 'zmumuRecEmb/recPt')

ratioPtEmb = divide(ptRecEmb, ptGenEmb)
etaGenEmb = get(sourceFile, 'zmumuGenEmb/genEta')
etaRecEmb = get(sourceFile, 'zmumuRecEmb_PtCorr/recEta')
if etaRecEmb: ratioEtaEmb = divide(etaRecEmb, etaGenEmb)

outFile = ROOT.TFile("ZmumuEvtSelEffCorrWeights.root", "RECREATE")
ratioPt.Write("ZmumuEvtSelEff_muMinusPt_vs_muPlusPt")
if etaRec: ratioEta.Write("ZmumuEvtSelEffCorr_muMinusEta_vs_muPlusEta")
ratioPtEmb.Write("EMB_ZmumuEvtSelEff_muMinusPt_vs_muPlusPt")
if etaRecEmb: ratioEtaEmb.Write("EMB_ZmumuEvtSelEffCorr_muMinusEta_vs_muPlusEta")
outFile.Close()
