#!/usr/bin/env python

import sys
import math
import plot
import ROOT

mc = plot.DatasetEmbTest('MC', 'out-embedtest_mutau_Summer12_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-v1-mutau')
emb = plot.DatasetEmbTest('RH embedding', 'out-embedtest_mutau_Summer12_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-v1_RHEmbedded_v5_NoRotation-mutau')
embPf = plot.DatasetEmbTest('PF embedding', 'out-embedtest_mutau_Summer12_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-v1_PFEmbedded_v1-mutau')

MUON_PT_BINNING = [17., 18., 20., 22., 24., 26., 28., 30., 35., 40., 45., 50., 60., 80., 100., 150.]
TAU_PT_BINNING = [20., 22., 24., 26., 28., 30., 35., 40., 45., 50., 60., 80., 100., 150.]
JET_PT_BINNING = [30., 35., 40., 45., 50., 55., 60., 65., 70., 80., 90., 100., 125., 150., 200., 300.]
MUON_ETA_BINNING = [-2.1, -1.9, -1.7, -1.5, -1.3, -1.1, -0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1]
TAU_ETA_BINNING = [-2.3, -2.1, -1.9, -1.7, -1.5, -1.3, -1.1, -0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3]
PHI_BINNING = 5
GEN_DITAU_MASS_BINNING = [20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 72., 74., 76., 78., 80., 82., 84., 86., 88., 90., 92., 94., 96., 98., 100., 102., 104., 106., 108., 110., 115., 120., 125., 130., 140., 150., 175., 200.]
SVFIT_MASS_BINNING = [40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100., 105., 110., 115., 120., 125., 130., 140., 150., 175., 200.]
GEN_DITAU_VISMASS_BINNING = [20., 30., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 100., 120., 150., 200.]
GEN_DITAU_PT_BINNING = [0., 5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 60., 70., 80., 90., 100., 125., 150.]
GEN_DITAU_ETA_BINNING = [-9.0, -8.0, -7.0, -6.0, -5.5, -5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0]
MET_BINNING = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 35.0, 40.0, 50.0, 60.0, 80.0, 100.0]
ISO_BINNING = [0.0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.20, 0.25, 0.30, 0.40, 0.50, 0.60, 0.80, 1.0]
NPV_BINNING = [0.5, 5.5, 8.5, 10.5, 12.5, 14.5, 16.5, 18.5, 20.5, 22.5, 24.5, 26.5, 28.5, 30.5, 35.5, 40.5, 45.5, 50.5]

PLOTS = {
	# Generator level
	'mutau-gen-pt-muon': {
		'xaxis': (17., 150., 'GeV'),
		'xlabel': 'Generator Muon p_{T}',
		'yaxis': (1e0, 1e5, True),
		'binning': MUON_PT_BINNING,
		'normalization': ('*', '*'),
		'ratio-range': (0.85, 1.15),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #tau_{#mu}^{}#tau_{h}^{}', 'inclusive', 'genVisMasscut/muonGenPt', ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', 'genVisMasscut/muonGenPt', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'genVisMasscut/muonGenPt', ROOT.kRed, 21, None],
		]
	}, 'mutau-gen-pt-tau': {
		'xaxis': (20., 150., 'GeV'),
		'xlabel': 'Generator Tau p_{T}',
		'yaxis': (1e0, 1e5, True),
		'binning': TAU_PT_BINNING,
		'normalization': ('*', '*'),
		'ratio-range': (0.85, 1.15),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #tau_{#mu}^{}#tau_{h}^{}', 'inclusive', 'genVisMasscut/tauGenPt', ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', 'genVisMasscut/tauGenPt', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'genVisMasscut/tauGenPt', ROOT.kRed, 21, None],
		]
	}, 'mutau-gen-eta-muon': {
		'xaxis': (-2.1, 2.1, ''),
		'xlabel': 'Generator Muon #eta',
		'yaxis': (0, 7e4, False),
		'binning': MUON_ETA_BINNING,
		'normalization': ('*', '*'),
		'ratio-range': (0.95, 1.05),
		'legend-pos': (0.4, 0.2),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #tau_{#mu}^{}#tau_{h}^{}', 'inclusive', 'genVisMasscut/muonGenEta', ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', 'genVisMasscut/muonGenEta', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'genVisMasscut/muonGenEta', ROOT.kRed, 21, None],
		]
	}, 'mutau-gen-eta-tau': {
		'xaxis': (-2.3, 2.3, ''),
		'xlabel': 'Generator Tau #eta',
		'yaxis': (0, 7e4, False),
		'binning': TAU_ETA_BINNING,
		'normalization': ('*', '*'),
		'ratio-range': (0.95, 1.05),
		'legend-pos': (0.4, 0.2),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #tau_{#mu}^{}#tau_{h}^{}', 'inclusive', 'genVisMasscut/tauGenEta', ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', 'genVisMasscut/tauGenEta', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'genVisMasscut/tauGenEta', ROOT.kRed, 21, None],
		]
	}, 'mutau-gen-phi-muon': {
		'xaxis': (-math.pi, math.pi, 'rad'),
		'xlabel': 'Generator Muon #phi',
		'yaxis': (0, 4e4, False),
		'binning': PHI_BINNING,
		'normalization': ('*', '*'),
		'ratio-range': (0.95, 1.05),
		'legend-pos': (0.4, 0.2),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #tau_{#mu}^{}#tau_{h}^{}', 'inclusive', 'genVisMasscut/muonGenPhi', ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', 'genVisMasscut/muonGenPhi', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'genVisMasscut/muonGenPhi', ROOT.kRed, 21, None],
		]
	}, 'mutau-gen-phi-tau': {
		'xaxis': (-math.pi, math.pi, 'rad'),
		'xlabel': 'Generator Tau #phi',
		'yaxis': (0, 4e4, False),
		'binning': PHI_BINNING,
		'normalization': ('*', '*'),
		'ratio-range': (0.95, 1.05),
		'legend-pos': (0.4, 0.2),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #tau_{#mu}^{}#tau_{h}^{}', 'inclusive', 'genVisMasscut/tauGenPhi', ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', 'genVisMasscut/tauGenPhi', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'genVisMasscut/tauGenPhi', ROOT.kRed, 21, None],
		]
	}, 'mutau-gen-npv': {
		'xaxis': (0.5, 50.5, ''),
		'yaxis': (1e0, 1e5, True),
		'binning': NPV_BINNING,
		'normalization': ('*', '*'),
		'ratio-range': (0.90, 1.10),
		'legend-pos': (0.4, 0.2),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #tau_{#mu}^{}#tau_{h}^{}', 'inclusive', 'genVisMasscut/nPV', ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', 'genVisMasscut/nPV', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'genVisMasscut/nPV', ROOT.kRed, 21, None],
		]
	}, 'mutau-gen-ditau-mass': {
		'xaxis': (50., 150., 'GeV'),
		'yaxis': (1e0, 1e5, True),
		'binning': GEN_DITAU_MASS_BINNING,
		'normalization': ('*', '*'),
		'ratio-range': (0.50, 2.00),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #tau_{#mu}^{}#tau_{h}^{}', 'inclusive', 'genVisMasscut/diTauGenFinalMass', ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', 'genVisMasscut/diTauGenFinalMass', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'genVisMasscut/diTauGenFinalMass', ROOT.kRed, 21, None],
		]
	}, 'mutau-gen-ditau-visible-mass': {
		'xaxis': (20., 150., 'GeV'),
		'yaxis': (1e0, 1e5, True),
		'binning': GEN_DITAU_VISMASS_BINNING,
		'normalization': ('*', '*'),
		'ratio-range': (0.90, 1.10),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #tau_{#mu}^{}#tau_{h}^{}', 'inclusive', 'genVisMasscut/diTauGenVisMass', ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', 'genVisMasscut/diTauGenVisMass', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'genVisMasscut/diTauGenVisMass', ROOT.kRed, 21, None],
		]
	}, 'mutau-gen-ditau-pt': {
		'xaxis': (0., 150., 'GeV'),
		'xlabel': 'Generator di-tau p_{T}',
		'yaxis': (1e1, 1e5, True),
		'binning': GEN_DITAU_PT_BINNING,
		'normalization': ('*', '*'),
		'ratio-range': (0.90, 1.10),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #tau_{#mu}^{}#tau_{h}^{}', 'inclusive', 'genVisMasscut/diTauGenPt', ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', 'genVisMasscut/diTauGenPt', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'genVisMasscut/diTauGenPt', ROOT.kRed, 21, None],
		]
	# Reco-level
	}, 'mutau-rec-pt-muon': {
		'xaxis': (17., 150., 'GeV'),
		'xlabel': 'Muon p_{T}',
		'yaxis': (1e0, 1e4, True),
		'binning': MUON_PT_BINNING,
		'normalization': ('*', '*'),
		'ratio-range': (0.8, 1.3),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #tau_{#mu}^{}#tau_{h}^{}', 'inclusive', 'selected/muonPt', ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', 'selected/muonPt', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'selected/muonPt', ROOT.kRed, 21, None],
		]
	}, 'mutau-rec-pt-tau': {
		'xaxis': (20., 150., 'GeV'),
		'xlabel': 'Tau p_{T}',
		'yaxis': (1e0, 1e4, True),
		'binning': TAU_PT_BINNING,
		'normalization': ('*', '*'),
		'ratio-range': (0.85, 1.15),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #tau_{#mu}^{}#tau_{h}^{}', 'inclusive', 'selected/tauPt', ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', 'selected/tauPt', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'selected/tauPt', ROOT.kRed, 21, None],
		]
	}, 'mutau-rec-eta-muon': {
		'xaxis': (-2.1, 2.1, ''),
		'xlabel': 'Muon #eta',
		'yaxis': (0, 3e4, False),
		'binning': MUON_ETA_BINNING,
		'normalization': ('*', '*'),
		'ratio-range': (0.85, 1.15),
		'legend-pos': (0.4, 0.2),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #tau_{#mu}^{}#tau_{h}^{}', 'inclusive', 'selected/muonEta', ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', 'selected/muonEta', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'selected/muonEta', ROOT.kRed, 21, None],
		]
	}, 'mutau-rec-eta-tau': {
		'xaxis': (-2.3, 2.3, ''),
		'xlabel': 'Tau #eta',
		'yaxis': (0, 3e4, False),
		'binning': TAU_ETA_BINNING,
		'normalization': ('*', '*'),
		'ratio-range': (0.85, 1.15),
		'legend-pos': (0.4, 0.2),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #tau_{#mu}^{}#tau_{h}^{}', 'inclusive', 'selected/tauEta', ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', 'selected/tauEta', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'selected/tauEta', ROOT.kRed, 21, None],
		]
	}, 'mutau-rec-phi-muon': {
		'xaxis': (-math.pi, math.pi, 'rad'),
		'xlabel': 'Muon #phi',
		'yaxis': (0, 2e4, False),
		'binning': PHI_BINNING,
		'normalization': ('*', '*'),
		'ratio-range': (0.90, 1.10),
		'legend-pos': (0.4, 0.2),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #tau_{#mu}^{}#tau_{h}^{}', 'inclusive', 'selected/muonPhi', ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', 'selected/muonPhi', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'selected/muonPhi', ROOT.kRed, 21, None],
		]
	}, 'mutau-rec-phi-tau': {
		'xaxis': (-math.pi, math.pi, 'rad'),
		'xlabel': 'Tau #phi',
		'yaxis': (0, 2e4, False),
		'binning': PHI_BINNING,
		'normalization': ('*', '*'),
		'ratio-range': (0.90, 1.10),
		'legend-pos': (0.4, 0.2),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #tau_{#mu}^{}#tau_{h}^{}', 'inclusive', 'selected/tauPhi', ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', 'selected/tauPhi', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'selected/tauPhi', ROOT.kRed, 21, None],
		]
	}, 'mutau-rec-npv': {
		'xaxis': (0.5, 50.5, ''),
		'yaxis': (1e0, 1e5, True),
		'binning': NPV_BINNING,
		'normalization': ('*', '*'),
		'ratio-range': (0.85, 1.15),
		'legend-pos': (0.4, 0.2),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #tau_{#mu}^{}#tau_{h}^{}', 'inclusive', 'selected/nPV', ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', 'selected/nPV', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'selected/nPV', ROOT.kRed, 21, None],
		]
	}, 'mutau-rec-ditau-visible-mass': {
		'xaxis': (20., 150., 'GeV'),
		'yaxis': (1e0, 1e4, True),
		'binning': GEN_DITAU_VISMASS_BINNING,
		'normalization': ('*', '*'),
		'ratio-range': (0.85, 1.15),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #tau_{#mu}^{}#tau_{h}^{}', 'inclusive', 'selected/visibleMass', ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', 'selected/visibleMass', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'selected/visibleMass', ROOT.kRed, 21, None],
		]
	}, 'mutau-rec-svfit-mass': {
		'xaxis': (40., 200., 'GeV'),
		'yaxis': (1e0, 1e4, True),
		'binning': SVFIT_MASS_BINNING,
		'normalization': ('*', '*'),
		'ratio-range': (0.85, 1.15),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #tau_{#mu}^{}#tau_{h}^{}', 'inclusive', 'selected/svfitMassHigh', ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', 'selected/svfitMassHigh', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'selected/svfitMassHigh', ROOT.kRed, 21, None],
		]
	}, 'mutau-rec-higgs-pt': {
		'xaxis': (0., 150., 'GeV'),
		'xlabel': 'Higgs p_{T}',
		'yaxis': (1e1, 1e4, True),
		'binning': GEN_DITAU_PT_BINNING,
		'normalization': ('*', '*'),
		'ratio-range': (0.85, 1.15),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #tau_{#mu}^{}#tau_{h}^{}', 'inclusive', 'selected/higgsPt', ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', 'selected/higgsPt', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'selected/higgsPt', ROOT.kRed, 21, None],
		]
	}, 'mutau-rec-njets30': {
		'xaxis': (-0.5, 10.5, ''),
		'yaxis': (1e-1, 1e5, True),
		'binning': [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 6.5, 10.5],
		'normalization': ('*', '*'),
		'ratio-range': (0.70, 1.10),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #tau_{#mu}^{}#tau_{h}^{}', 'inclusive', 'selected/nJets30', ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', 'selected/nJets30', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'selected/nJets30', ROOT.kRed, 21, None],
		]
	}, 'mutau-rec-leadjet-pt': {
		'xaxis': (30., 150., 'GeV'),
		'xlabel': 'Jet p_{T}',
		'yaxis': (1e1, 1e4, True),
		'binning': JET_PT_BINNING,
		'normalization': ('*', '*'),
		'ratio-range': (0.85, 1.15),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #tau_{#mu}^{}#tau_{h}^{}', 'inclusive', 'selected/leadJetPt', ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', 'selected/leadJetPt', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'selected/leadJetPt', ROOT.kRed, 21, None],
		]
	}, 'mutau-rec-pfmet': {
		'xaxis': (0., 100., 'GeV'),
		'yaxis': (1e1, 1e4, True),
		'binning': MET_BINNING,
		'normalization': ('*', '*'),
		'ratio-range': (0.90, 1.10),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #tau_{#mu}^{}#tau_{h}^{}', 'inclusive', 'selected/pfMet', ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', 'selected/pfMet', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'selected/pfMet', ROOT.kRed, 21, None],
		]
	}, 'mutau-rec-calomet': {
		'xaxis': (0., 100., 'GeV'),
		'yaxis': (1e1, 1e4, True),
		'binning': MET_BINNING,
		'normalization': ('*', '*'),
		'ratio-range': (0.80, 1.20),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #tau_{#mu}^{}#tau_{h}^{}', 'inclusive', 'selected/caloMet', ROOT.kBlack, 20, None],
#			[embPf, 'PF embedding', 'inclusive', 'selected/caloMet', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'selected/caloMet', ROOT.kRed, 21, None],
		]
	# Muon Isolation
	}, 'mutau-rec-muon-charged-hadron-iso': {
		'xaxis': (0., 1., ''),
		'yaxis': (1e1, 1e8, True),
		'binning': ISO_BINNING,
		'normalization': ('*', '*'),
		'ratio-range': (0.00, 1.10),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #tau_{#mu}^{}#tau_{h}^{}', 'inclusive', 'lep_identified/muonChargedHadronIso', ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', 'lep_identified/muonChargedHadronIso', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'lep_identified/muonChargedHadronIso', ROOT.kRed, 21, None],
		]
	}, 'mutau-rec-muon-neutral-hadron-iso': {
		'xaxis': (0., 1., ''),
		'yaxis': (1e1, 1e8, True),
		'binning': ISO_BINNING,
		'normalization': ('*', '*'),
		'ratio-range': (0.50, 2.00),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #tau_{#mu}^{}#tau_{h}^{}', 'inclusive', 'lep_identified/muonNeutralHadronIso', ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', 'lep_identified/muonNeutralHadronIso', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'lep_identified/muonNeutralHadronIso', ROOT.kRed, 21, None],
		]

	}, 'mutau-rec-muon-photon-iso': {
		'xaxis': (0., 1., ''),
		'yaxis': (1e2, 1e7, True),
		'binning': ISO_BINNING,
		'normalization': ('*', '*'),
		'ratio-range': (0.80, 3.00),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #tau_{#mu}^{}#tau_{h}^{}', 'inclusive', 'lep_identified/muonPhotonIso', ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', 'lep_identified/muonPhotonIso', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'lep_identified/muonPhotonIso', ROOT.kRed, 21, None],
		]

	}, 'mutau-rec-muon-pu-iso': {
		'xaxis': (0., 1., ''),
		'yaxis': (1e1, 1e8, True),
		'binning': ISO_BINNING,
		'normalization': ('*', '*'),
		'ratio-range': (0.90, 1.10),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #tau_{#mu}^{}#tau_{h}^{}', 'inclusive', 'lep_identified/muonPUIso', ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', 'lep_identified/muonPUIso', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'lep_identified/muonPUIso', ROOT.kRed, 21, None],
		]
	}, 'mutau-rec-muon-combined-iso': {
		'xaxis': (0., 1., ''),
		'yaxis': (1e2, 1e7, True),
		'binning': ISO_BINNING,
		'normalization': ('*', '*'),
		'ratio-range': (0.50, 1.60),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #tau_{#mu}^{}#tau_{h}^{}', 'inclusive', 'lep_identified/muonIso', ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', 'lep_identified/muonIso', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'lep_identified/muonIso', ROOT.kRed, 21, None],
		]
	# Efficiency plots
	}, 'mutau-eff-muonid-pt': {
		'xaxis': (17., 150., 'GeV'),
		'xlabel': 'Generator Muon p_{T}',
		'yaxis': (0.8, 1.0, False),
		'ylabel': 'Muon Identification Efficiency',
		'binning': MUON_PT_BINNING,
		'normalization': None,
		'ratio-range': (0.90, 1.10),
		'legend-pos': (0.2, 0.2),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #tau_{#mu}^{}#tau_{h}^{}', 'inclusive', ('genVisMasscut/muonGenPt', 'lep_identified/muonGenPt'), ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', ('genVisMasscut/muonGenPt', 'lep_identified/muonGenPt'), ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', ('genVisMasscut/muonGenPt', 'lep_identified/muonGenPt'), ROOT.kRed, 21, None],
		]
	}, 'mutau-eff-muonid-eta': {
		'xaxis': (-2.1, 2.1, ''),
		'xlabel': 'Generator Muon #eta',
		'yaxis': (0.75, 1.0, False),
		'ylabel': 'Muon Identification Efficiency',
		'binning': MUON_ETA_BINNING,
		'normalization': None,
		'ratio-range': (0.95, 1.05),
		'legend-pos': (0.2, 0.2),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #tau_{#mu}^{}#tau_{h}^{}', 'inclusive', ('genVisMasscut/muonGenEta', 'lep_identified/muonGenEta'), ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', ('genVisMasscut/muonGenEta', 'lep_identified/muonGenEta'), ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', ('genVisMasscut/muonGenEta', 'lep_identified/muonGenEta'), ROOT.kRed, 21, None],
		]
	}, 'mutau-eff-muoniso-pt': {
		'xaxis': (17., 150., 'GeV'),
		'xlabel': 'Generator Muon p_{T}',
		'yaxis': (0.7, 1.0, False),
		'ylabel': 'Muon Isolation Efficiency',
		'binning': MUON_PT_BINNING,
		'normalization': None,
		'ratio-range': (0.95, 1.05),
		'legend-pos': (0.5, 0.2),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #tau_{#mu}^{}#tau_{h}^{}', 'inclusive', ('lep_identified/muonGenPt', 'lep_isolated/muonGenPt'), ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', ('lep_identified/muonGenPt', 'lep_isolated/muonGenPt'), ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', ('lep_identified/muonGenPt', 'lep_isolated/muonGenPt'), ROOT.kRed, 21, None],
		]
	}, 'mutau-eff-muoniso-eta': {
		'xaxis': (-2.1, 2.1, ''),
		'xlabel': 'Generator Muon #eta',
		'yaxis': (0.7, 1.0, False),
		'ylabel': 'Muon Isolation Efficiency',
		'binning': MUON_ETA_BINNING,
		'normalization': None,
		'ratio-range': (0.95, 1.05),
		'legend-pos': (0.4, 0.6),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #tau_{#mu}^{}#tau_{h}^{}', 'inclusive', ('lep_identified/muonGenEta', 'lep_isolated/muonGenEta'), ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', ('lep_identified/muonGenEta', 'lep_isolated/muonGenEta'), ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', ('lep_identified/muonGenEta', 'lep_isolated/muonGenEta'), ROOT.kRed, 21, None],
		]
	}, 'mutau-eff-tauid-pt': {
		'xaxis': (20., 150., 'GeV'),
		'xlabel': 'Generator Tau p_{T}',
		'yaxis': (0.0, 0.6, False),
		'ylabel': 'Tau Identification Efficiency',
		'binning': TAU_PT_BINNING,
		'normalization': None,
		'ratio-range': (0.9, 1.2),
		'legend-pos': (0.2, 0.2),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #tau_{#mu}^{}#tau_{h}^{}', 'inclusive', ('genVisMasscut/tauGenPt', 'tau_identified/tauGenPt'), ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', ('genVisMasscut/tauGenPt', 'tau_identified/tauGenPt'), ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', ('genVisMasscut/tauGenPt', 'tau_identified/tauGenPt'), ROOT.kRed, 21, None],
		]
	}, 'mutau-eff-tauid-eta': {
		'xaxis': (-2.3, 2.3, ''),
		'xlabel': 'Generator Tau #eta',
		'yaxis': (0.0, 0.6, False),
		'ylabel': 'Tau Identification Efficiency',
		'binning': TAU_ETA_BINNING,
		'normalization': None,
		'ratio-range': (0.9, 1.2),
		'legend-pos': (0.2, 0.2),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #tau_{#mu}^{}#tau_{h}^{}', 'inclusive', ('genVisMasscut/tauGenEta', 'tau_identified/tauGenEta'), ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', ('genVisMasscut/tauGenEta', 'tau_identified/tauGenEta'), ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', ('genVisMasscut/tauGenEta', 'tau_identified/tauGenEta'), ROOT.kRed, 21, None],
		]

	}, 'mutau-eff-rec-npv': {
		'xaxis': (0.5, 50.5, ''),
		'yaxis': (0.0, 0.5, False),
		'ylabel': 'Reconstruction Efficiency',
		'binning': NPV_BINNING,
		'normalization': None,
		'ratio-range': (0.90, 1.30),
		'legend-pos': (0.4, 0.2),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #tau_{#mu}^{}#tau_{h}^{}', 'inclusive', ('genVisMasscut/nPV', 'selected/nPV'), ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', ('genVisMasscut/nPV', 'selected/nPV'), ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', ('genVisMasscut/nPV', 'selected/nPV'), ROOT.kRed, 21, None],
		]
	}
}

plot.create_plots(PLOTS)
