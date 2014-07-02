#!/usr/bin/env python

import sys
import math
import plot
import ROOT

mc = plot.DatasetEmbTest('MC', 'out-embedtest_Summer12_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-v1-filter13')
emb = plot.DatasetEmbTest('RH embedding', 'out-embedtest_Summer12_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-v1_RHEmbedded_v5_NoRotation-filter13')
embPf = plot.DatasetEmbTest('PF embedding', 'out-embedtest_Summer12_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-v1_PFEmbedded_v1-filter13')

PT_BINNING = [10., 14, 17., 20., 23., 26., 29., 32., 34., 36., 38., 40., 41., 42., 43., 44., 45., 46., 47., 48., 49., 50., 52., 54., 56., 60., 70., 80., 100., 150.]
JET_PT_BINNING = [30., 35., 40., 45., 50., 55., 60., 65., 70., 80., 90., 100., 125., 150., 200., 300.]
GEN_ETA_BINNING = [-2.4, -2.3, -2.1, -1.9, -1.7, -1.5, -1.3, -1.1, -0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.4]
REC_ETA_BINNING = [-2.1, -1.9, -1.7, -1.5, -1.3, -1.1, -0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1]
DIMUON_MASS_BINNING = [20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 72., 74., 76., 78., 80., 82., 84., 85., 86., 87., 88., 89., 90., 91., 92., 93., 94., 95., 96., 98., 100., 102., 104., 106., 108., 110., 115., 120., 125., 130., 140., 150., 175., 200.]
DIMUON_PT_BINNING = [0., 5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 60., 70., 80., 90., 100., 125., 150.]
DIMUON_ETA_BINNING = [-9.0, -8.0, -7.0, -6.0, -5.5, -5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0]
MET_BINNING = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 35.0, 40.0, 50.0, 60.0, 80.0, 100.0]
NPV_BINNING = [0.5, 5.5, 8.5, 10.5, 12.5, 14.5, 16.5, 18.5, 20.5, 22.5, 24.5, 26.5, 28.5, 30.5, 35.5, 40.5, 45.5, 50.5]

PLOTS = {
	'gen-pt-pos': {
		'xaxis': (10, 100, 'GeV'),
		'xlabel': 'Generator Muon p_{T}',
		'yaxis': (1e3, 1e6, True),
		'binning': PT_BINNING,
		'normalization': ('*', '*'),
		'ratio-range': (0.85, 1.15),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #mu#mu', 'inclusive', 'genFinalMasscut/posGenMuonPt', ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', 'genFinalMasscut/posGenMuonPt', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'genFinalMasscut/posGenMuonPt', ROOT.kRed, 21, None],
		]
	}, 'gen-eta-pos': {
		'xaxis': (-2.4, 2.4, ''),
		'xlabel': 'Generator Muon #eta',
		'yaxis': (0, 3e6, False),
		'binning': GEN_ETA_BINNING,
		'normalization': ('*', '*'),
		'ratio-range': (0.95, 1.05),
		'legend-pos': (0.4, 0.2),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #mu#mu', 'inclusive', 'genFinalMasscut/posGenMuonEta', ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', 'genFinalMasscut/posGenMuonEta', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'genFinalMasscut/posGenMuonEta', ROOT.kRed, 21, None],
		]
	}, 'gen-phi-pos': {
		'xaxis': (-math.pi, math.pi, 'rad'),
		'xlabel': 'Generator Muon #phi',
		'yaxis': (0, 2e6, False),
		'binning': 5,
		'ratio-range': (0.95, 1.05),
		'normalization': ('*', '*'),
		'legend-pos': (0.4, 0.2),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #mu#mu', 'inclusive', 'genFinalMasscut/posGenMuonPhi', ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', 'genFinalMasscut/posGenMuonPhi', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'genFinalMasscut/posGenMuonPhi', ROOT.kRed, 21, None],
		]
	}, 'gen-dimuon-mass': {
		'xaxis': (60., 120., 'GeV'),
		'yaxis': (1e3, 1e7, True),
		'binning': DIMUON_MASS_BINNING,
		'normalization': ('*', '*'),
		'ratio-range': (0.6, 1.8),
		'legend-pos': (0.17, 0.70),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #mu#mu', 'inclusive', 'genFinalMasscut/genFinalDiMuonMass', ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', 'genFinalMasscut/genFinalDiMuonMass', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'genFinalMasscut/genFinalDiMuonMass', ROOT.kRed, 21, None],
		]
	}, 'gen-dimuon-pt': {
		'xaxis': (0., 120., 'GeV'),
		'xlabel': 'Generator Di-muon p_{T}',
		'yaxis': (1e3, 1e6, True),
		'binning': DIMUON_PT_BINNING,
		'normalization': ('*', '*'),
		'ratio-range': (0.85, 1.15),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #mu#mu', 'inclusive', 'genFinalMasscut/genDiMuonPt', ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', 'genFinalMasscut/genDiMuonPt', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'genFinalMasscut/genDiMuonPt', ROOT.kRed, 21, None],
		]
	}, 'gen-dimuon-eta': {
		'xaxis': (-9.0, 9.0, ''),
		'xlabel': 'Generator Di-muon #eta',
		'yaxis': (1e1, 1e6, True),
		'binning': DIMUON_ETA_BINNING,
		'normalization': ('*', '*'),
		'legend-pos': (0.4, 0.2),
		'ratio-range': (0.6, 1.4),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #mu#mu', 'inclusive', 'genFinalMasscut/genDiMuonEta', ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', 'genFinalMasscut/genDiMuonEta', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'genFinalMasscut/genDiMuonEta', ROOT.kRed, 21, None],
		]
	}, 'rec-pt-pos': {
		'xaxis': (10, 100, 'GeV'),
		'xlabel': 'Muon p_{T}',
		'yaxis': (1e2, 1e6, True),
		'binning': PT_BINNING,
		'normalization': ('*', '*'),
		'ratio-range': (0.9, 1.3),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #mu#mu', 'inclusive', 'recMasscut/posMuonPt', ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', 'recMasscut/posMuonPt', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'recMasscut/posMuonPt', ROOT.kRed, 21, None],
		]
	}, 'rec-pt-pos-no-iso': {
		'xaxis': (10, 100, 'GeV'),
		'xlabel': 'Muon p_{T}',
		'yaxis': (1e2, 1e6, True),
		'binning': PT_BINNING,
		'normalization': ('*', '*'),
		'ratio-range': (0.9, 1.05),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #mu#mu', 'inclusive', 'recMasscutNoIsolation/posMuonPt', ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', 'recMasscutNoIsolation/posMuonPt', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'recMasscutNoIsolation/posMuonPt', ROOT.kRed, 21, None],
		]
	}, 'rec-eta-pos': {
		'xaxis': (-2.1, 2.1, ''),
		'xlabel': 'Muon #eta',
		'yaxis': (0, 1e6, False),
		'binning': REC_ETA_BINNING,
		'ratio-range': (0.93, 1.07),
		'legend-pos': (0.4, 0.2),
		'normalization': ('*', '*'),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #mu#mu', 'inclusive', 'recMasscut/posMuonEta', ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', 'recMasscut/posMuonEta', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'recMasscut/posMuonEta', ROOT.kRed, 21, None],
		]
	}, 'rec-phi-pos': {
		'xaxis': (-math.pi, math.pi, 'rad'),
		'xlabel': 'Muon #phi',
		'yaxis': (0, 6e5, False),
		'binning': 5,
		'ratio-range': (0.98, 1.02),
		'legend-pos': (0.60, 0.20),
		'normalization': ('*', '*'),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #mu#mu', 'inclusive', 'recMasscut/posMuonPhi', ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', 'recMasscut/posMuonPhi', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'recMasscut/posMuonPhi', ROOT.kRed, 21, None],
		]
	}, 'rec-dimuon-mass': {
		'xaxis': (60., 120., 'GeV'),
		'yaxis': (1e3, 1e6, True),
		'binning': DIMUON_MASS_BINNING,
		'ratio-range': (0.8, 1.3),
		'legend-pos': (0.17, 0.70),
		'normalization': ('*', '*'),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #mu#mu', 'inclusive', 'recMasscut/diMuonMass', ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', 'recMasscut/diMuonMass', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'recMasscut/diMuonMass', ROOT.kRed, 21, None],
		]
	}, 'rec-dimuon-pt': {
		'xaxis': (0., 120., 'GeV'),
		'xlabel': 'Di-muon p_{T}',
		'yaxis': (1e2, 1e6, True),
		'ratio-range': (0.85, 1.15),
		'binning': DIMUON_PT_BINNING,
		'normalization': ('*', '*'),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #mu#mu', 'inclusive', 'recMasscut/diMuonPt', ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', 'recMasscut/diMuonPt', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'recMasscut/diMuonPt', ROOT.kRed, 21, None],
		]
	}, 'rec-dimuon-eta': {
		'xaxis': (-9.0, 9.0, ''),
		'xlabel': 'Di-muon #eta',
		'yaxis': (1e1, 1e6, True),
		'binning': DIMUON_ETA_BINNING,
		'legend-pos': (0.4, 0.2),
		'normalization': ('*', '*'),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #mu#mu', 'inclusive', 'recMasscut/diMuonEta', ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', 'recMasscut/diMuonEta', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'recMasscut/diMuonEta', ROOT.kRed, 21, None],
		]
	}, 'gen-npv': {
		'xaxis': (0.5, 49.5, ''),
		'yaxis': (1e1, 1e6, True),
		'binning': NPV_BINNING,
		'ratio-range': (0.8, 1.2),
		'legend-pos': (0.3, 0.2),
		'normalization': ('*', '*'),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #mu#mu', 'inclusive', 'genFinalMasscut/nPV', ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', 'genFinalMasscut/nPV', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'genFinalMasscut/nPV', ROOT.kRed, 21, None],
		]
	}, 'rec-npv-no-iso': {
		'xaxis': (0.5, 49.5, ''),
		'yaxis': (1e1, 1e6, True),
		'binning': NPV_BINNING,
		'ratio-range': (0.8, 1.2),
		'legend-pos': (0.3, 0.2),
		'normalization': ('*', '*'),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #mu#mu', 'inclusive', 'recMasscutNoIsolation/nPV', ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', 'recMasscutNoIsolation/nPV', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'recMasscutNoIsolation/nPV', ROOT.kRed, 21, None],
		]
	}, 'rec-npv': {
		'xaxis': (0.5, 49.5, ''),
		'yaxis': (1e1, 1e6, True),
		'binning': NPV_BINNING,
		'ratio-range': (0.8, 1.2),
		'legend-pos': (0.3, 0.2),
		'normalization': ('*', '*'),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #mu#mu', 'inclusive', 'recMasscut/nPV', ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', 'recMasscut/nPV', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'recMasscut/nPV', ROOT.kRed, 21, None],
		]
	}, 'rec-njets30': {
		'xaxis': (-0.5, 10.5, ''),
		'yaxis': (1., 1e7, True),
		'binning': [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 6.5, 10.5],
		'normalization': ('*', '*'),
		'ratio-range': (0.6, 1.5),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #mu#mu', 'inclusive', 'recMasscut/nJets30', ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', 'recMasscut/nJets30', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'recMasscut/nJets30', ROOT.kRed, 21, None],
		]
	}, 'rec-leadjet-pt': {
		'xaxis': (30., 300., ''),
		'xlabel': 'Jet p_{T}',
		'yaxis': (1., 1e5, True),
		'binning': JET_PT_BINNING,
		'ratio-range': (0.85, 1.15),
		'normalization': ('*', '*'),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #mu#mu', 'inclusive', 'recMasscut/leadJetPt', ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', 'recMasscut/leadJetPt', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'recMasscut/leadJetPt', ROOT.kRed, 21, None],
		]
	}, 'rec-pfmet': {
		'xaxis': (0., 100.0, 'GeV'),
		'yaxis': (1e1, 1e6, True),
		'binning': MET_BINNING,
		'ratio-range': (0.9, 1.3),
		'normalization': ('*', '*'),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #mu#mu', 'inclusive', 'recMasscut/pfMet', ROOT.kBlack, 20, None],
			[embPf, 'PF embedding', 'inclusive', 'recMasscut/pfMet', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'recMasscut/pfMet', ROOT.kRed, 21, None],
		]
	}, 'rec-calomet': {
		'xaxis': (0., 100.0, 'GeV'),
		'yaxis': (1e2, 1e6, True),
		'binning': MET_BINNING,
		'ratio-range': (0.94, 1.06),
		'normalization': ('*', '*'),
		'samples': [
			[mc, 'MC Z/#gamma_{}^{}* #rightarrow #mu#mu', 'inclusive', 'recMasscut/caloMet', ROOT.kBlack, 20, None],
#			[embPf, 'PF embedding', 'inclusive', 'recMasscut/caloMet', ROOT.kGreen+3, 23, None],
			[emb, 'RH embedding', 'inclusive', 'recMasscut/caloMet', ROOT.kRed, 21, None],
		]
	}
}

plot.create_plots(PLOTS)
