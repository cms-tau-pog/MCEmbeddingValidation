#!/usr/bin/env python

import sys
import plot
import ROOT

refZ = plot.DatasetEmbTest('Z', 'ReferenceSinglePionHistograms/DY')
refH = plot.DatasetEmbTest('H', 'ReferenceSinglePionHistograms/GluGluH125')

X1_BINNING = [0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.60, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0]
ZS_BINNING = [-0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.15, -0.10, -0.05, 0., 0.05, 0.1, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]
VISMASS_BINNING = [0., 5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100.]
#PT_BINNING = [0., 5., 10., 14, 17., 20., 23., 26., 29., 32., 34., 36., 38., 40., 41., 42., 43., 44., 45., 46., 47., 48., 49., 50., 52., 54., 56., 60., 70., 80., 100., 150.]
PT_BINNING = [0., 5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 60., 70., 80., 100., 150.]

PLOTS = {
	# z_s variable
	'zs': {
		'xaxis': (-0.5, 0.5, ''),
		'yaxis': (0, 6e4, False),
		'binning': ZS_BINNING,
		'normalization': ('*', '*'),
		'ratio-range': (0.7, 1.30),
		'legend-pos': (0.4, 0.2),
		'samples': [
			[refZ, 'Z/#gamma_{}^{}* #rightarrow #tau_{h}^{}#tau_{h}^{}', 'inclusiveNoRadiation', 'all/zs_rf', ROOT.kBlack, 20, None],
			[refH, 'H_{}^{} #rightarrow #tau_{h}^{}#tau_{h}^{}', 'inclusiveNoRadiation', 'all/zs_rf', ROOT.kRed, 21, None]
		]
	# visible mass
	}, 'vismass': {
		'xaxis': (0., 100., 'GeV'),
		'yaxis': (0, 7e2, False),
		'binning': VISMASS_BINNING,
		'normalization': ('*', '*'),
		'ratio-range': (0.7, 1.30),
		'legend-pos': (0.25, 0.2),
		'samples': [
			[refZ, 'Z/#gamma_{}^{}* #rightarrow #tau_{h}^{}#tau_{h}^{}', 'inclusiveNoRadiation', 'all/vismass', ROOT.kBlack, 20, None],
			[refH, 'H_{}^{} #rightarrow #tau_{h}^{}#tau_{h}^{}', 'inclusiveNoRadiation', 'all/vismass', ROOT.kRed, 21, None]
		]
	}
	# more plots can be added here...
}

plot.create_plots(PLOTS)
