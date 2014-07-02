#!/usr/bin/env python

import sys
sys.argv.append('-b')

import copy
import math
import ROOT
ROOT.gStyle.SetLineWidth(2)
ROOT.gStyle.SetCanvasBorderSize(2)
ROOT.gStyle.SetPadBorderSize(2)
ROOT.gStyle.SetHistLineWidth(2)
ROOT.gStyle.SetFrameLineWidth(2)

ROOT.gROOT.ProcessLineSync('.x tdrStyle.C')
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)

import numpy
#import stattest

def SingleIntegral(histo):
	# Include underflow and overflow bins
	sum = sumw2 = 0
	for x in range(0, histo.GetNbinsX()+2):
		sum += histo.GetBinContent(x)
		sumw2 += histo.GetBinError(x)*histo.GetBinError(x)
	return sum, sumw2

class Dataset:
	def __init__(self, label):
		self.label = label

	def getLuminosity(self):
		return 19712. # Determined w/ pixelLumiCalc

	def getLabel(self):
		return self.label

	def getHistogram(self, histogram):
		return self.getHistogramWithConfiguration(histogram, '')

class DatasetFile(Dataset):
	def __init__(self, label, filename):
		Dataset.__init__(self, label)
		self.filename = filename
		self.files = {}

	def getHistogramWithConfiguration(self, histogram, configuration):
		if configuration not in self.files:
			configuration_string = configuration
			if configuration_string != '':
				configuration_string = '-' + configuration_string

			if type(self.filename) == list:
				self.files[configuration] = [ROOT.TFile('%s%s.root' % (filename, configuration_string), 'READ') for filename in self.filename]
			else:
				self.files[configuration] = [ROOT.TFile('%s%s.root' % (self.filename, configuration_string), 'READ')]

		histo = None
		for file in self.files[configuration]:
			cur_histo = file.Get(histogram).Clone()
			if histo == None:
				histo = cur_histo
			else:
				histo.Add(cur_histo)

		histo.GetXaxis().SetLabelOffset(0.005)
		histo.GetXaxis().SetTitleOffset(1.0)
		histo.GetXaxis().SetTitleFont(62)
		histo.GetXaxis().SetTitleSize(0.04)
		histo.GetXaxis().SetLabelFont(62)
		histo.GetXaxis().SetLabelSize(0.04)
		histo.GetYaxis().SetTitleFont(62)
		histo.GetYaxis().SetTitleSize(0.04)
		histo.GetYaxis().SetLabelFont(62)
		histo.GetYaxis().SetLabelSize(0.04)

		return histo

class DatasetMC(DatasetFile):
	def __init__(self, label, filename, color, xsection, scale = 1.0):
		DatasetFile.__init__(self, label, filename)
		self.color = color
		self.xsection = xsection
		self.nevents = -1
		self.scale = scale

	def getHistogramWithConfiguration(self, histogram, configuration):
		if self.nevents == -1:
			entries = DatasetFile.getHistogramWithConfiguration(self, 'nEvents', configuration)
			self.nevents = entries.GetBinContent(1)

		histo = DatasetFile.getHistogramWithConfiguration(self, histogram, configuration)
		histo.Scale(self.scale * self.getLuminosity() * self.xsection / self.nevents)
		histo.SetFillColor(self.color)
		histo.SetLineColor(ROOT.kBlack)
		histo.SetLineWidth(2)
		return histo

	def getSystUncertainties(self):
		return {
			'lumi': 0.026,
			'theory': 0.043,
		}

class DatasetMCCombined(Dataset):
	def __init__(self, label, color, dataset_list):
		Dataset.__init__(self, label)
		self.dataset_list = dataset_list
		self.color = color

	def getHistogramWithConfiguration(self, histogram, configuration):
		histo = self.dataset_list[0].getHistogramWithConfiguration(histogram, configuration)
		for dataset in self.dataset_list[1:]:
			histo.Add(dataset.getHistogramWithConfiguration(histogram, configuration))
		histo.SetFillColor(self.color)
		return histo.Clone()

	def getSystUncertainties(self):
		return {
			'lumi': 0.026,
			'theory': max([x.getSystUncertainties()['theory'] for x in self.dataset_list]),
		}

class DatasetQCD(Dataset):
	def __init__(self, label, color, data, bkg):
		Dataset.__init__(self, label)

		self.color = color

		self.data = data
		self.bkg = bkg
		self.ss_os_factor = -1.

	def getHistogram(self, histogram):
		if self.ss_os_factor < 0:
			os_events = SingleIntegral(self.data.getHistogramWithConfiguration('inclusive/weakly_isolated/diMuonMass', ''))[0]
			ss_events = SingleIntegral(self.data.getHistogramWithConfiguration('inclusive/weakly_isolated/diMuonMass', 'ss'))[0]
			bkg_os_events = SingleIntegral(self.bkg.getHistogramWithConfiguration('inclusive/weakly_isolated/diMuonMass', ''))[0]
			bkg_ss_events = SingleIntegral(self.bkg.getHistogramWithConfiguration('inclusive/weakly_isolated/diMuonMass', 'ss'))[0]
			self.ss_os_factor = (os_events - bkg_os_events) / (ss_events - bkg_ss_events)
			print 'OS/SS factor = %g' % self.ss_os_factor

		histo = self.data.getHistogramWithConfiguration(histogram, 'ss')
		bkg_histo = self.bkg.getHistogramWithConfiguration(histogram, 'ss')
		histo.Add(bkg_histo, -1.)
		histo.Scale(self.ss_os_factor)
		histo.SetFillColor(self.color)
		histo.SetLineColor(ROOT.kBlack)
		histo.SetLineWidth(2)
		return histo.Clone()

	def getSystUncertainties(self):
		return {
			'QCD': 0.20
		}

class DatasetData(DatasetFile):
	def __init__(self, label, filename):
		DatasetFile.__init__(self, label, filename)

	def getHistogramWithConfiguration(self, histogram, configuration):
		histo = DatasetFile.getHistogramWithConfiguration(self, histogram, configuration)
		histo.SetMarkerSize(1.0)
		histo.SetLineWidth(2)
		histo.SetMarkerStyle(20)
		histo.SetMarkerColor(ROOT.kBlack)
		histo.SetLineColor(ROOT.kBlack)
		return histo

mc_samples = [
	DatasetMC('WW', 'out-embedtest_Summer12_WW_TuneZ2star_8TeV_pythia6_tauola-v1', ROOT.kPink-7, 54.838),
	DatasetMC('ZZ', 'out-embedtest_Summer12_ZZ_TuneZ2star_8TeV_pythia6_tauola-v1', ROOT.kRed-7, 8.06),
	DatasetMC('WZ', 'out-embedtest_Summer12_WZ_TuneZ2star_8TeV_pythia6_tauola-v1', ROOT.kRed-9, 33.21),
	DatasetMC('t#bar{t} + Jets', 'out-embedtest_Summer12_TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola-v2', ROOT.kYellow-7, 234.0, 1.07),
	DatasetMCCombined('Z/#gamma* #rightarrow #tau#tau', ROOT.kCyan-7, [
		DatasetMC('Z/#gamma* #rightarrow #tau#tau (M_{ll} > 50 GeV)', 'out-embedtest_Summer12_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-v1-filter15', ROOT.kCyan-7, 3503.71),
		DatasetMC('Z/#gamma* #rightarrow #tau#tau (10 GeV < M_{ll} < 50 GeV)', 'out-embedtest_Summer12_DYJetsToLL_M-10To50filter_TuneZ2Star_8TeV-madgraph-v1-filter15', ROOT.kCyan-7, 905.6)]),
	DatasetMCCombined('Z/#gamma* #rightarrow #mu#mu', ROOT.kAzure-7, [
		DatasetMC('Z/#gamma* #rightarrow #mu#mu (M_{ll} > 50 GeV)', 'out-embedtest_Summer12_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-v1-filter13', ROOT.kAzure-7, 3503.71),
		DatasetMC('Z/#gamma* #rightarrow #mu#mu (10 GeV < M_{ll} < 50 GeV)', 'out-embedtest_Summer12_DYJetsToLL_M-10To50filter_TuneZ2Star_8TeV-madgraph-v1-filter13', ROOT.kAzure-7, 905.6)])
]

# W+Jets does not contribute to the significantly, but it is part of the background in the
# anti-isolated region used for the QCD estimate, so we need it as well:
wjets_sample = [
	DatasetMC('WJets', 'out-embedtest_Summer12_WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball-v2', ROOT.kBlack, 36257.2)
]

#data_sample = DatasetData('Run 2012A', ['out-embedtest_Data_DoubleMu_Run2012A_22Jan2013-v1'])
data_sample = DatasetData('19.7 fb^{-1}', ['out-embedtest_Data_DoubleMu_Run2012A_22Jan2013-v1', 'out-embedtest_Data_DoubleMu_Run2012B_22Jan2013-v1', 'out-embedtest_Data_DoubleMu_Run2012C_22Jan2013-v1', 'out-embedtest_Data_DoubleMu_Run2012D_22Jan2013-v1'])
qcd_sample = DatasetQCD('QCD Multijet', ROOT.kGray, data_sample, DatasetMCCombined('Non-QCD Background', ROOT.kCyan-7, mc_samples + wjets_sample))

background_samples = mc_samples[:]
background_samples.insert(3, qcd_sample)

cnv = ROOT.TCanvas() #'bla', 'blub', 0, 0, 768 + 4, 768 + 28)
cnv.SetFillColor(ROOT.kWhite)
cnv.SetBorderSize(0)
cnv.SetLineWidth(2)

HISTOGRAMS = [
#	('nPV', (0, 50, ''), (1, 1e7, True), 1),
#	('nPV_u', (0, 50, ''), (1, 1e7, True), 1),
#	('leadMuonPt', (17, 100, 'GeV'), (1, 1e7, True), 20),
#	('leadMuonEta', (-2.5, 2.5, ''), (1, 1e7, True), 1),
#	('leadMuonPhi', (-math.pi, math.pi, ''), (1, 1e7, True), 1),
	('leadMuonIso', (0.0, 1.0, ''), (1e4, 1e9, True), 2),
#	('subMuonPt', (8, 100, 'GeV'), (1, 1e7, True), 10),
#	('subMuonEta', (-2.5, 2.5, ''), (1, 1e7, True), 1),
#	('subMuonPhi', (-math.pi, math.pi, ''), (1, 1e7, True), 1),
	('subMuonIso', (0.0, 1.0, ''), (1e4, 1e9, True), 2),
	('diMuonPt', (0, 100, 'GeV'), (1, 1e7, True), 25),
#	('diMuonEta', (-9.0, 9.0, ''), (1, 1e7, True), 1),
#	('diMuonDeltaPhi', (0, math.pi, ''), (1, 1e7, True), 1),
#	('diMuonDeltaR', (0, 6.0, ''), (1, 1e7, True), 1),
	('diMuonMass', (50., 200., 'GeV'), (1e1, 2e6, True), 50),
#	('nTracks5', (0, 50, ''), (1, 1e7, True), 1),
#	('nTracks10', (0, 40, ''), (1, 1e7, True), 1),
#	('nTracks20', (0, 30, ''), (1, 1e7, True), 1),
#	('nTracks30', (0, 20, ''), (1, 1e7, True), 1),
#	('nTracks40', (0, 10, ''), (1, 1e7, True), 1),
#	('nJets30', (0, 10, ''), (1, 1e7, True), 1),
#	('bJetPt', (20, 100, 'GeV'), (1e-1, 1e5, True), 20),
#	('nGlobalMuons', (0, 10, ''), (1, 1e7, True), 1),
#	('nStandaloneMuons', (0, 10, ''), (1, 1e7, True), 1),
#	('nPFMuons', (0, 10, ''), (1, 1e7, True), 1),
#	('pfMet', (0, 100, 'GeV'), (1, 1e7, True), 20),
#	('pfMetType1', (0, 100, 'GeV'), (1, 1e7, True), 20),
]

CATEGORIES = ['inclusive', 'btag'] #, 'vbf']
CATEGORY_SCALE = { 'inclusive': 1.0, 'btag': 0.1, 'vbf': 0.01 }
CATEGORY_REBIN = { 'inclusive': 1, 'btag': 2, 'vbf': 4 }

makeRatio = False
for category in CATEGORIES:
	for name, (xmin, xmax, unit), (ymin, ymax, log), rebin in HISTOGRAMS:
		mumu_histos = {}
		for selection in ['ided', 'isolated']:
			# Only take isolation from ided
			if selection == 'isolated' and (name == 'subMuonIso' or name == 'leadMuonIso'): continue
			if selection == 'ided' and name != 'subMuonIso' and name != 'leadMuonIso' and name != 'diMuonMass': continue

			# Create the histos
			histos = []
			for sample in background_samples + [data_sample]:
				histo = sample.getHistogram('%s/%s/%s' % (category, selection, name))
				histo.Rebin(rebin * CATEGORY_REBIN[category])
				histo.GetXaxis().SetRangeUser(xmin+0.001, xmax-0.001)
				histo.GetYaxis().SetRangeUser(ymin * CATEGORY_SCALE[category], ymax * CATEGORY_SCALE[category])
				histo.GetYaxis().SetTitle("Events/%g %s" % (histo.GetBinWidth(1), unit))
				histo.GetYaxis().SetTitleOffset(1.4)
				histo.SetStats(0)

				# Divide by bin width and adapt title of Y axis
				for x in range(histo.GetNbinsX()):
					histo.SetBinContent(x+1, histo.GetBinContent(x+1) / histo.GetBinWidth(x+1))
					histo.SetBinError(x+1, histo.GetBinError(x+1) / histo.GetBinWidth(x+1))

				if unit != '':
					histo.GetYaxis().SetTitle('dN/d%s [1/%s]' % (histo.GetXaxis().GetTitle(), unit))
					histo.GetXaxis().SetTitle(histo.GetXaxis().GetTitle() + ' [%s]' % unit)
				else:
					histo.GetYaxis().SetTitle('dN/d%s' % (histo.GetXaxis().GetTitle()))

				title = histo.GetXaxis().GetTitle()
				if makeRatio:
					histo.GetXaxis().SetLabelSize(0)
					histo.GetXaxis().SetTitle("")

				if sample.getLabel() == 'Z/#gamma* #rightarrow #mu#mu':
					mumu_histos[selection] = histo.Clone()

				if sample != data_sample:
					if len(histos) > 0:
						histo.Add(histos[-1][0])
						histo.GetYaxis().SetRangeUser(ymin * CATEGORY_SCALE[category], ymax * CATEGORY_SCALE[category]) # .Add() messes with yrange, so fix it up

					histos.append((histo, sample.getLabel()))
				else:
					data_histo = (histo, sample.getLabel())

			if makeRatio:
				# Prepare pad
				pad = ROOT.TPad("bla","blub",0,0.25,1.0,1.0);
			else:
				pad = cnv
			pad.SetFillColor(0)
			pad.SetFillStyle(0)
			if log: pad.SetLogy(1)
			else: pad.SetLogy(0)
			pad.SetGrid()
			pad.Draw()
			pad.cd()

			# Draw main plot
			prevHisto = None
			for histo, label in reversed(histos):
				if prevHisto is not None:
					histo.DrawCopy('hist same')
				else:
					histo.DrawCopy('hist')
				prevHisto = histo
			data_histo[0].Draw('esame')

			leg = ROOT.TLegend(0.60, 0.58, 0.92, 0.93)
			leg.SetFillColor(ROOT.kWhite)
			leg.SetLineWidth(0)
			leg.SetBorderSize(0)
			leg.AddEntry(data_histo[0], data_histo[1], 'LP')
			for histo, label in reversed(histos):
				leg.AddEntry(histo, label, 'F')
			leg.Draw()
			pad.RedrawAxis()

			if makeRatio:
				# Draw ratio plot
				cnv.cd()
				pad2 = ROOT.TPad("bla","blub",0,0.05,1.0,0.325)
				pad2.SetTopMargin(0.0)
				pad2.SetFillColor(0)
				pad2.SetFillStyle(0)
				pad2.SetGrid()
				pad2.Draw()
				pad2.cd()

				# TODO: Use a TAsymmErrorGraph and PC errors
				ratio_histo = data_histo[0].Clone('%s_Ratio' % data_histo[0].GetName())
				ratio_histo.Divide(histos[-1][0])
				#ratio_histo.Add(-1)
				ratio_histo.SetTitle('')
				ratio_histo.SetMarkerStyle(20)
				ratio_histo.SetMarkerSize(1)
				ratio_histo.SetMarkerColor(ROOT.kBlack)
				ratio_histo.SetLineColor(ROOT.kBlack)
				ratio_histo.SetStats(0)
				ratio_histo.GetXaxis().SetTitle(title)
				ratio_histo.GetXaxis().SetLabelSize(0.10)
				ratio_histo.GetXaxis().SetTitleSize(0.10)
				ratio_histo.GetXaxis().SetRangeUser(xmin, xmax-0.001)
				ratio_histo.GetYaxis().SetLabelSize(0.10)
				ratio_histo.GetYaxis().SetTitleSize(0.12)
				ratio_histo.GetYaxis().SetTitleOffset(0.4)
				ratio_histo.GetYaxis().SetTitle('Data / MC')
				ratio_histo.GetYaxis().SetRangeUser(0.5,1.99)

				line = ROOT.TLine(xmin, 1, xmax, 1)
				line.SetLineColor(ROOT.kRed)
				line.SetLineWidth(2)

				# Draw syst. uncertainties
				nbins = data_histo[0].GetNbinsX()
				syst_uncertainties = [[None for x in range(len(background_samples))] for y in range(nbins)]
				for index, sample in enumerate(background_samples):
					systematics = sample.getSystUncertainties()
					for i in range(0, nbins):
						entries = histos[index][0].GetBinContent(i+1)
						if index > 0: entries -= histos[index-1][0].GetBinContent(i+1)
						syst_uncertainties[i][index] = dict([ (source, value*entries) for source, value in systematics.iteritems()])

				sources = set(item for sublist in [x.keys() for x in syst_uncertainties[0]] for item in sublist)
				systematics_added = [dict([(source, reduce(lambda y,z: y+z, [x.get(source, 0.0) for x in unc])) for source in sources]) for unc in syst_uncertainties]
				systematics_combined = [reduce(lambda y,z: (y**2+z**2)**0.5, unc.values()) for unc in systematics_added]

				def save_div(x, y):
					if y == 0.0: return 0
					return x/y

				x_ = [data_histo[0].GetBinCenter(i+1) for i in range(0, nbins)]
				y_ = numpy.ones(nbins)
				xe_ = [0.5*(data_histo[0].GetXaxis().GetBinUpEdge(i+1) - data_histo[0].GetXaxis().GetBinLowEdge(i+1)) for i in range(0, nbins)]
				ye_ = [save_div(x, histos[-1][0].GetBinContent(i+1)) for i,x in enumerate(systematics_combined)]
				graph = ROOT.TGraphErrors(len(x_), numpy.array(x_), y_, numpy.array(xe_), numpy.array(ye_))
				graph.SetFillColor(ROOT.kGray)
				graph.SetTitle("")

				ratio_histo.Draw("p") # I would like this to draw nothing at all, only the axis with labels, but I have no idea how to do that. "AXIG" does not draw the axis labels
				graph.Draw("2")
				pad2.RedrawAxis("g")
				line.Draw()
				ratio_histo.Draw("esame")
				pad2.RedrawAxis()

				graph.GetXaxis().SetLabelSize(0.10)
				graph.GetXaxis().SetTitleSize(0.10)
				graph.GetXaxis().SetRangeUser(xmin, xmax-0.001)
				graph.GetYaxis().SetLabelSize(0.10)
				graph.GetYaxis().SetTitleSize(0.12)
				graph.GetYaxis().SetTitleOffset(0.3)
				graph.GetYaxis().SetTitle('Data / MC')
				graph.GetYaxis().SetRangeUser(0.5,1.99)

			cnv.cd()
			cnv.SaveAs('plots/%s-%s-%s.png' % (category, selection, name))
			cnv.SaveAs('plots-pdf/%s-%s-%s.pdf' % (category, selection, name))
			cnv.Clear()

"""
			# TODO: Needs to be fixed!
			all = all_err2 = 0
			for histo, label in reversed(histos):
				int, int_err2 = SingleIntegral(histo)
				all += int
				all_err2 += int_err2

			for histo, label in reversed(histos):
				int, int_err2 = SingleIntegral(histo)
				#d, d_err2 = stattest.getPC(int, all, int_err2, all_err2) TODO
				d = int / all
				d_err2 = 0
				print label + ':', d, '+/-', d_err2**0.5
"""

"""
		for num, denom in [('masscut', 'ided'), ('isolated', 'ided'), ('isolated_masscut', 'masscut')]:
			pad = ROOT.TPad("bla","blub",0,0.25,1.0,1.0);
			pad.SetFillColor(0)
			pad.SetFillStyle(0)
			if log: pad.SetLogy(1)
			else: pad.SetLogy(0)
			pad.Draw()
			pad.cd()

			# Let's go and make a isolated vs. non-isolated comparison plot of this, whatever it is :)
			mumu_histos[num].SetMarkerStyle(21)
			mumu_histos[num].SetLineColor(ROOT.kRed)
			mumu_histos[num].SetMarkerColor(ROOT.kRed)
			mumu_histos[num].Scale(1/SingleIntegral(mumu_histos[num])[0])
			mumu_histos[num].GetYaxis().SetRangeUser(1e-5, 1e0)
			mumu_histos[num].Draw('e')
			mumu_histos[denom].SetMarkerStyle(22)
			mumu_histos[denom].SetLineColor(ROOT.kBlue)
			mumu_histos[denom].SetMarkerColor(ROOT.kBlue)
			mumu_histos[denom].Scale(1/SingleIntegral(mumu_histos[denom])[0])
			mumu_histos[denom].Draw("esame")

			leg = ROOT.TLegend(0.60, 0.60, 0.90, 0.85)
			leg.SetFillColor(ROOT.kWhite)
			leg.SetLineWidth(0)
			leg.SetBorderSize(0)
			leg.AddEntry(mumu_histos[num], '%s Z/#gamma* #rightarrow #mu#mu' % num, 'L')
			leg.AddEntry(mumu_histos[denom], '%s Z/#gamma* #rightarrow #mu#mu' % denom, 'L')
			leg.Draw()

			cnv.cd()
			pad2 = ROOT.TPad("bla","blub",0,0.02,1.0,0.25)
			pad2.SetFillColor(0)
			pad2.SetFillStyle(0)
			pad2.SetGridx()
			pad2.SetGridy()
			pad2.Draw()
			pad2.cd()

			# TODO: Use a TAsymmErrorGraph and PC errors
			ratio_histo = mumu_histos[num].Clone()
			ratio_histo.Divide(mumu_histos[denom].Clone())
			ratio_histo.SetName('Ratio')
			ratio_histo.SetMarkerStyle(20)
			ratio_histo.SetMarkerSize(1)
			ratio_histo.SetMarkerColor(ROOT.kBlack)
			ratio_histo.SetLineColor(ROOT.kBlack)
			ratio_histo.GetXaxis().SetLabelSize(0.12)
			ratio_histo.GetXaxis().SetTitle('')
			ratio_histo.GetYaxis().SetLabelSize(0.10)
			ratio_histo.GetYaxis().SetTitle('%s/%s' % (num, denom))
			ratio_histo.GetYaxis().SetTitleSize(0.10)
			ratio_histo.GetYaxis().SetTitleOffset(0.5)
			ratio_histo.GetYaxis().SetRangeUser(0.7,1.3)
			ratio_histo.SetStats(0)
			ratio_histo.Draw("e")

			line = ROOT.TLine(ratio_histo.GetXaxis().GetBinLowEdge(1),1,ratio_histo.GetXaxis().GetBinUpEdge(ratio_histo.GetXaxis().GetLast()), 1)
			line.SetLineStyle(2)
			line.SetLineColor(ROOT.kRed)
			line.Draw()

			cnv.cd()
			cnv.SaveAs('plots/%s-%s_vs_%s-%s.png' % (category, num, denom, name))
			cnv.Clear()
"""
