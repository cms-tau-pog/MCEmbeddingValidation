#!/usr/bin/env python

import sys
sys.argv.append('-b')

import ROOT
ROOT.gStyle.SetLineWidth(2)
ROOT.gStyle.SetCanvasBorderSize(2)
ROOT.gStyle.SetPadBorderSize(2)
ROOT.gStyle.SetHistLineWidth(2)
ROOT.gStyle.SetFrameLineWidth(2)
import numpy
import math
import pc

ROOT.gROOT.ProcessLineSync('.x tdrStyle.C')
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)

import errno
import os

def prepare_histo(histo):
	histo.GetXaxis().SetLabelOffset(0.005)
	histo.GetXaxis().SetTitleOffset(1.0)
	histo.GetXaxis().SetTitleFont(62)
	histo.GetXaxis().SetTitleSize(0.04)
	histo.GetXaxis().SetLabelFont(62)
	histo.GetXaxis().SetLabelSize(0.04)
	histo.GetYaxis().SetTitleFont(62)
	histo.GetYaxis().SetTitleSize(0.05)
	histo.GetYaxis().SetLabelFont(62)
	histo.GetYaxis().SetLabelSize(0.05)

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

	def getLabel(self):
		return self.label

class DatasetFile(Dataset):
	def __init__(self, label, filename):
		Dataset.__init__(self, label)
		self.filename = filename
		self.files = None

	def getHistogram(self, histogram):
		if self.files is None:
			if type(self.filename) == list:
				self.files = [ROOT.TFile('%s.root' % (filename), 'READ') for filename in self.filename]
			else:
				self.files = [ROOT.TFile('%s.root' % (self.filename), 'READ')]

		histo = None
		for file in self.files:
			file_histo = file.Get(histogram)
			if not file_histo:
				raise Exception('Histogram \"%s\" does not exist in file \"%s\"' % (histogram, file.GetName()))

			cur_histo = file_histo.Clone()
			if histo == None:
				histo = cur_histo
			else:
				histo.Add(cur_histo)

		return histo

class DatasetEmbTest(DatasetFile):
	def __init__(self, label, filename):
		DatasetFile.__init__(self, label, filename)

	def getHistogram(self, histogram):
		histo = DatasetFile.getHistogram(self, histogram)
		histo.SetLineWidth(2)
		prepare_histo(histo)
		return histo

def extract(plots, name):
	return {name: plots[name]}

def mkdir_p(name):
	try:
		os.mkdir(name)
	except OSError, e:
		if e.errno != errno.EEXIST:
			raise e

def create_plots(plots):
	cnv = ROOT.TCanvas()
	cnv.SetFillColor(ROOT.kWhite)
	cnv.SetBorderSize(0)
	cnv.SetLineWidth(2)

	for plotName in plots:
		print '%s...' % plotName
		plot = plots[plotName]
		xmin, xmax, unit = plot['xaxis']
		ymin, ymax, ylog = plot['yaxis']
		binning = plot['binning']
		samples = plot['samples']

		if plot['normalization'] is not None:
			normalizationCategory, normalizationHistogramName = plot['normalization']

		def get_norm_factor(dataset, category, histogramName):
			def get_normalization_category():
				if normalizationCategory == '*': return category
				return normalizationCategory

			def get_normalization_histogram_name():
				if normalizationHistogramName == '*': return histogramName
				return normalizationHistogramName

			# Get histogram
			normCategory = get_normalization_category()
			normHistogramName = get_normalization_histogram_name()
			histo = dataset.getHistogram('%s/%s' % (normCategory, normHistogramName))

			bin1 = histo.FindBin(xmin+1e-3)
			bin2 = histo.FindBin(xmax-1e-3)
			return histo.Integral(bin1, bin2)

#			if histogramName == normHistogramName:
#				if type(binning) == list:
#					histo = histo.Rebin(len(binning)-1, "%s_rebinned_norm" % histo.GetName(), numpy.array(binning))
#				else:
#					histo.Rebin(binning)

			# Divide by bin width
#			before = SingleIntegral(histo)[0]
#			for x in range(histo.GetNbinsX()):
#				histo.SetBinContent(x+1, histo.GetBinContent(x+1) / histo.GetBinWidth(x+1))
#				histo.SetBinError(x+1, histo.GetBinError(x+1) / histo.GetBinWidth(x+1))
#			return SingleIntegral(histo)[0]

		normCache = {}
		def get_norm_ratio_for_sample(sample):
			sampleNorm = sample[6]
			if sampleNorm is None:
				(normNum, normDenom) = (0, samples.index(sample))
			else:
				(normNum, normDenom) = sampleNorm

			if normNum not in normCache:
				if normNum != -1:
					normCache[normNum] = get_norm_factor(samples[normNum][0], samples[normNum][2], samples[normNum][3])
				else:
					normCache[normNum] = 1.

			if normDenom not in normCache:
				if normDenom != -1:
					normCache[normDenom] = get_norm_factor(samples[normDenom][0], samples[normDenom][2], samples[normDenom][3])
				else:
					normCache[normDenom] = 1.

			return normCache[normNum] / normCache[normDenom]

		def make_histo_from_sample(sample):
			(dataset, label, category, histogramName, color, style, norm) = sample
			def get_single_histogram(histoName):
				if type(histoName) == list:
					histo = dataset.getHistogram('%s/%s' % (category, histoName[0]))
					for name in histoName[1:]:
						histo.Add(dataset.getHistogram('%s/%s' % (category, name)))
				else:
					histo = dataset.getHistogram('%s/%s' % (category, histoName))

				# rebin it
				if type(binning) == list:
					histo = histo.Rebin(len(binning)-1, "%s_rebinned" % histo.GetName(), numpy.array(binning))
				else:
					histo.Rebin(binning)

				# Divide by bin width
				for x in range(histo.GetNbinsX()):
					histo.SetBinContent(x+1, histo.GetBinContent(x+1) / histo.GetBinWidth(x+1))
					histo.SetBinError(x+1, histo.GetBinError(x+1) / histo.GetBinWidth(x+1))

				# Normalize
				if plot['normalization'] is not None:
					normRatio = get_norm_ratio_for_sample(sample)
					histo.Scale(normRatio)

				return histo

			isEfficiency = (type(histogramName) == tuple)
			if isEfficiency:
				denomHistoName, numHistoName = histogramName
				numHisto = get_single_histogram(numHistoName)
				denomHisto = get_single_histogram(denomHistoName)

				histo = pc.makeEfficiencyGraph(numHisto, denomHisto)
				histo.SetLineWidth(denomHisto.GetLineWidth())
				prepare_histo(histo)
			else:
				histo = get_single_histogram(histogramName)

			# Setup histogram graphics
			histo.GetXaxis().SetRangeUser(xmin+1e-3, xmax-1e-3)
			histo.GetYaxis().SetRangeUser(ymin-ymin/1e-5, ymax)

			xlabel = histo.GetXaxis().GetTitle()
			ylabel = "dN/d%s" % xlabel # note this uses the original version of the x label!!
			if 'xlabel' in plot:
				xlabel = plot['xlabel']
			if 'ylabel' in plot:
				ylabel = plot['ylabel']

			# Set up axis title
			if isEfficiency:
				histo.GetYaxis().SetTitle(ylabel)
			else:
				if unit != '':
					histo.GetYaxis().SetTitle("%s [1/%s]" % (ylabel, unit))
				else:
					histo.GetYaxis().SetTitle("%s" % ylabel)
			if unit != '':
				histo.GetXaxis().SetTitle("%s [%s]" % (xlabel, unit))
			else:
				histo.GetXaxis().SetTitle("%s" % xlabel)

			histo.GetYaxis().SetTitleOffset(1.4)
			histo.GetXaxis().SetLabelSize(0)
			if type(histo) != ROOT.TGraphAsymmErrors: histo.SetStats(0)

			histo.SetTitle(label)
			histo.SetLineColor(color)
			histo.SetMarkerColor(color)
			histo.SetMarkerStyle(style)

			histo.GetYaxis().SetRangeUser(ymin-ymin/1e5, ymax) # reset, since it is overwritten
			return histo

		def make_ratio_histo(histo, histo0):
			if type(histo) == ROOT.TGraphAsymmErrors:
				ratio_histo = pc.makeEfficiencyGraphRatio(histo, histo0)
				prepare_histo(ratio_histo)
				ratio_histo.SetLineWidth(histo.GetLineWidth())
				ratio_histo.SetLineColor(histo.GetLineColor())
				ratio_histo.SetMarkerColor(histo.GetMarkerColor())
			else:
				ratio_histo = histo.Clone('%s_Ratio' % histo.GetName())
				ratio_histo.Divide(histo0)

			ratio_histo.SetTitle('')
			ratio_histo.SetMarkerSize(1)
			if type(ratio_histo) != ROOT.TGraphAsymmErrors: ratio_histo.SetStats(0)
			ratio_histo.GetXaxis().SetLabelSize(0.11)
			ratio_histo.GetXaxis().SetLabelOffset(0.02)
			ratio_histo.GetXaxis().SetTitleOffset(1.25)
			ratio_histo.GetXaxis().SetTitleSize(0.11)
			#ratio_histo.GetXaxis().SetTitle('')
			ratio_histo.GetXaxis().SetRangeUser(xmin+0.001, xmax-0.001)
			ratio_histo.GetYaxis().SetLabelSize(0.08)
			ratio_histo.GetYaxis().SetTitleSize(0.11)
			ratio_histo.GetYaxis().SetTitleOffset(0.65)

			ratiolabel = 'Embedded / MC'
			if 'ratiolabel' in plot:
				ratiolabel = plot['ratiolabel']
			ratio_histo.GetYaxis().SetTitle(ratiolabel)

			ratio_range = (0.7, 1.3)
			if 'ratio-range' in plot:
				ratio_range = plot['ratio-range']

#			if plotName == 'rad-rec-ditau-svfit-mass-talk':
#				f = ROOT.TF1('f', '[0]+[1]*x', 60., 120.)
#				f.SetLineWidth(2)
#				if ratio_histo.GetLineColor() == ROOT.kBlue:
#					f.SetLineColor(ROOT.kAzure+7)
#				else:
#					f.SetLineColor(ROOT.kOrange+7)
#				ratio_histo.Fit(f)

			ratio_histo.GetYaxis().SetRangeUser(ratio_range[0]+1e-3,ratio_range[1]-1e-3)
			return ratio_histo

		def make_plot(samples, histos):
			# Prepare pad
			pad = ROOT.TPad("bla","blub",0,0.25,1.0,1.0);
			pad.SetFillColor(0)
			pad.SetFillStyle(0)
			if ylog: pad.SetLogy(1)
			else: pad.SetLogy(0)
			pad.SetGrid()
			pad.Draw()
			pad.cd()

			ratio_histos = [make_ratio_histo(histo, histos[0]) for histo in histos[1:]]
			histos[0].GetXaxis().SetTitle('')

#			if plotName == 'bug' and False:
#				f = ROOT.TF1('f', '[0]+[1]*x', 0.1, 1.0)
#				f.SetLineWidth(2)
#				if histos[0].GetLineColor() == ROOT.kBlue:
#					f.SetLineColor(ROOT.kAzure+7)
#				else:
#					f.SetLineColor(ROOT.kOrange+7)
#				histos[0].Fit(f)
#				area = 0.5 * f.GetParameter(1) + f.GetParameter(0)
#				print 'Slope: %g' % ( (f.Eval(1.0)/area - f.Eval(0.1)/area) / (1.0 - 0.1))

			# Draw main plot
			if type(histos[0]) == ROOT.TGraphAsymmErrors:
				histos[0].Draw("AP")
				for hist in histos[1:]:
					hist.Draw('psame')
			else:
				histos[0].DrawCopy()
				for hist in histos[1:]:
					hist.DrawCopy('same')

			legend_width = 0.35
			legend_height = 2./30. * len(histos)
			if 'legend-width' in plot:
				legend_width = plot['legend-width']
			if 'legend-height' in plot:
				legend_height = plot['legend-height']
			legend_pos = (0.95 - legend_width, 0.92 - legend_height)
			if 'legend-pos' in plot:
				legend_pos = plot['legend-pos']

			leg = ROOT.TLegend(legend_pos[0], legend_pos[1], legend_pos[0] + legend_width, legend_pos[1] + legend_height)
			leg.SetFillColor(ROOT.kWhite)
			leg.SetLineWidth(0)
			leg.SetBorderSize(0)
			for histo in histos:
				leg.AddEntry(histo, histo.GetTitle(), 'LP')
			leg.Draw()
			pad.RedrawAxis()

			# Draw ratio plot
			cnv.cd()
			#cnv.SetBottomMargin(-0.50)
			pad2 = ROOT.TPad("bla","blub",0,0.0,1.0,0.335)
			pad2.SetBottomMargin(0.35) # makes x axis label visible in PDF
			pad2.SetTopMargin(0.0)
			pad2.SetFillColor(0)
			pad2.SetFillStyle(0)
			pad2.SetGrid()
			pad2.Draw()
			pad2.cd()

			line = ROOT.TLine(xmin, 1, xmax, 1)
			line.SetLineColor(ROOT.kBlack)
			line.SetLineStyle(2)
			line.SetLineWidth(2)

			if type(ratio_histos[0]) == ROOT.TGraphAsymmErrors:
				ratio_histos[0].Draw("AP")
				for histo in ratio_histos[1:]:
					histo.Draw("psame")
			else:
				ratio_histos[0].Draw()
				for histo in ratio_histos[1:]:
					histo.Draw("same")
			line.Draw()

			mkdir_p('plots')
			mkdir_p('plots/png')
			mkdir_p('plots/pdf')

			cnv.cd()
			cnv.SaveAs('plots/png/%s.png' % (plotName))
			cnv.SaveAs('plots/pdf/%s.pdf' % (plotName))
			cnv.Clear()

		histos = [make_histo_from_sample(sample) for sample in samples]
		make_plot(samples, histos)
