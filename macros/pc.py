import ROOT
import numpy

def getPCConfInt(k, N, alpha = 0.32):
	if N == 0: return (0., 0.)

	if k == N: return (ROOT.Math.beta_quantile(.5 * alpha, k, N + 1 - k), 1.0)
	if k == 0: return (0., ROOT.Math.beta_quantile(1. - .5 * alpha, k + 1, N - k))

	result = float(k) / float(N)
	return (ROOT.Math.beta_quantile(.5 * alpha, k, N + 1 - k), ROOT.Math.beta_quantile(1. - .5 * alpha, k + 1, N - k))

def makeEfficiencyGraph(numerator, denominator):
	assert numerator.GetNbinsX() == denominator.GetNbinsX()

	x = []
	y = []
	exl = []
	exh = []
	eyl = []
	eyh = []
	for i in range(1, numerator.GetNbinsX()+1):
		numEff = numerator.GetBinContent(i)**2 / numerator.GetBinError(i)**2
		denomEff = denominator.GetBinContent(i)**2 / denominator.GetBinError(i)**2
		xval = numerator.GetXaxis().GetBinCenter(i)
		yval = yvalEff = 0.

		if denominator.GetBinContent(i) > 0.:
			yval = numerator.GetBinContent(i) / denominator.GetBinContent(i)
			yvalEff = numEff / denomEff

		ylowEff, yhighEff = getPCConfInt(numEff, denomEff)

		if yvalEff > 0.:
			ylowRel = ylowEff / yvalEff
		else:
			ylowRel = 0.

		if yvalEff < 1.:
			yhighRel = (1. - yhighEff) / (1. - yvalEff)
		else:
			yhighRel = 0.

		ylow = yval * ylowRel
		yhigh = 1. - yhighRel * (1. - yval)

		x.append(xval)
		y.append(yval)
		exl.append(xval - numerator.GetXaxis().GetBinLowEdge(i))
		exh.append(numerator.GetXaxis().GetBinUpEdge(i) - xval)
		eyl.append(yval - ylow)
		eyh.append(yhigh - yval)

	graph = ROOT.TGraphAsymmErrors(len(x), numpy.array(x), numpy.array(y), numpy.array(exl), numpy.array(exh), numpy.array(eyl), numpy.array(eyh))
	graph.GetXaxis().SetTitle(numerator.GetXaxis().GetTitle())
	graph.GetYaxis().SetTitle(numerator.GetYaxis().GetTitle())
	return graph

def makeEfficiencyGraphRatio(graph1, graph2):
	assert(graph1.GetN() == graph2.GetN())

	x = []
	y = []
	exl = []
	exh = []
	eyl = []
	eyh = []

	for i in range(graph1.GetN()):
		assert(graph1.GetX()[i] == graph2.GetX()[i])
		assert(graph1.GetEXlow()[i] == graph2.GetEXlow()[i])
		assert(graph1.GetEXhigh()[i] == graph2.GetEXhigh()[i])

		y1 = graph1.GetY()[i]
		y2 = graph2.GetY()[i]
		y1err = max(graph1.GetEYlow()[i], graph1.GetEYhigh()[i])
		y2err = max(graph2.GetEYlow()[i], graph2.GetEYhigh()[i])

		if y2 > 0.:
			yerr = ((y1err/y2)**2 + ((y1*y2err)/y2**2)**2)**0.5

			x.append(graph1.GetX()[i])
			exl.append(graph1.GetEXlow()[i])
			exh.append(graph1.GetEXhigh()[i])
			y.append(y1 / y2)
			eyl.append(yerr)
			eyh.append(yerr)

	graph = ROOT.TGraphAsymmErrors(len(x), numpy.array(x), numpy.array(y), numpy.array(exl), numpy.array(exh), numpy.array(eyl), numpy.array(eyh))
	graph.GetXaxis().SetTitle(graph1.GetXaxis().GetTitle())
	graph.GetYaxis().SetTitle(graph1.GetYaxis().GetTitle())
	return graph
