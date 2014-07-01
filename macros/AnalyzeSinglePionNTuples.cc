#include <TTree.h>
#include <TFile.h>
#include <TMath.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <Math/VectorUtil.h>
#include <sstream>
#include <memory>

static TLorentzVector TLorentzVector_PtEtaPhiE(double pt, double eta, double phi, double e)
{
	TLorentzVector v;
	v.SetPtEtaPhiE(pt, eta, phi, e);
	return v;
}

const bool lorentzEqual(const TLorentzVector& first, const TLorentzVector& second)
{
	return fabs(first.Px() - second.Px()) < 1e-3 && fabs(first.Py() - second.Py()) < 1e-3 && fabs(first.Pz() - second.Pz()) < 1e-3 && fabs(first.E() - second.E()) < 1e-3;
}

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > RotLorentzVector;

std::pair<double, double> xInRF(const RotLorentzVector& posPi, const RotLorentzVector& posTau,
                                const RotLorentzVector& negPi, const RotLorentzVector& negTau)
{
  RotLorentzVector zP4_lab = posTau + negTau;
  ROOT::Math::Boost boost_to_rf(zP4_lab.BoostToCM());

  RotLorentzVector posTau_rf = boost_to_rf(posTau);
  RotLorentzVector negTau_rf = boost_to_rf(negTau);
  RotLorentzVector posPi_rf = boost_to_rf(posPi);
  RotLorentzVector negPi_rf = boost_to_rf(negPi);

  return std::make_pair(posPi_rf.P() / posTau_rf.E(), negPi_rf.P() / negTau_rf.E());
}

std::pair<double, double> xInRF(const TLorentzVector& posPi, const TLorentzVector& posTau,
                                const TLorentzVector& negPi, const TLorentzVector& negTau)
{
  RotLorentzVector rPosPi(posPi.Px(), posPi.Py(), posPi.Pz(), posPi.E());
  RotLorentzVector rPosTau(posTau.Px(), posTau.Py(), posTau.Pz(), posTau.E());
  RotLorentzVector rNegPi(negPi.Px(), negPi.Py(), negPi.Pz(), negPi.E());
  RotLorentzVector rNegTau(negTau.Px(), negTau.Py(), negTau.Pz(), negTau.E());
  return xInRF(rPosPi, rPosTau, rNegPi, rNegTau);
}

static int findBin(TAxis* axis, double x)
{
	int bin = axis->FindFixBin(x);
	if(bin < 1) bin = 1;
	if(bin > axis->GetNbins()) bin = axis->GetNbins();
	return bin;
}

struct TreeVars {
	TreeVars() {
		run = lumi = event = ~0u;
		minVisPtFilterWeight = -1.0f;
		tauSpinWeight = -1.0f;
		zmumuEvtSelEffCorrWeight = -1.0f;
		nTrueInteractions = -1.0f;

		PosTau.GenPt = PosTau.GenEta = PosTau.GenPhi = PosTau.GenE = -1.0f;
		PosTau.GenFinalPt = PosTau.GenFinalEta = PosTau.GenFinalPhi = PosTau.GenFinalE = -1.0f;
		PosTau.GenVisPt = PosTau.GenVisEta = PosTau.GenVisPhi = PosTau.GenVisE = -1.0f;
		PosTau.OrigGenPt = PosTau.OrigGenEta = PosTau.OrigGenPhi = PosTau.OrigGenE = -1.0f;
		PosTau.OrigGenFinalPt = PosTau.OrigGenFinalEta = PosTau.OrigGenFinalPhi = PosTau.OrigGenFinalE = -1.0f;
		PosTau.GenVisPDG = PosTau.OrigGenChargedLepPDG = 0;

		NegTau.GenPt = NegTau.GenEta = NegTau.GenPhi = NegTau.GenE = -1.0f;
		NegTau.GenFinalPt = NegTau.GenFinalEta = NegTau.GenFinalPhi = NegTau.GenFinalE = -1.0f;
		NegTau.GenVisPt = NegTau.GenVisEta = NegTau.GenVisPhi = NegTau.GenVisE = -1.0f;
		NegTau.OrigGenPt = NegTau.OrigGenEta = NegTau.OrigGenPhi = NegTau.OrigGenE = -1.0f;
		NegTau.OrigGenFinalPt = NegTau.OrigGenFinalEta = NegTau.OrigGenFinalPhi = NegTau.OrigGenFinalE = -1.0f;
		NegTau.GenVisPDG = NegTau.OrigGenChargedLepPDG = 0;
	}

	unsigned int run;
	unsigned int lumi;
	unsigned int event;

	float minVisPtFilterWeight;
	float tauSpinWeight;
	float zmumuEvtSelEffCorrWeight;
	float nTrueInteractions;

	struct Tau {
		float GenPt;
		float GenEta;
		float GenPhi;
		float GenE;

		float GenFinalPt;
		float GenFinalEta;
		float GenFinalPhi;
		float GenFinalE;

		float GenVisPt;
		float GenVisEta;
		float GenVisPhi;
		float GenVisE;

		float OrigGenPt;
		float OrigGenEta;
		float OrigGenPhi;
		float OrigGenE;

		float OrigGenFinalPt;
		float OrigGenFinalEta;
		float OrigGenFinalPhi;
		float OrigGenFinalE;

		int OrigGenChargedLepPDG;
		int GenVisPDG;
	};

	Tau PosTau;
	Tau NegTau;
};

class Tau
{
public:
	Tau(const TreeVars::Tau& tau, int charge):
		hasOrig(tau.OrigGenFinalPt > 0.),
		hasVis(tau.GenVisPt > 0.),
		p4genInitial(TLorentzVector_PtEtaPhiE(tau.GenPt, tau.GenEta, tau.GenPhi, tau.GenE)),
		p4genFinal(TLorentzVector_PtEtaPhiE(tau.GenFinalPt, tau.GenFinalEta, tau.GenFinalPhi, tau.GenFinalE)),
		p4genVis(TLorentzVector_PtEtaPhiE(tau.GenVisPt, tau.GenVisEta, tau.GenVisPhi, tau.GenVisE)),
		p4origGenInitial(TLorentzVector_PtEtaPhiE(tau.OrigGenPt, tau.OrigGenEta, tau.OrigGenPhi, tau.OrigGenE)),
		p4origGenFinal(TLorentzVector_PtEtaPhiE(tau.OrigGenFinalPt, tau.OrigGenFinalEta, tau.OrigGenFinalPhi, tau.OrigGenFinalE)),
		genCharge(charge)
	{
	}

	float x() const
	{
		const double x = p4genVis.E() / p4genInitial.E();
		//assert(x >= 0. && x <= 1.);
		return x;
	}

	const bool hasOrig;
	const bool hasVis;

	const TLorentzVector p4genInitial;
	const TLorentzVector p4genFinal;
	const TLorentzVector p4genVis;
	const TLorentzVector p4origGenInitial;
	const TLorentzVector p4origGenFinal;
	const int genCharge;
};

class ZmumuEvtSelEffCorrWeight
{
public:
	ZmumuEvtSelEffCorrWeight(const char* filename)
	{
		file = new TFile(filename, "READ");
		if(file->IsZombie()) throw std::runtime_error( ("Could not load " + std::string(filename)).c_str());

		muMinusPt_vs_muPlusPt = dynamic_cast<TH2D*>(file->Get("EMB_ZmumuEvtSelEff_muMinusPt_vs_muPlusPt"));
		muMinusEta_vs_muPlusEta = dynamic_cast<TH2D*>(file->Get("EMB_ZmumuEvtSelEffCorr_muMinusEta_vs_muPlusEta"));
		if(!muMinusPt_vs_muPlusPt || !muMinusEta_vs_muPlusEta)
			throw std::runtime_error("Failed to load ZmumuEvtSelEff histograms");
	}

	double getWeight(const TLorentzVector& plus, const TLorentzVector& minus) const
	{
		const int binPtX = findBin(muMinusPt_vs_muPlusPt->GetXaxis(), plus.Pt());
		const int binPtY = findBin(muMinusPt_vs_muPlusPt->GetYaxis(), minus.Pt());
		double effPt = muMinusPt_vs_muPlusPt->GetBinContent(binPtX, binPtY);

		const int binEtaX = findBin(muMinusEta_vs_muPlusEta->GetXaxis(), plus.Eta());
		const int binEtaY = findBin(muMinusEta_vs_muPlusEta->GetYaxis(), minus.Eta());
		double effEta = muMinusEta_vs_muPlusEta->GetBinContent(binEtaX, binEtaY);

		if(effPt < 1e-1) { /*std::cout << "etaplus=" << plus.Eta() << ", etaminus=" << minus.Eta() << std::endl; */return 1.; }
		if(effEta < 1e-1) { /*std::cout << "etaplus=" << plus.Eta() << ", etaminus=" << minus.Eta() << std::endl; */return 1.; }
		return 1./(effPt * effEta);
	}

private:
	TFile* file;
	TH2D* muMinusPt_vs_muPlusPt;
	TH2D* muMinusEta_vs_muPlusEta;
};

struct Acceptance
{
	bool acceptEmb(const TLorentzVector& first, const TLorentzVector& second)
	{
		// Embedding cuts (very loose)
		if(first.Pt() < 8 || second.Pt() < 8) return false;
		else if(first.Pt() < 17 && second.Pt() < 17) return false;
		else if(fabs(first.Eta()) > 2.5) return false;
		else if(fabs(second.Eta()) > 2.5) return false;
		else if((first + second).M() < 50.0) return false;
		return true;
	}

	bool acceptRec(const TLorentzVector& first, const TLorentzVector& second)
	{
		// Tighter cuts with well-defined efficiency away from the trigger threshold
		if(first.Pt() < 10 || second.Pt() < 10) return false;
		else if(first.Pt() < 20 && second.Pt() < 20) return false;
		else if(fabs(first.Eta()) > 2.1) return false;
		else if(fabs(second.Eta()) > 2.1) return false;
		else if((first + second).M() < 50.0) return false;
		return true;
	}

	Acceptance(const Tau& first, const Tau& second):
		embAccepted(acceptEmb(first.p4genInitial, second.p4genInitial)),
		recAccepted(acceptRec(first.p4genInitial, second.p4genInitial))
	{
	}

	const bool embAccepted;
	const bool recAccepted;
};

struct RadiationHistograms {
	RadiationHistograms(TDirectory* dir, const char* name) {
		TDirectory* newdir = dir->mkdir(name);
		newdir->cd();

		const double relBinning[] = { 0., 0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.004, 0.006, 0.01, 0.015, 0.02, 0.03, 0.05, 0.1, 1.0 };
		const double absBinning[] = { 0., 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 100.0, 1000., 10000. };

		const int relNBins = sizeof(relBinning)/sizeof(relBinning[0]) - 1;
		const int absNBins = sizeof(absBinning)/sizeof(absBinning[0]) - 1;

		radAbs = new TH1D("radAbs", "abs;E;Entries", absNBins, absBinning);
		radPosAbs = new TH1D("radPosAbs", "abs;E;Entries", absNBins, absBinning);

		dir->cd();
	}

	void fill(const Tau& posTau, const Tau& negTau, float weight)
	{
		const double Einitial = posTau.p4genInitial.E() + negTau.p4genInitial.E();
		const double Efinal = posTau.p4genFinal.E() + negTau.p4genFinal.E();
		const double Eloss = Einitial - Efinal;
		radAbs->Fill(std::max(0., Eloss), weight);
		radPosAbs->Fill(std::max(0., posTau.p4genInitial.E() - posTau.p4genFinal.E()));
	}

	TH1D* radAbs;
	TH1D* radPosAbs;
};

struct Histograms {
	Histograms(TDirectory* dir, const char* name, bool is_embedded) {
		TDirectory* newdir = dir->mkdir(name);
		newdir->cd();

		posInitialPt = new TH1D("posInitialPt", "posInitialPt;p_{T} [GeV];Entries", 200, 0., 200.);
		posFinalPt = new TH1D("posFinalPt", "posFinalPt;p_{T} [GeV];Entries", 200, 0., 200.);

		x1x2 = new TH2D("x1x2", "x1x2;x_{1};x_{2}", 102, -0.01, 1.01, 102, -0.01, 1.01);
		vismass = new TH1D("vismass", "Visible mass;m_{vis};Entries/GeV", 150, 0., 150.);

		zs_rf = new TH1D("zs_rf", "zs_rf;z_{s};Entries/0.01", 100, -0.5, 0.5);
		x1_rf = new TH1D("x1_rf", "x1_rf;x_{1};Entries/0.01", 100, 0., 1.);
		x1_rf_h = new TH1D("x1_rf_h", "x1_rf_h;x_{1};Entries/0.01", 100, 0., 1.);

		x1_rf_80_85 = new TH1D("x1_rf_80_85", "x1_rf;x_{1};Entries/0.01", 100, 0., 1.);
		x1_rf_85_90 = new TH1D("x1_rf_85_90", "x1_rf;x_{1};Entries/0.01", 100, 0., 1.);
		x1_rf_90_95 = new TH1D("x1_rf_90_95", "x1_rf;x_{1};Entries/0.01", 100, 0., 1.);
		x1_rf_95_100 = new TH1D("x1_rf_95_100", "x1_rf;x_{1};Entries/0.01", 100, 0., 1.);
		x1_rf_100_105 = new TH1D("x1_rf_100_105", "x1_rf;x_{1};Entries/0.01", 100, 0., 1.);

		dir->cd();
	}

	double calczs(double xplus, double xminus)
	{
		const double a = xplus - xminus;
		if(a > 0.) return 0.5*(1. - (1.-a)*(1.-a));
		else return 0.5*((1.+a)*(1.+a) - 1.);
	}

	void fill(const Tau& posTau, const Tau& negTau, const TreeVars& vars, float weight)
	{
		std::pair<double, double> rfx = xInRF(posTau.p4genVis, posTau.p4genInitial, negTau.p4genVis, negTau.p4genInitial);

		if(posTau.hasVis)
		{
			posInitialPt->Fill(posTau.p4genInitial.Pt());
			posFinalPt->Fill(posTau.p4genFinal.Pt());
		}

		if(negTau.hasVis && posTau.hasVis)
		{
			x1x2->Fill(rfx.first, rfx.second, weight);
			vismass->Fill( (posTau.p4genVis + negTau.p4genVis).M(), weight);
			zs_rf->Fill(calczs(rfx.first, rfx.second), weight);
		}

		// TODO: Fill the same histos for negTau to double statistics
		if(posTau.hasVis)
		{
			const double M = (posTau.p4genInitial + negTau.p4genInitial).M();
			if(fabs((posTau.p4genInitial + negTau.p4genInitial).M() - 91.19) < 3.)
				x1_rf->Fill(rfx.first, weight);
			if(fabs((posTau.p4genInitial + negTau.p4genInitial).M() - 125.) < 3.)
				x1_rf_h->Fill(rfx.first, weight);

			if(M >= 80. && M <= 85.) x1_rf_80_85->Fill(rfx.first, weight);
			if(M >= 85. && M <= 90.) x1_rf_85_90->Fill(rfx.first, weight);
			if(M >= 90. && M <= 95.) x1_rf_90_95->Fill(rfx.first, weight);
			if(M >= 95. && M <= 100.) x1_rf_95_100->Fill(rfx.first, weight);
			if(M >= 100. && M <= 105.) x1_rf_100_105->Fill(rfx.first, weight);
		}

		if(negTau.hasVis)
		{
			const double M = (posTau.p4genInitial + negTau.p4genInitial).M();
			if(fabs((posTau.p4genInitial + negTau.p4genInitial).M() - 91.19) < 3.)
				x1_rf->Fill(rfx.second, weight);
			if(fabs((posTau.p4genInitial + negTau.p4genInitial).M() - 125.) < 3.)
				x1_rf_h->Fill(rfx.second, weight);

			if(M >= 80. && M <= 85.) x1_rf_80_85->Fill(rfx.second, weight);
			if(M >= 85. && M <= 90.) x1_rf_85_90->Fill(rfx.second, weight);
			if(M >= 90. && M <= 95.) x1_rf_90_95->Fill(rfx.second, weight);
			if(M >= 95. && M <= 100.) x1_rf_95_100->Fill(rfx.second, weight);
			if(M >= 100. && M <= 105.) x1_rf_100_105->Fill(rfx.second, weight);
		}
	}

	TH1D* posInitialPt;
	TH1D* posFinalPt;

	TH2D* x1x2;
	TH1D* vismass;

	TH1D* zs_rf;
	TH1D* x1_rf;
	TH1D* x1_rf_h;

	TH1D* x1_rf_80_85;
	TH1D* x1_rf_85_90;
	TH1D* x1_rf_90_95;
	TH1D* x1_rf_95_100;
	TH1D* x1_rf_100_105;
};

struct Category
{
	Category(TDirectory* dir, const char* name, bool is_embedded):
		isEmbedded(is_embedded)
	{
		TDirectory* newdir = dir->mkdir(name);
		newdir->cd();

		all.reset(new Histograms(newdir, "all", is_embedded));
		embAccepted.reset(new Histograms(newdir, "embAccepted", is_embedded));
		recAccepted.reset(new Histograms(newdir, "recAccepted", is_embedded));

		dir->cd();
	}

	void fill(const Tau& posTau, const Tau& negTau, const TreeVars& vars, float weight)
	{
		const Acceptance acceptance(posTau, negTau);
		all->fill(posTau, negTau, vars, weight);
		if(acceptance.embAccepted) embAccepted->fill(posTau, negTau, vars, weight);
		if(acceptance.recAccepted) recAccepted->fill(posTau, negTau, vars, weight);
	}

	const bool isEmbedded;

	std::auto_ptr<Histograms> all;
	std::auto_ptr<Histograms> embAccepted;
	std::auto_ptr<Histograms> recAccepted;
};

struct Out
{
	Out(TDirectory* dir, bool is_embedded):
		isEmbedded(is_embedded),
		radiation(dir, "radiation"),
		inclusive(dir, "inclusive", is_embedded)
	{
		inclusiveNoRadiation.reset(new Category(dir, "inclusiveNoRadiation", is_embedded));
		inclusiveNoOrigRadiation.reset(new Category(dir, "inclusiveNoOrigRadiation", is_embedded));
		inclusiveNoGenRadiation.reset(new Category(dir, "inclusiveNoGenRadiation", is_embedded));

		if(is_embedded)
		{
			inclusiveSpinned.reset(new Category(dir, "inclusiveSpinned", is_embedded));
			inclusiveSpinnedNoRadiation.reset(new Category(dir, "inclusiveSpinnedNoRadiation", is_embedded));
			inclusiveSpinnedNoOrigRadiation.reset(new Category(dir, "inclusiveSpinnedNoOrigRadiation", is_embedded));
			inclusiveSpinnedNoGenRadiation.reset(new Category(dir, "inclusiveSpinnedNoGenRadiation", is_embedded));
		}

		nEvents = new TH1D("nEvents", "Total number of events", 1, 0, 1);
		nTrueInteractions = new TH1D("nTrueInteractions", "Number of True Interactions", 100, -0.5, 99.5);
	}

	bool isOrigRadiation(const Tau& posTau, const Tau& negTau)
	{
		if(!isEmbedded) return false;

		// Just to be consistent with what is done below:
		if(posTau.p4origGenInitial.E() - posTau.p4origGenFinal.E() < 0 && negTau.p4origGenInitial.E() - negTau.p4origGenFinal.E() < 0)
			return false;

		if(!lorentzEqual(posTau.p4origGenInitial, posTau.p4origGenFinal) || !lorentzEqual(negTau.p4origGenInitial, negTau.p4origGenFinal)) return true;
		return false;
	}

	bool isGenRadiation(const Tau& posTau, const Tau& negTau)
	{
		// For some reason, for the non-radiation case, this often becomes negative. Maybe there is a problem in the skimming code?
		if(posTau.p4genInitial.E() - posTau.p4genFinal.E() < 0 && negTau.p4genInitial.E() - negTau.p4genFinal.E() < 0)
			return false;

		if(!lorentzEqual(posTau.p4genInitial, posTau.p4genFinal) || !lorentzEqual(negTau.p4genInitial, negTau.p4genFinal))
			return true;
		return false;
	}

	void fill(const Tau& posTau, const Tau& negTau, const TreeVars& vars, float weight)
	{
		if(abs(vars.PosTau.GenVisPDG) != 211 || abs(vars.NegTau.GenVisPDG) != 211) return;

		const bool origRadiation = isOrigRadiation(posTau, negTau);
		const bool genRadiation = isGenRadiation(posTau, negTau);
		const bool radiation = origRadiation || genRadiation;

		this->radiation.fill(posTau, negTau, weight);
		inclusive.fill(posTau, negTau, vars, weight);
		if(!radiation && inclusiveNoRadiation.get()) inclusiveNoRadiation->fill(posTau, negTau, vars, weight);
		if(!origRadiation && inclusiveNoOrigRadiation.get()) inclusiveNoOrigRadiation->fill(posTau, negTau, vars, weight);
		if(!genRadiation && inclusiveNoGenRadiation.get()) inclusiveNoGenRadiation->fill(posTau, negTau, vars, weight);

		if(inclusiveSpinned.get()) inclusiveSpinned->fill(posTau, negTau, vars, weight*vars.tauSpinWeight);
		if(!radiation && inclusiveSpinnedNoRadiation.get()) inclusiveSpinnedNoRadiation->fill(posTau, negTau, vars, weight*vars.tauSpinWeight);
		if(!origRadiation && inclusiveSpinnedNoOrigRadiation.get()) inclusiveSpinnedNoOrigRadiation->fill(posTau, negTau, vars, weight*vars.tauSpinWeight);
		if(!genRadiation && inclusiveSpinnedNoGenRadiation.get()) inclusiveSpinnedNoGenRadiation->fill(posTau, negTau, vars, weight*vars.tauSpinWeight);
	}

	void add_events(unsigned int n)
	{
		nEvents->SetBinContent(1, static_cast<double>(nEvents->GetBinContent(1) + n));
	}

	void add_interactions(TH1F* nTrueInteractions)
	{
		this->nTrueInteractions->Add(nTrueInteractions);
	}

	const bool isEmbedded;

	RadiationHistograms radiation;
	Category inclusive;
	std::auto_ptr<Category> inclusiveNoRadiation;
	std::auto_ptr<Category> inclusiveNoOrigRadiation;
	std::auto_ptr<Category> inclusiveNoGenRadiation;
	std::auto_ptr<Category> inclusiveSpinned;
	std::auto_ptr<Category> inclusiveSpinnedNoRadiation;
	std::auto_ptr<Category> inclusiveSpinnedNoOrigRadiation;
	std::auto_ptr<Category> inclusiveSpinnedNoGenRadiation;

	TH1D* nEvents;
	TH1D* nTrueInteractions;
};

void process_file(const char* filename, Out& out, bool is_embedded, const ZmumuEvtSelEffCorrWeight* zmumuWeight)
{
	TreeVars vars;

	std::string directory = "spinTest";
	TFile* file = new TFile(filename, "READ");
	TTree* tree = dynamic_cast<TTree*>(file->Get((directory + "/TauTree").c_str()));

	tree->SetBranchAddress("Run", &vars.run);
	tree->SetBranchAddress("Lumi", &vars.lumi);
	tree->SetBranchAddress("Event", &vars.event);

	if(is_embedded)
	{
		tree->SetBranchAddress("MinVisPtFilterWeight", &vars.minVisPtFilterWeight);
		tree->SetBranchAddress("TauSpinWeight", &vars.tauSpinWeight);
		tree->SetBranchAddress("ZmumuEvtSelEffCorrWeight", &vars.zmumuEvtSelEffCorrWeight);

		tree->SetBranchAddress("OrigGenPosChargedLepPt", &vars.PosTau.OrigGenPt);
		tree->SetBranchAddress("OrigGenPosChargedLepEta", &vars.PosTau.OrigGenEta);
		tree->SetBranchAddress("OrigGenPosChargedLepPhi", &vars.PosTau.OrigGenPhi);
		tree->SetBranchAddress("OrigGenPosChargedLepE", &vars.PosTau.OrigGenE);
		tree->SetBranchAddress("OrigGenPosChargedLepFinalPt", &vars.PosTau.OrigGenFinalPt);
		tree->SetBranchAddress("OrigGenPosChargedLepFinalEta", &vars.PosTau.OrigGenFinalEta);
		tree->SetBranchAddress("OrigGenPosChargedLepFinalPhi", &vars.PosTau.OrigGenFinalPhi);
		tree->SetBranchAddress("OrigGenPosChargedLepFinalE", &vars.PosTau.OrigGenFinalE);
		tree->SetBranchAddress("OrigGenPosChargedLepPDG", &vars.PosTau.OrigGenChargedLepPDG);

		tree->SetBranchAddress("OrigGenNegChargedLepPt", &vars.NegTau.OrigGenPt);
		tree->SetBranchAddress("OrigGenNegChargedLepEta", &vars.NegTau.OrigGenEta);
		tree->SetBranchAddress("OrigGenNegChargedLepPhi", &vars.NegTau.OrigGenPhi);
		tree->SetBranchAddress("OrigGenNegChargedLepE", &vars.NegTau.OrigGenE);
		tree->SetBranchAddress("OrigGenNegChargedLepFinalPt", &vars.NegTau.OrigGenFinalPt);
		tree->SetBranchAddress("OrigGenNegChargedLepFinalEta", &vars.NegTau.OrigGenFinalEta);
		tree->SetBranchAddress("OrigGenNegChargedLepFinalPhi", &vars.NegTau.OrigGenFinalPhi);
		tree->SetBranchAddress("OrigGenNegChargedLepFinalE", &vars.NegTau.OrigGenFinalE);
		tree->SetBranchAddress("OrigGenNegChargedLepPDG", &vars.NegTau.OrigGenChargedLepPDG);
	}

	tree->SetBranchAddress("NTrueInteractions", &vars.nTrueInteractions);

	tree->SetBranchAddress("GenPosTauPt", &vars.PosTau.GenPt);
	tree->SetBranchAddress("GenPosTauEta", &vars.PosTau.GenEta);
	tree->SetBranchAddress("GenPosTauPhi", &vars.PosTau.GenPhi);
	tree->SetBranchAddress("GenPosTauE", &vars.PosTau.GenE);
	tree->SetBranchAddress("GenPosTauFinalPt", &vars.PosTau.GenFinalPt);
	tree->SetBranchAddress("GenPosTauFinalEta", &vars.PosTau.GenFinalEta);
	tree->SetBranchAddress("GenPosTauFinalPhi", &vars.PosTau.GenFinalPhi);
	tree->SetBranchAddress("GenPosTauFinalE", &vars.PosTau.GenFinalE);
	tree->SetBranchAddress("GenPosTauVisPt", &vars.PosTau.GenVisPt);
	tree->SetBranchAddress("GenPosTauVisEta", &vars.PosTau.GenVisEta);
	tree->SetBranchAddress("GenPosTauVisPhi", &vars.PosTau.GenVisPhi);
	tree->SetBranchAddress("GenPosTauVisE", &vars.PosTau.GenVisE);
	tree->SetBranchAddress("GenPosTauDecayPDG", &vars.PosTau.GenVisPDG);

	tree->SetBranchAddress("GenNegTauPt", &vars.NegTau.GenPt);
	tree->SetBranchAddress("GenNegTauEta", &vars.NegTau.GenEta);
	tree->SetBranchAddress("GenNegTauPhi", &vars.NegTau.GenPhi);
	tree->SetBranchAddress("GenNegTauE", &vars.NegTau.GenE);
	tree->SetBranchAddress("GenNegTauFinalPt", &vars.NegTau.GenFinalPt);
	tree->SetBranchAddress("GenNegTauFinalEta", &vars.NegTau.GenFinalEta);
	tree->SetBranchAddress("GenNegTauFinalPhi", &vars.NegTau.GenFinalPhi);
	tree->SetBranchAddress("GenNegTauFinalE", &vars.NegTau.GenFinalE);
	tree->SetBranchAddress("GenNegTauVisPt", &vars.NegTau.GenVisPt);
	tree->SetBranchAddress("GenNegTauVisEta", &vars.NegTau.GenVisEta);
	tree->SetBranchAddress("GenNegTauVisPhi", &vars.NegTau.GenVisPhi);
	tree->SetBranchAddress("GenNegTauVisE", &vars.NegTau.GenVisE);
	tree->SetBranchAddress("GenNegTauDecayPDG", &vars.NegTau.GenVisPDG);

	TH1F* h_nTrueInteractions = dynamic_cast<TH1F*>(file->Get((directory + "/h_nTrueInteractions").c_str()));
	out.add_events(static_cast<unsigned int>(h_nTrueInteractions->GetEntries()));
	out.add_interactions(h_nTrueInteractions);

	const unsigned int CountTotal = tree->GetEntries();
	for(int i = 0; i < CountTotal; ++i)
	{
		tree->GetEntry(i);

		// PDG Filter
		if(is_embedded)
			if(abs(vars.PosTau.OrigGenChargedLepPDG) != 13 || abs(vars.NegTau.OrigGenChargedLepPDG) != 13)
				continue;

		// Taus
		Tau posTau(vars.PosTau, +1);
		Tau negTau(vars.NegTau, -1);

		double weight = 1.0f;
		if(zmumuWeight) weight *= zmumuWeight->getWeight(posTau.p4genInitial, negTau.p4genInitial);

		// Output
		out.fill(posTau, negTau, vars, weight);
	}

	delete file;
}

void process_dataset(const char* dataset, const std::vector<std::string>& files)
{
	const bool is_embedded = (strstr(dataset, "Embed") != NULL);
	const bool is_gen_embedded = (strstr(dataset, "GenEmbed") != NULL);

	std::stringstream out_filename_stream;
	out_filename_stream << "histos-" << dataset;
	out_filename_stream << ".root";
	const std::string out_filename = out_filename_stream.str();

	std::cout << dataset << "... => " << out_filename << std::endl;
	std::cout << "\tEmbedded: " << is_embedded << std::endl;
	std::cout << "\tGenEmbedded: " << is_gen_embedded << std::endl;

	std::auto_ptr<ZmumuEvtSelEffCorrWeight> zmumuWeight;
	if(is_embedded && !is_gen_embedded) zmumuWeight.reset(new ZmumuEvtSelEffCorrWeight("../data/ZmumuEvtSelEffCorrWeights.root"));

	TFile* outfile = new TFile(out_filename.c_str(), "RECREATE");
	Out out(outfile, is_embedded);

	for(std::vector<std::string>::const_iterator iter = files.begin(); iter != files.end(); ++iter)
	{
		std::cout << "\t[" << (iter - files.begin() + 1) << "/" << files.size() << "] " << *iter << "..." << std::endl;
		process_file(iter->c_str(), out, is_embedded, zmumuWeight.get());
	}

	outfile->Write();
	delete outfile;
}

int main(int argc, char* argv[])
{
	TH1::SetDefaultSumw2(true);

	std::vector<std::string> files;
	for(int i = 1; i < argc; ++i)
	{
		files.push_back(argv[i]);
	}

	typedef std::map<std::string, std::vector<std::string> > DatasetMap;
	DatasetMap datasets;
	for(std::vector<std::string>::const_iterator iter = files.begin(); iter != files.end(); ++iter)
	{
		std::string::size_type ipos = iter->rfind('/');
		if(ipos == std::string::npos) ipos = 0; else ++ipos;

		std::string basename = iter->substr(ipos);
		if(basename.rfind(".root") != std::string::npos)
			basename.erase(basename.length() - 5 , 5);

		std::string dataset = basename;
		std::string::size_type epos = basename.rfind('_');
		if(epos != std::string::npos)
			dataset = basename.substr(0, epos);

		datasets[dataset].push_back(*iter);
	}

	std::cout << "Datasets:" << std::endl;
	for(DatasetMap::const_iterator iter = datasets.begin(); iter != datasets.end(); ++iter)
		std::cout << "\t" << iter->first << std::endl;

	std::cout << std::endl;
	for(DatasetMap::const_iterator iter = datasets.begin(); iter != datasets.end(); ++iter)
		process_dataset(iter->first.c_str(), iter->second);

	return 0;
}
