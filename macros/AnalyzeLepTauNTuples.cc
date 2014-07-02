#include <TTree.h>
#include <TFile.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
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

static double calcMt(double px1, double py1, double px2, double py2)
{
	double pt1 = TMath::Sqrt(px1*px1 + py1*py1);
	double pt2 = TMath::Sqrt(px2*px2 + py2*py2);

	double p1Dotp2 = px1*px2 + py1*py2;
	double cosAlpha = p1Dotp2/(pt1*pt2);

	return TMath::Sqrt(2*pt1*pt2*(1 - cosAlpha));
}

static int findBin(TAxis* axis, double x)
{
	int bin = axis->FindBin(x);
	if(bin < 1) bin = 1;
	if(bin > axis->GetNbins()) bin = axis->GetNbins();
	return bin;
}

struct TreeVars {
	TreeVars() {
		run = lumi = event = ~0u;
		minVisPtFilterWeight = -1.0f;
		tauSpinWeight = -1.0f;
		muonRadiationFilter = true;
		muonRadiationFilter2Sel1 = muonRadiationFilter2Sel2 = muonRadiationFilter2Sel3 = true;

		nPV = 0; nTrueInteractions = -1.0f;

		origHltMu17Mu8 = false;
		origPosMuonPt = origPosMuonEta = origPosMuonPhi = origPosMuonE = -1.0f;
		origGenPosMuonPt = origGenPosMuonEta = origGenPosMuonPhi = origGenPosMuonE = -1.0f;
		origGenPosMuonFinalPt = origGenPosMuonFinalEta = origGenPosMuonFinalPhi = origGenPosMuonFinalE = -1.0f;
		origNegMuonPt = origNegMuonEta = origNegMuonPhi = origNegMuonE = -1.0f;
		origGenNegMuonPt = origGenNegMuonEta = origGenNegMuonPhi = origGenNegMuonE = -1.0f;
		origGenNegMuonFinalPt = origGenNegMuonFinalEta = origGenNegMuonFinalPhi = origGenNegMuonFinalE = -1.0f;
		origGenPosMuonPDG = origGenNegMuonPDG = 0;

		caloMet = caloMetPhi = -1.0f;
		pfMet = pfMetPhi = -1.0f;
		mvaPfMet = mvaPfMetPhi = -1.0f;

		electronPt = electronEta = electronPhi = electronE = -1.0f;
		electronCharge = 0;
		electronChargedParticlePfIso04 = electronChargedHadronPfIso04 = electronNeutralHadronPfIso04 = electronPhotonPfIso04 = electronPUPfIso04 = -1.0f;
		electronMvaId = -1.0f;
		electronQuality = false;
		electronHasConversion = false;
		electronGenPt = electronGenEta = electronGenPhi = electronGenE = -1.0f;
		electronGenFinalPt = electronGenFinalEta = electronGenFinalPhi = electronGenFinalE = -1.0f;
		electronGenCharge = 0;
		electronGenPtVis = electronGenEtaVis = electronGenPhiVis = electronGenEVis = -1.0f;

		muonPt = muonEta = muonPhi = muonE = -1.0f;
		muonCharge = 0;
		muonChargedParticlePfIso04 = muonChargedHadronPfIso04 = muonNeutralHadronPfIso04 = muonPhotonPfIso04 = muonPUPfIso04 = -1.0f;
		muonQuality = false;
		muonGenPt = muonGenEta = muonGenPhi = muonGenE = -1.0f;
		muonGenFinalPt = muonGenFinalEta = muonGenFinalPhi = muonGenFinalE = -1.0f;
		muonGenCharge = 0;
		muonGenPtVis = muonGenEtaVis = muonGenPhiVis = muonGenEVis = -1.0f;

		tauPt = tauEta = tauPhi = tauE = -1.0f;
		tauCharge = 0;
		tauDecayModeFinding = false;
		tauQuality = false;
		tauIsolation3Hits = -1.0f;
		tauLooseIsolation3Hits = tauMediumIsolation3Hits = tauTightIsolation3Hits = false;
		tauAgainstMuonLoose = tauAgainstMuonTight = false;
		tauAgainstElectronLoose = tauAgainstElectronLooseMVA3 = tauAgainstElectronMediumMVA3 = tauAgainstElectronTightMVA3 = false;
		tauGenPt = tauGenEta = tauGenPhi = tauGenE = -1.0f;
		tauGenFinalPt = tauGenFinalEta = tauGenFinalPhi = tauGenFinalE = -1.0f;
		tauGenCharge = 0;
		tauGenPtVis = tauGenEtaVis = tauGenPhiVis = tauGenEVis = -1.0f;

		nJets = 0;

		svfitMass = 0.0f;
	}

	unsigned int run;
	unsigned int lumi;
	unsigned int event;

	float minVisPtFilterWeight;
	float tauSpinWeight;
	bool muonRadiationFilter;
	bool muonRadiationFilter2Sel1;
	bool muonRadiationFilter2Sel2;
	bool muonRadiationFilter2Sel3;

	unsigned int nPV;
	float nTrueInteractions;

	bool origHltMu17Mu8;
	float origPosMuonPt;
	float origPosMuonEta;
	float origPosMuonPhi;
	float origPosMuonE;
	float origGenPosMuonPt;
	float origGenPosMuonEta;
	float origGenPosMuonPhi;
	float origGenPosMuonE;
	float origGenPosMuonFinalPt;
	float origGenPosMuonFinalEta;
	float origGenPosMuonFinalPhi;
	float origGenPosMuonFinalE;
	int origGenPosMuonPDG;
	float origNegMuonPt;
	float origNegMuonEta;
	float origNegMuonPhi;
	float origNegMuonE;
	float origGenNegMuonPt;
	float origGenNegMuonEta;
	float origGenNegMuonPhi;
	float origGenNegMuonE;
	float origGenNegMuonFinalPt;
	float origGenNegMuonFinalEta;
	float origGenNegMuonFinalPhi;
	float origGenNegMuonFinalE;
	int origGenNegMuonPDG;

	float caloMet;
	float caloMetPhi;
	float pfMet;
	float pfMetPhi;
	float mvaPfMet;
	float mvaPfMetPhi;

	float electronPt;
	float electronEta;
	float electronPhi;
	float electronE;
	int electronCharge;
	float electronChargedParticlePfIso04;
	float electronChargedHadronPfIso04;
	float electronNeutralHadronPfIso04;
	float electronPhotonPfIso04;
	float electronPUPfIso04;
	float electronMvaId;
	bool electronQuality;
	bool electronHasConversion;
	float electronGenPt;
	float electronGenEta;
	float electronGenPhi;
	float electronGenE;
	float electronGenFinalPt;
	float electronGenFinalEta;
	float electronGenFinalPhi;
	float electronGenFinalE;
	int electronGenCharge;
	float electronGenPtVis;
	float electronGenEtaVis;
	float electronGenPhiVis;
	float electronGenEVis;

	float muonPt;
	float muonEta;
	float muonPhi;
	float muonE;
	int muonCharge;
	float muonChargedParticlePfIso04;
	float muonChargedHadronPfIso04;
	float muonNeutralHadronPfIso04;
	float muonPhotonPfIso04;
	float muonPUPfIso04;
	bool muonQuality;
	float muonGenPt;
	float muonGenEta;
	float muonGenPhi;
	float muonGenE;
	float muonGenFinalPt;
	float muonGenFinalEta;
	float muonGenFinalPhi;
	float muonGenFinalE;
	int muonGenCharge;
	float muonGenPtVis;
	float muonGenEtaVis;
	float muonGenPhiVis;
	float muonGenEVis;

	float tauPt;
	float tauEta;
	float tauPhi;
	float tauE;
	int tauCharge;
	bool tauDecayModeFinding;
	bool tauQuality;
	float tauIsolation3Hits;
	bool tauLooseIsolation3Hits;
	bool tauMediumIsolation3Hits;
	bool tauTightIsolation3Hits;
	bool tauAgainstMuonLoose;
	bool tauAgainstMuonTight;
	bool tauAgainstElectronLoose;
	bool tauAgainstElectronLooseMVA3;
	bool tauAgainstElectronMediumMVA3;
	bool tauAgainstElectronTightMVA3;
	float tauGenPt;
	float tauGenEta;
	float tauGenPhi;
	float tauGenE;
	float tauGenFinalPt;
	float tauGenFinalEta;
	float tauGenFinalPhi;
	float tauGenFinalE;
	int tauGenCharge;
	float tauGenPtVis;
	float tauGenEtaVis;
	float tauGenPhiVis;
	float tauGenEVis;

	unsigned int nJets;
	float jetPt[20];
	float jetEta[20];
	float jetPhi[20];
	float jetEnergy[20];
	float jetBTag[20];

	float svfitMass;
};

TLorentzVector getOrigVector(const TreeVars& vars, int charge)
{
	if(charge > 0)
		return TLorentzVector_PtEtaPhiE(vars.origPosMuonPt, vars.origPosMuonEta, vars.origPosMuonPhi, vars.origPosMuonE);
	else
		return TLorentzVector_PtEtaPhiE(vars.origNegMuonPt, vars.origNegMuonEta, vars.origNegMuonPhi, vars.origNegMuonE);
}

TLorentzVector getOrigGenInitialVector(const TreeVars& vars, int charge)
{
	if(charge > 0)
		return TLorentzVector_PtEtaPhiE(vars.origGenPosMuonPt, vars.origGenPosMuonEta, vars.origGenPosMuonPhi, vars.origGenPosMuonE);
	else
		return TLorentzVector_PtEtaPhiE(vars.origGenNegMuonPt, vars.origGenNegMuonEta, vars.origGenNegMuonPhi, vars.origGenNegMuonE);
}

TLorentzVector getOrigGenFinalVector(const TreeVars& vars, int charge)
{
	if(charge > 0)
		return TLorentzVector_PtEtaPhiE(vars.origGenPosMuonFinalPt, vars.origGenPosMuonFinalEta, vars.origGenPosMuonFinalPhi, vars.origGenPosMuonFinalE);
	else
		return TLorentzVector_PtEtaPhiE(vars.origGenNegMuonFinalPt, vars.origGenNegMuonFinalEta, vars.origGenNegMuonFinalPhi, vars.origGenNegMuonFinalE);
}

class Electron
{
public:
	Electron(const TreeVars& vars):
		hasGen(vars.electronGenPt > 0.),
		hasRec(vars.electronPt > 0.),
		hasOrigGen(vars.electronCharge > 0 ? (vars.origGenPosMuonPt > 0.) : (vars.origGenNegMuonPt > 0.)),
		hasOrigRec(vars.electronCharge > 0 ? (vars.origPosMuonPt > 0.) : (vars.origNegMuonPt > 0.)),
		p4(TLorentzVector_PtEtaPhiE(vars.electronPt, vars.electronEta, vars.electronPhi, vars.electronE)),
		p4genInitial(TLorentzVector_PtEtaPhiE(vars.electronGenPt, vars.electronGenEta, vars.electronGenPhi, vars.electronGenE)),
		p4genFinal(TLorentzVector_PtEtaPhiE(vars.electronGenFinalPt, vars.electronGenFinalEta, vars.electronGenFinalPhi, vars.electronGenFinalE)),
		p4genVis(TLorentzVector_PtEtaPhiE(vars.electronGenPtVis, vars.electronGenEtaVis, vars.electronGenPhiVis, vars.electronGenEVis)),
		p4orig(getOrigVector(vars, vars.electronGenCharge)),
		p4origGenInitial(getOrigGenInitialVector(vars, vars.electronGenCharge)),
		p4origGenFinal(getOrigGenFinalVector(vars, vars.electronGenCharge)),
		genCharge(vars.electronGenCharge),
		charge(vars.electronCharge),
		quality(vars.electronQuality),
		hasConversion(vars.electronHasConversion),
		mvaId(vars.electronMvaId),
		chargedParticleIsolation(vars.electronChargedParticlePfIso04/vars.electronPt),
		chargedHadronIsolation(vars.electronChargedHadronPfIso04/vars.electronPt),
		neutralHadronIsolation(vars.electronNeutralHadronPfIso04/vars.electronPt),
		photonIsolation(vars.electronPhotonPfIso04/vars.electronPt),
		puIsolation(vars.electronPUPfIso04/vars.electronPt),
		isolation((vars.electronChargedParticlePfIso04 + std::max(0., vars.electronNeutralHadronPfIso04 + vars.electronPhotonPfIso04 - 0.5 * vars.electronPUPfIso04))/vars.electronPt)
	{
	}

	bool in_acceptance(const TLorentzVector& p4) const { return p4.Pt() > 22. && fabs(p4.Eta()) < 2.1; }
	bool isolated() const { return isolation < (p4.Eta() < 1.479 ? 0.15 : 0.10); }

	const bool looseID() const
	{
		const double pt = p4.Pt();
		const double eta = fabs(p4.Eta());

		if(pt < 20. && eta < 0.8) return mvaId > 0.925;
		else if(pt < 20. && eta >= 0.8 && eta < 1.479) return mvaId > 0.915;
		else if(pt < 20. && eta >= 1.479) return mvaId > 0.965;
		else if(pt >= 20. && eta < 0.8) return mvaId > 0.905;
		else if(pt >= 20. && eta >= 0.8 && eta < 1.479) return mvaId > 0.955;
		else if(pt >= 20. && eta >= 1.479) return mvaId > 0.975;
		else { assert(false); return false; }
	}

	const bool id() const { return looseID() && !hasConversion; }

	const bool hasGen;
	const bool hasRec;
	const bool hasOrigGen;
	const bool hasOrigRec;
	const TLorentzVector p4;
	const TLorentzVector p4genInitial;
	const TLorentzVector p4genFinal;
	const TLorentzVector p4genVis;
	const TLorentzVector p4orig;
	const TLorentzVector p4origGenInitial;
	const TLorentzVector p4origGenFinal;
	const int genCharge;
	const int charge;
	const bool quality;
	const bool hasConversion;
	const float mvaId;
	const float chargedParticleIsolation;
	const float chargedHadronIsolation;
	const float neutralHadronIsolation;
	const float photonIsolation;
	const float puIsolation;
	const float isolation;
};

class Muon
{
public:
	Muon(const TreeVars& vars):
		hasGen(vars.muonGenPt > 0.),
		hasRec(vars.muonPt > 0.),
		hasOrigGen(vars.muonCharge > 0 ? (vars.origGenPosMuonPt > 0.) : (vars.origGenNegMuonPt > 0.)),
		hasOrigRec(vars.muonCharge > 0 ? (vars.origPosMuonPt > 0.) : (vars.origNegMuonPt > 0.)),
		p4(TLorentzVector_PtEtaPhiE(vars.muonPt, vars.muonEta, vars.muonPhi, vars.muonE)),
		p4genInitial(TLorentzVector_PtEtaPhiE(vars.muonGenPt, vars.muonGenEta, vars.muonGenPhi, vars.muonGenE)),
		p4genFinal(TLorentzVector_PtEtaPhiE(vars.muonGenFinalPt, vars.muonGenFinalEta, vars.muonGenFinalPhi, vars.muonGenFinalE)),
		p4genVis(TLorentzVector_PtEtaPhiE(vars.muonGenPtVis, vars.muonGenEtaVis, vars.muonGenPhiVis, vars.muonGenEVis)),
		p4orig(getOrigVector(vars, vars.muonGenCharge)),
		p4origGenInitial(getOrigGenInitialVector(vars, vars.muonGenCharge)),
		p4origGenFinal(getOrigGenFinalVector(vars, vars.muonGenCharge)),
		genCharge(vars.muonGenCharge),
		charge(vars.muonCharge),
		quality(vars.muonQuality),
		chargedParticleIsolation(vars.muonChargedParticlePfIso04/vars.muonPt),
		chargedHadronIsolation(vars.muonChargedHadronPfIso04/vars.muonPt),
		neutralHadronIsolation(vars.muonNeutralHadronPfIso04/vars.muonPt),
		photonIsolation(vars.muonPhotonPfIso04/vars.muonPt),
		puIsolation(vars.muonPUPfIso04/vars.muonPt),
		isolation((vars.muonChargedParticlePfIso04 + std::max(0., vars.muonNeutralHadronPfIso04 + vars.muonPhotonPfIso04 - 0.5 * vars.muonPUPfIso04))/vars.muonPt)
	{
	}

	bool in_acceptance(const TLorentzVector& p4) const { return p4.Pt() > 17. && fabs(p4.Eta()) < 2.1; }
	bool isolated() const { return isolation < 0.1; } 

	const bool id() const { return quality; }

	const bool hasGen;
	const bool hasRec;
	const bool hasOrigGen;
	const bool hasOrigRec;
	const TLorentzVector p4;
	const TLorentzVector p4genInitial;
	const TLorentzVector p4genFinal;
	const TLorentzVector p4genVis;
	const TLorentzVector p4orig;
	const TLorentzVector p4origGenInitial;
	const TLorentzVector p4origGenFinal;
	const int genCharge;
	const int charge;
	const bool quality;
	const float chargedParticleIsolation;
	const float chargedHadronIsolation;
	const float neutralHadronIsolation;
	const float photonIsolation;
	const float puIsolation;
	const float isolation;
};

class Tau
{
public:
	Tau(const TreeVars& vars):
		hasGen(vars.tauGenPt > 0.),
		hasRec(vars.tauPt > 0.),
		hasOrigGen(vars.tauCharge > 0 ? (vars.origGenPosMuonPt > 0.) : (vars.origGenNegMuonPt > 0.)),
		hasOrigRec(vars.tauCharge > 0 ? (vars.origPosMuonPt > 0.) : (vars.origNegMuonPt > 0.)),
		p4(TLorentzVector_PtEtaPhiE(vars.tauPt, vars.tauEta, vars.tauPhi, vars.tauE)),
		p4genInitial(TLorentzVector_PtEtaPhiE(vars.tauGenPt, vars.tauGenEta, vars.tauGenPhi, vars.tauGenE)),
		p4genFinal(TLorentzVector_PtEtaPhiE(vars.tauGenFinalPt, vars.tauGenFinalEta, vars.tauGenFinalPhi, vars.tauGenFinalE)),
		p4genVis(TLorentzVector_PtEtaPhiE(vars.tauGenPtVis, vars.tauGenEtaVis, vars.tauGenPhiVis, vars.tauGenEVis)),
		p4orig(getOrigVector(vars, vars.tauGenCharge)),
		p4origGenInitial(getOrigGenInitialVector(vars, vars.tauGenCharge)),
		p4origGenFinal(getOrigGenFinalVector(vars, vars.tauGenCharge)),
		genCharge(vars.tauGenCharge),
		charge(vars.tauCharge),
		decayModeFinding(vars.tauDecayModeFinding),
		quality(vars.tauQuality),
		isolation3Hits(vars.tauIsolation3Hits),
		looseIsolation3Hits(vars.tauLooseIsolation3Hits),
		tightIsolation3Hits(vars.tauTightIsolation3Hits),
		againstMuonTight(vars.tauAgainstMuonTight),
		againstElectronLoose(vars.tauAgainstElectronLoose)
	{
	}

	bool in_acceptance(const TLorentzVector& p4) const { return p4.Pt() > 20. && fabs(p4.Eta()) < 2.3; }

	const bool hasGen;
	const bool hasRec;
	const bool hasOrigGen;
	const bool hasOrigRec;
	const TLorentzVector p4;
	const TLorentzVector p4genInitial;
	const TLorentzVector p4genFinal;
	const TLorentzVector p4genVis;
	const TLorentzVector p4orig;
	const TLorentzVector p4origGenInitial;
	const TLorentzVector p4origGenFinal;
	const int genCharge;
	const int charge;
	const bool decayModeFinding;
	const bool quality;
	const float isolation3Hits;
	const bool looseIsolation3Hits;
	const bool tightIsolation3Hits;
	const bool againstMuonTight;
	const bool againstElectronLoose;
};

class ZmumuEvtSelEffCorrWeight
{
public:
	ZmumuEvtSelEffCorrWeight(const char* filename)
	{
		file = new TFile(filename, "READ");
		if(file->IsZombie()) throw std::runtime_error(("Could not load " + std::string(filename)).c_str());

		muMinusPt_vs_muPlusPt = dynamic_cast<TH2D*>(file->Get("EMB_ZmumuEvtSelEff_muMinusPt_vs_muPlusPt"));
		muMinusEta_vs_muPlusEta = dynamic_cast<TH2D*>(file->Get("EMB_ZmumuEvtSelEffCorr_muMinusEta_vs_muPlusEta"));
		if(!muMinusPt_vs_muPlusPt || !muMinusEta_vs_muPlusEta)
			throw std::runtime_error("Failed to load ZmumuEvtSelEff histograms");
	}

	template<typename TLepton>
	double getWeight(const TLepton& lepton, const Tau& tau) const
	{
		assert(lepton.hasOrigRec && tau.hasOrigRec);
		const TLorentzVector& plus = lepton.genCharge > 0 ? lepton.p4orig : tau.p4orig;
		const TLorentzVector& minus = lepton.genCharge < 0 ? lepton.p4orig : tau.p4orig;

		const int binPtX = findBin(muMinusPt_vs_muPlusPt->GetXaxis(), plus.Pt());
		const int binPtY = findBin(muMinusPt_vs_muPlusPt->GetYaxis(), minus.Pt());
		double effPt = muMinusPt_vs_muPlusPt->GetBinContent(binPtX, binPtY);

		const int binEtaX = findBin(muMinusEta_vs_muPlusEta->GetXaxis(), plus.Eta());
		const int binEtaY = findBin(muMinusEta_vs_muPlusEta->GetYaxis(), minus.Eta());
		double effEta = muMinusEta_vs_muPlusEta->GetBinContent(binEtaX, binEtaY);

		if(effPt < 1e-1) { return 1.; }//std::cout << "etaplus=" << plus.Eta() << ", etaminus=" << minus.Eta() << std::endl; return 1.; }
		if(effEta < 1e-1) { return 1.; }//std::cout << "etaplus=" << plus.Eta() << ", etaminus=" << minus.Eta() << std::endl; return 1.; }
		return 1./(effPt * effEta);
	}

private:
	TFile* file;
	TH2D* muMinusPt_vs_muPlusPt;
	TH2D* muMinusEta_vs_muPlusEta;
};

template<typename TLepton>
struct EventSelector
{
	bool acceptEmb(const TLorentzVector& first, const TLorentzVector& second)
	{
		// Embedding cuts (very loose). Choose 10/20 for pT since this is used to generate PF embedded.
		if(first.Pt() < 10 || second.Pt() < 10) return false;
		else if(first.Pt() < 20 && second.Pt() < 20) return false;
		else if(fabs(first.Eta()) > 2.5) return false;
		else if(fabs(second.Eta()) > 2.5) return false;
		else if((first + second).M() < 50.0) return false;
		return true;
	}

	EventSelector(const TreeVars& vars, const TLepton& lepton, const Tau& tau):
		origAccepted(lepton.hasOrigRec && tau.hasOrigRec && vars.origHltMu17Mu8 && acceptEmb(lepton.p4orig, tau.p4orig)),
		genInitialAccepted(lepton.hasGen && tau.hasGen && acceptEmb(lepton.p4genInitial, tau.p4genInitial) && fabs(lepton.p4genInitial.Eta()) < 2.4 && fabs(tau.p4genInitial.Eta()) < 2.4),
		genInitialMasscut((lepton.p4genInitial + tau.p4genInitial).M() > 50),
		genFinalAccepted(lepton.hasGen && tau.hasGen && acceptEmb(lepton.p4genFinal, tau.p4genFinal) && fabs(lepton.p4genFinal.Eta()) < 2.4 && fabs(tau.p4genFinal.Eta()) < 2.4),
		genFinalMasscut((lepton.p4genFinal + tau.p4genFinal).M() > 50),
		genVisAccepted(lepton.hasGen && tau.hasGen && lepton.in_acceptance(lepton.p4genVis) && tau.in_acceptance(tau.p4genVis)),
		lep_matched(lepton.hasGen && lepton.hasRec && ROOT::Math::VectorUtil::DeltaR(lepton.p4, lepton.p4genVis) < 0.3),
		tau_matched(tau.hasGen && tau.hasRec && ROOT::Math::VectorUtil::DeltaR(tau.p4, tau.p4genVis) < 0.3),
		lep_identified(lepton.hasRec && lepton.in_acceptance(lepton.p4) && lepton.id() && lepton.quality),
		tau_identified(tau.hasRec && tau.in_acceptance(tau.p4) && tau.decayModeFinding && tau.quality && tau.againstMuonTight && tau.againstElectronLoose && tau.looseIsolation3Hits),
		lep_isolated(lepton.hasRec && lepton.isolated()),
		selected(lep_identified && tau_identified && lep_isolated && ROOT::Math::VectorUtil::DeltaR(lepton.p4, tau.p4) > 0.3 && calcMt(lepton.p4.Px(), lepton.p4.Py(), vars.mvaPfMet * cos(vars.mvaPfMetPhi), vars.mvaPfMet * sin(vars.mvaPfMetPhi)) < 30.)
	{
	}

	const bool origAccepted;
	const bool genInitialAccepted;
	const bool genInitialMasscut;
	const bool genFinalAccepted;
	const bool genFinalMasscut;
	const bool genVisAccepted;
	const bool lep_matched;
	const bool tau_matched;
	const bool lep_identified;
	const bool tau_identified;
	const bool lep_isolated;
	const bool selected;
};

template<typename TLepton>
struct RadiationHistograms
{
	RadiationHistograms(TDirectory* dir, const char* name) {
		TDirectory* newdir = dir->mkdir(name);
		newdir->cd();

		const double relBinning[] = { 0., 0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.004, 0.006, 0.01, 0.015, 0.02, 0.03, 0.05, 0.1, 1.0 };
		const double absBinning[] = { 0., 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 100.0 };

		const int relNBins = sizeof(relBinning)/sizeof(relBinning[0]) - 1;
		const int absNBins = sizeof(absBinning)/sizeof(absBinning[0]) - 1;

		// TODO: These might need variable binning...!
		radAbs = new TH1D("radAbs", "Radiation;E_{loss}^{rad} [GeV];Entries", absNBins, absBinning);
		radRel = new TH1D("radRel", "Relative radiation;E_{loss}^{rad} / E_{total}", relNBins, relBinning);

		dir->cd();
	}

	void fill(const TreeVars& vars, const TLepton& lepton, const Tau& tau, float weight)
	{
		assert(lepton.hasOrigGen && tau.hasOrigGen);

		const double Einitial = lepton.p4origGenInitial.E() + tau.p4origGenInitial.E();
		const double Efinal = lepton.p4origGenFinal.E() + tau.p4origGenFinal.E();

		const double Eloss = Einitial - Efinal;
		assert(Eloss >= -1e-3);
		assert(Eloss <= Einitial+1e-3);

		const double ElossRel = Eloss / Einitial;

		radAbs->Fill(Eloss, weight);
		radRel->Fill(ElossRel, weight);
	}

	TH1D* radAbs;
	TH1D* radRel;

	TProfile* genMuonPt;
	TProfile* genDiMuonMass;

	TProfile* recMuonPt;
	TProfile* recDiMuonMass;

	TProfile* radiationFlag;
};

template<typename TLepton>
struct OrigHistograms
{
	OrigHistograms(TDirectory* dir, const char* name) {
		TDirectory* newdir = dir->mkdir(name);
		newdir->cd();

		origMuon1PtPlus = new TH1D("origMuon1PtPlus", "p_{T};p_{T} [GeV];Entries", 100, 0., 100.);
		origMuon1PtMinus = new TH1D("origMuon1PtMinus", "p_{T};p_{T} [GeV];Entries", 100, 0., 100.);
		origMuon1EtaPlus = new TH1D("origMuon1EtaPlus", "#eta;#eta;Entries", 50, -2.5, +2.5);
		origMuon1EtaMinus = new TH1D("origMuon1EtaMinus", "#eta;#eta;Entries", 50, -2.5, +2.5);
		origMuon2PtPlus = new TH1D("origMuon2PtPlus", "p_{T};p_{T} [GeV];Entries", 100, 0., 100.);
		origMuon2PtMinus = new TH1D("origMuon2PtMinus", "p_{T};p_{T} [GeV];Entries", 100, 0., 100.);
		origMuon2EtaPlus = new TH1D("origMuon2EtaPlus", "#eta;#eta;Entries", 50, -2.5, +2.5);
		origMuon2EtaMinus = new TH1D("origMuon2EtaMinus", "#eta;#eta;Entries", 50, -2.5, +2.5);

		embedMuon1PtPlus = new TH1D("embedMuon1PtPlus", "p_{T};p_{T} [GeV];Entries", 100, 0., 100.);
		embedMuon1PtMinus = new TH1D("embedMuon1PtMinus", "p_{T};p_{T} [GeV];Entries", 100, 0., 100.);
		embedMuon1EtaPlus = new TH1D("embedMuon1EtaPlus", "#eta;#eta;Entries", 50, -2.5, +2.5);
		embedMuon1EtaMinus = new TH1D("embedMuon1EtaMinus", "#eta;#eta;Entries", 50, -2.5, +2.5);
		embedMuon2PtPlus = new TH1D("embedMuon2PtPlus", "p_{T};p_{T} [GeV];Entries", 100, 0., 100.);
		embedMuon2PtMinus = new TH1D("embedMuon2PtMinus", "p_{T};p_{T} [GeV];Entries", 100, 0., 100.);
		embedMuon2EtaPlus = new TH1D("embedMuon2EtaPlus", "#eta;#eta;Entries", 50, -2.5, +2.5);
		embedMuon2EtaMinus = new TH1D("embedMuon2EtaMinus", "#eta;#eta;Entries", 50, -2.5, +2.5);

		deltaLepPt = new TH1D("deltaLepPt", "#Delta#p_{T};#p_{T}^{after} - p_{T}^{before} [GeV];Entries", 200, -100., +100.);
		deltaLepEta = new TH1D("deltaLepEta", "#Delta#eta;#eta_{after} - #eta_{before};Entries", 100, -5.0, +5.0);
		deltaLepPhi = new TH1D("deltaLepPhi", "#Delta#phi;#phi_{after} - #phi_{before};Entries", 100, -M_PI, +M_PI);

		deltaLepEtaVsEta = new TProfile("deltaLepEtaVsEta", "#Delta#eta vs. #eta;#eta_{before};#eta_{after} - #eta_{before}", 50, -2.5, +2.5);
		deltaTauEtaVsEta = new TProfile("deltaTauEtaVsEta", "#Delta#eta vs. #eta;#eta_{before};#eta_{after} - #eta_{before}", 50, -2.5, +2.5);

		dir->cd();
	}

	void fill(const TreeVars& vars, const TLepton& lepton, const Tau& tau, float weight)
	{
		if(lepton.genCharge > 0)
		{
			origMuon1PtPlus->Fill(lepton.p4orig.Pt(), weight);
			origMuon1EtaPlus->Fill(lepton.p4orig.Eta(), weight);
			origMuon2PtMinus->Fill(tau.p4orig.Pt(), weight);
			origMuon2EtaMinus->Fill(tau.p4orig.Eta(), weight);
			embedMuon1PtPlus->Fill(lepton.p4genInitial.Pt(), weight);
			embedMuon1EtaPlus->Fill(lepton.p4genInitial.Eta(), weight);
			embedMuon2PtMinus->Fill(tau.p4genInitial.Pt(), weight);
			embedMuon2EtaMinus->Fill(tau.p4genInitial.Eta(), weight);
		}
		else
		{
			origMuon1PtMinus->Fill(lepton.p4orig.Pt(), weight);
			origMuon1EtaMinus->Fill(lepton.p4orig.Eta(), weight);
			origMuon2PtPlus->Fill(tau.p4orig.Pt(), weight);
			origMuon2EtaPlus->Fill(tau.p4orig.Eta(), weight);
			embedMuon1PtMinus->Fill(lepton.p4genInitial.Pt(), weight);
			embedMuon1EtaMinus->Fill(lepton.p4genInitial.Eta(), weight);
			embedMuon2PtPlus->Fill(tau.p4genInitial.Pt(), weight);
			embedMuon2EtaPlus->Fill(tau.p4genInitial.Eta(), weight);
		}

		deltaLepPt->Fill(lepton.p4genInitial.Pt() - lepton.p4orig.Pt(), weight);
		deltaLepEta->Fill(lepton.p4genInitial.Eta() - lepton.p4orig.Eta(), weight);
		deltaLepPhi->Fill(ROOT::Math::VectorUtil::DeltaPhi(lepton.p4genInitial, lepton.p4orig), weight);

		deltaLepEtaVsEta->Fill(lepton.p4orig.Eta(), lepton.p4genInitial.Eta() - lepton.p4orig.Eta(), weight);
		deltaTauEtaVsEta->Fill(tau.p4orig.Eta(), tau.p4genInitial.Eta() - tau.p4orig.Eta(), weight);
	}

	TH1D* origMuon1PtPlus;
	TH1D* origMuon1PtMinus;
	TH1D* origMuon1EtaPlus;
	TH1D* origMuon1EtaMinus;
	TH1D* origMuon2PtPlus;
	TH1D* origMuon2PtMinus;
	TH1D* origMuon2EtaPlus;
	TH1D* origMuon2EtaMinus;

	TH1D* embedMuon1PtPlus;
	TH1D* embedMuon1PtMinus;
	TH1D* embedMuon1EtaPlus;
	TH1D* embedMuon1EtaMinus;
	TH1D* embedMuon2PtPlus;
	TH1D* embedMuon2PtMinus;
	TH1D* embedMuon2EtaPlus;
	TH1D* embedMuon2EtaMinus;

	TH1D* deltaLepPt;
	TH1D* deltaLepEta;
	TH1D* deltaLepPhi;
	TProfile* deltaLepEtaVsEta;
	TProfile* deltaTauEtaVsEta;
};

struct FillInfo
{
	template<typename TLepton>
	FillInfo(const TreeVars& vars, const TLepton& lep, const Tau& tau):
		lepHasGen(lep.hasGen),
		lepHasRec(lep.hasRec),
		tauHasGen(tau.hasGen),
		tauHasRec(tau.hasRec)
	{
		TLorentzVector met;
		met.SetPtEtaPhiE(vars.mvaPfMet, 0., vars.mvaPfMetPhi, vars.mvaPfMet);

		if(tauHasGen && lepHasGen)
		{
			const TLorentzVector diTauGenInitial = (tau.p4genInitial + lep.p4genInitial);
			const TLorentzVector diTauGenFinal = (tau.p4genFinal + lep.p4genFinal);
			const TLorentzVector diTauGenVis = (tau.p4genVis + lep.p4genVis);

			diTauGenInitialMass = diTauGenInitial.M();
			diTauGenFinalMass = diTauGenFinal.M();
			diTauGenVisMass = diTauGenVis.M();
			diTauGenPt = diTauGenFinal.Pt();
			diTauGenEta = diTauGenFinal.Eta();
		}

		if(lepHasGen)
		{
			lepGenPt = lep.p4genVis.Pt();
			lepGenEta = lep.p4genVis.Eta();
			lepGenPhi = lep.p4genVis.Phi();
		}

		if(lepHasRec)
		{
			lepPt = lep.p4.Pt();
			lepEta = lep.p4.Eta();
			lepPhi = lep.p4.Phi();
			lepIso = lep.isolation;

			lepChargedParticleIso = lep.chargedParticleIsolation;
			lepChargedHadronIso = lep.chargedHadronIsolation;
			lepNeutralHadronIso = lep.neutralHadronIsolation;
			lepPhotonIso = lep.photonIsolation;
			lepPUIso = lep.puIsolation;

			if(lepHasGen)
				lepMomentumResolution = (lep.p4.Pt() - lep.p4genVis.Pt()) / lep.p4genVis.Pt();
		}

		if(tauHasGen)
		{
			tauGenPt = tau.p4genVis.Pt();
			tauGenEta = tau.p4genVis.Eta();
			tauGenPhi = tau.p4genVis.Phi();
		}

		if(tauHasRec)
		{
			tauPt = tau.p4.Pt();
			tauEta = tau.p4.Eta();
			tauPhi = tau.p4.Phi();
			tauIso = tau.isolation3Hits;

			if(tauHasGen)
				tauMomentumResolution = (tau.p4.Pt() - tau.p4genVis.Pt()) / tau.p4genVis.Pt();

			const double x = lep.p4.Px() * tau.p4.Py() - lep.p4.Py() * tau.p4.Px();
			x1 = x / (x + tau.p4.Py() * met.Px() - tau.p4.Px() * met.Py());
			x2 = x / (x + lep.p4.Px() * met.Py() - lep.p4.Py() * met.Px());
		}

		if(tauHasRec && lepHasRec)
		{
			visibleMass = (lep.p4 + tau.p4).M();
			collinearMass = -1.0f;
			if(x1 >= 0. && x2 >= 0. && x1 <= 1. && x2 <= 1.) collinearMass = (lep.p4 * (1.0/x1) + tau.p4 * (1.0/x2)).M();
			svfitMass = vars.svfitMass;
			higgsPt = (lep.p4 + tau.p4 + met).Pt();
		}

		caloMet = vars.caloMet;
		pfMet = vars.pfMet;
		mvaPfMet = vars.mvaPfMet;
		nPV = vars.nPV;
		nJets30 = 0;
		nBJets20 = 0;
		leadJetPt = -1.0f;
		for(unsigned int i = 0; i < vars.nJets; ++i)
		{
			if(nJets30 == 0 || vars.jetPt[i] > leadJetPt)
				leadJetPt = vars.jetPt[i];
			if(vars.jetPt[i] > 30.)
				++nJets30;
			if(vars.jetPt[i] > 20. && vars.jetBTag[i] > 0.679)
				++nBJets20;
		}
	}

	const bool lepHasGen;
	const bool lepHasRec;
	const bool tauHasGen;
	const bool tauHasRec;

	// ditau gen
	float diTauGenInitialMass;
	float diTauGenFinalMass;
	float diTauGenVisMass;
	float diTauGenPt;
	float diTauGenEta;

	// lepton
	float lepGenPt;
	float lepGenEta;
	float lepGenPhi;
	float lepPt;
	float lepEta;
	float lepPhi;
	float lepIso;

	float lepChargedParticleIso;
	float lepChargedHadronIso;
	float lepNeutralHadronIso;
	float lepPhotonIso;
	float lepPUIso;

	float lepMomentumResolution;

	// tau
	float tauGenPt;
	float tauGenEta;
	float tauGenPhi;
	float tauPt;
	float tauEta;
	float tauPhi;
	float tauIso;

	float tauMomentumResolution;

	// di-tau
	float x1;
	float x2;

	float visibleMass;
	float collinearMass;
	float svfitMass;
	float higgsPt;

	// event
	float caloMet;
	float pfMet;
	float mvaPfMet;
	unsigned int nPV;
	unsigned int nJets30;
	unsigned int nBJets20;
	float leadJetPt;
};

template<typename TLepton>
struct LeptonHistograms {};

template<>
struct LeptonHistograms<Electron>
{
	LeptonHistograms() {
		electronGenPt = new TH1D("electronGenPt", "Generated electron transverse momentum;p_{T};Entries/0.1 GeV", 200, 0.0, 200.0);
		electronGenEta = new TH1D("electronGenEta", "Generated electron pseudo-rapidity;#eta;Entries/0.1", 180, -9.0, 9.0);
		electronGenPhi = new TH1D("electronGenPhi", "Generated electron azimuthal angle;#phi;Entries", 100, -M_PI, M_PI);
		electronPt = new TH1D("electronPt", "Electron transverse momentum;p_{T};Entries/0.1 GeV", 200, 0.0, 200.0);
		electronEta = new TH1D("electronEta", "Electron pseudo-rapidity;#eta;Entries/0.1", 180, -9.0, 9.0);
		electronPhi = new TH1D("electronPhi", "Electron azimuthal angle;#phi;Entries", 100, -M_PI, M_PI);
		electronIso = new TH1D("electronIso", "Electron relative PF isolation;I_{rel}^{PF}", 100, 0.0, 1.0);
		electronChargedParticleIso = new TH1D("electronChargedParticleIso", "Electron relative charged particle isolation;I_{rel}^{PF}", 100, 0.0, 1.0);
		electronChargedHadronIso = new TH1D("electronChargedHadronIso", "Electron relative charged hadron isolation;I_{rel}^{PF}", 100, 0.0, 1.0);
		electronNeutralHadronIso = new TH1D("electronNeutralHadronIso", "Electron relative neutral hadron isolation;I_{rel}^{PF}", 100, 0.0, 1.0);
		electronPhotonIso = new TH1D("electronPhotonIso", "Electron relative photon isolation;I_{rel}^{PF}", 100, 0.0, 1.0);
		electronPUIso = new TH1D("electronPUIso", "Electron relative PU isolation;I_{rel}^{PF}", 100, 0.0, 1.0);
		electronMomentumResolution = new TH1D("electronMomentumResolution", "Electron p_{T} resolution;(p_{T}^{rec} - p_{T}^{gen}) / p_{T}^{gen};Entries", 400, -0.2, 0.2);
		electronMomentumResolutionVsPt = new TProfile("electronMomentumResolutionVsPt", "Electron p_{T} resolution;p_{T}^{gen};(p_{T}^{rec} - p_{T}^{gen}) / p_{T}^{gen}", 100, 0.0, 100.);
	}

	void fill(const FillInfo& info, float weight)
	{
		if(info.lepHasGen)
		{
			electronGenPt->Fill(info.lepGenPt, weight);
			electronGenEta->Fill(info.lepGenEta, weight);
			electronGenPhi->Fill(info.lepGenPhi, weight);
		}

		if(info.lepHasRec)
		{
			electronPt->Fill(info.lepPt, weight);
			electronEta->Fill(info.lepEta, weight);
			electronPhi->Fill(info.lepPhi, weight);
			electronIso->Fill(info.lepIso, weight);
			electronChargedParticleIso->Fill(info.lepChargedParticleIso, weight);
			electronChargedHadronIso->Fill(info.lepChargedHadronIso, weight);
			electronNeutralHadronIso->Fill(info.lepNeutralHadronIso, weight);
			electronPhotonIso->Fill(info.lepPhotonIso, weight);
			electronPUIso->Fill(info.lepPUIso, weight);

			if(info.lepHasGen)
			{
				electronMomentumResolution->Fill(info.lepMomentumResolution, weight);
				electronMomentumResolutionVsPt->Fill(info.lepGenPt, info.lepMomentumResolution, weight);
			}
		}
	}

	TH1D* electronGenPt;
	TH1D* electronGenEta;
	TH1D* electronGenPhi;
	TH1D* electronPt;
	TH1D* electronEta;
	TH1D* electronPhi;
	TH1D* electronIso;

	TH1D* electronChargedParticleIso;
	TH1D* electronChargedHadronIso;
	TH1D* electronNeutralHadronIso;
	TH1D* electronPhotonIso;
	TH1D* electronPUIso;

	TH1D* electronMomentumResolution;
	TProfile* electronMomentumResolutionVsPt;
};

template<>
struct LeptonHistograms<Muon>
{
	LeptonHistograms() {
		muonGenPt = new TH1D("muonGenPt", "Generated muon transverse momentum;p_{T};Entries/0.1 GeV", 200, 0.0, 200.0);
		muonGenEta = new TH1D("muonGenEta", "Generated muon pseudo-rapidity;#eta;Entries/0.1", 180, -9.0, 9.0);
		muonGenPhi = new TH1D("muonGenPhi", "Generated muon azimuthal angle;#phi;Entries", 100, -M_PI, M_PI);
		muonPt = new TH1D("muonPt", "Muon transverse momentum;p_{T};Entries/0.1 GeV", 200, 0.0, 200.0);
		muonEta = new TH1D("muonEta", "Muon pseudo-rapidity;#eta;Entries/0.1", 180, -9.0, 9.0);
		muonPhi = new TH1D("muonPhi", "Muon azimuthal angle;#phi;Entries", 100, -M_PI, M_PI);
		muonIso = new TH1D("muonIso", "Muon relative PF isolation;I_{rel}^{PF}", 100, 0.0, 1.0);
		muonChargedParticleIso = new TH1D("muonChargedParticleIso", "Muon relative charged particle isolation;I_{rel}^{PF}", 100, 0.0, 1.0);
		muonChargedHadronIso = new TH1D("muonChargedHadronIso", "Muon relative charged hadron isolation;I_{rel}^{PF}", 100, 0.0, 1.0);
		muonNeutralHadronIso = new TH1D("muonNeutralHadronIso", "Muon relative neutral hadron isolation;I_{rel}^{PF}", 100, 0.0, 1.0);
		muonPhotonIso = new TH1D("muonPhotonIso", "Muon relative photon isolation;I_{rel}^{PF}", 100, 0.0, 1.0);
		muonPUIso = new TH1D("muonPUIso", "Muon relative PU isolation;I_{rel}^{PF}", 100, 0.0, 1.0);
		muonMomentumResolution = new TH1D("muonMomentumResolution", "Muon p_{T} resolution;(p_{T}^{rec} - p_{T}^{gen}) / p_{T}^{gen};Entries", 200, -0.1, 0.1);
		muonMomentumResolutionVsPt = new TProfile("muonMomentumResolutionVsPt", "Muon p_{T} resolution;p_{T}^{gen};(p_{T}^{rec} - p_{T}^{gen}) / p_{T}^{gen}", 100, 0.0, 100.);
	}

	void fill(const FillInfo& info, float weight)
	{
		if(info.lepHasGen)
		{
			muonGenPt->Fill(info.lepGenPt, weight);
			muonGenEta->Fill(info.lepGenEta, weight);
			muonGenPhi->Fill(info.lepGenPhi, weight);
		}

		if(info.lepHasRec)
		{
			muonPt->Fill(info.lepPt, weight);
			muonEta->Fill(info.lepEta, weight);
			muonPhi->Fill(info.lepPhi, weight);
			muonIso->Fill(info.lepIso, weight);
			muonChargedParticleIso->Fill(info.lepChargedParticleIso, weight);
			muonChargedHadronIso->Fill(info.lepChargedHadronIso, weight);
			muonNeutralHadronIso->Fill(info.lepNeutralHadronIso, weight);
			muonPhotonIso->Fill(info.lepPhotonIso, weight);
			muonPUIso->Fill(info.lepPUIso, weight);

			if(info.lepHasGen)
			{
				muonMomentumResolution->Fill(info.lepMomentumResolution, weight);
				muonMomentumResolutionVsPt->Fill(info.lepGenPt, info.lepMomentumResolution, weight);
			}
		}
	}

	TH1D* muonGenPt;
	TH1D* muonGenEta;
	TH1D* muonGenPhi;
	TH1D* muonPt;
	TH1D* muonEta;
	TH1D* muonPhi;
	TH1D* muonIso;

	TH1D* muonChargedParticleIso;
	TH1D* muonChargedHadronIso;
	TH1D* muonNeutralHadronIso;
	TH1D* muonPhotonIso;
	TH1D* muonPUIso;

	TH1D* muonMomentumResolution;
	TProfile* muonMomentumResolutionVsPt;
};

template<typename TLepton>
struct Histograms
{
	Histograms(TDirectory* dir, const char* name) {
		TDirectory* newdir = dir->mkdir(name);
		newdir->cd();

		nPV = new TH1D("nPV", "Number of primary vertices;N_{PV};Entries", 100, -0.5, 99.5);
		nPV_u = new TH1D("nPV_u", "Number of primary vertices, unweighted;N_{PV};Entries", 100, -0.5, 99.5);

		leptonHistograms.reset(new LeptonHistograms<TLepton>());

		diTauGenInitialMass = new TH1D("diTauGenInitialMass", "Generated di-tau mass;M_{#tau#tau};Entries", 200, 0., 200.);
		diTauGenFinalMass = new TH1D("diTauGenFinalMass", "Generated di-tau mass;M_{#tau#tau};Entries", 200, 0., 200.);
		diTauGenVisMass = new TH1D("diTauGenVisMass", "Generated di-tau visible mass;M_{#tau#tau}^{vis};Entries", 200, 0., 200.);
		diTauGenPt = new TH1D("diTauGenPt", "Generated di-tau transverse momentum;p_{T};Entries", 150, 0., 150.);
		diTauGenEta = new TH1D("diTauGenEta", "Generated di-tau pseudorapidity;#eta;Entries/0.1", 180, -9.0, 9.0);

		tauGenPt = new TH1D("tauGenPt", "Generated tau transverse momentum;p_{T};Entries/0.1 GeV", 200, 0.0, 200.0);
		tauGenEta = new TH1D("tauGenEta", "Generated tau pseudo-rapidity;#eta;Entries/0.1", 180, -9.0, 9.0);
		tauGenPhi = new TH1D("tauGenPhi", "Generated tau azimuthal angle;#phi;Entries", 100, -M_PI, M_PI);
		tauPt = new TH1D("tauPt", "Tau transverse momentum;p_{T};Entries/0.1 GeV", 200, 0.0, 200.0);
		tauEta = new TH1D("tauEta", "Tau pseudo-rapidity;#eta;Entries/0.1", 180, -9.0, 9.0);
		tauPhi = new TH1D("tauPhi", "Tau azimuthal angle;#phi;Entries", 100, -M_PI, M_PI);
		tauIso = new TH1D("tauIso", "Tau absolute PF isolation (3hits);I_{3hits}^{PF}", 100, 0.0, 10.0);
		tauMomentumResolution = new TH1D("tauMomentumResolution", "Tau p_{T} resolution;(p_{T}^{rec} - p_{T}^{gen}) / p_{T}^{gen};Entries", 900, -0.6, 0.3);
		tauMomentumResolutionVsPt = new TProfile("tauMomentumResolutionVsPt", "Tau p_{T} resolution;p_{T}^{gen};(p_{T}^{rec} - p_{T}^{gen}) / p_{T}^{gen}", 100, 0.0, 100.);

		x1x2 = new TH2D("x1x2", "x_{1} vs. x_{2}", 150, -1., +2., 150, -1., +2.);

		visibleMass = new TH1D("visibleMass", "Visible mass;M_{#tau#tau}^{vis};Entries/GeV", 400, 0., 400.);
		collinearMass = new TH1D("collinearMass", "Collinear mass;M_{#tau#tau}^{col};Entries/GeV", 400, 0., 400.);
		svfitMass = new TH1D("svfitMass", "SVfit mass;M_{#tau#tau}^{svfit};Entries/GeV", 400, 0., 400.);
		svfitMassHigh = new TH1D("svfitMassHigh", "SVfit mass;M_{#tau#tau}^{svfit};Entries/GeV", 400, 0., 400.);
		higgsPt = new TH1D("higgsPt", "Higgs p_{T};p_{T}^{H};Entries/GeV", 200, 0., 200.);
		caloMet = new TH1D("caloMet", "Calorimeter-based MET;E_{T}^{miss};Entries/GeV", 200, 0., 200.);
		pfMet = new TH1D("pfMet", "PFlow-based MET;E_{T}^{miss};Entries/GeV", 200, 0., 200.);
		mvaPfMet = new TH1D("mvaPfMet", "MVA based MET;E_{T}^{miss};Entries/GeV", 200, 0., 200.);
		nJets30 = new TH1D("nJets30", "Jet multiplicity;n_{Jets};Entries", 20, -0.5, 19.5);
		nBJets20 = new TH1D("nBJets20", "B-Jet multiplicity;n_{Jets};Entries", 20, -0.5, 19.5);
		leadJetPt = new TH1D("leadJetPt", "Leading jet transverse momentum;p_{T};Entries", 150, 0., 150.);

		dir->cd();
	}

	void fill(const FillInfo& info, float weight)
	{
		nPV->Fill(info.nPV, weight);
		nPV_u->Fill(info.nPV);

		leptonHistograms->fill(info, weight);

		if(info.lepHasGen && info.tauHasGen)
		{
			diTauGenInitialMass->Fill(info.diTauGenInitialMass, weight);
			diTauGenFinalMass->Fill(info.diTauGenFinalMass, weight);
			diTauGenVisMass->Fill(info.diTauGenVisMass, weight);
			diTauGenPt->Fill(info.diTauGenPt, weight);
			diTauGenEta->Fill(info.diTauGenEta, weight);
		}

		if(info.tauHasGen)
		{
			tauGenPt->Fill(info.tauGenPt, weight);
			tauGenEta->Fill(info.tauGenEta, weight);
			tauGenPhi->Fill(info.tauGenPhi, weight);
		}

		if(info.tauHasRec)
		{
			tauPt->Fill(info.tauPt, weight);
			tauEta->Fill(info.tauEta, weight);
			tauPhi->Fill(info.tauPhi, weight);
			tauIso->Fill(info.tauIso, weight);

			if(info.tauHasGen)
			{
				tauMomentumResolution->Fill(info.tauMomentumResolution, weight);
				tauMomentumResolutionVsPt->Fill(info.tauGenPt, info.tauMomentumResolution, weight);
			}
		}

		if(info.lepHasRec && info.tauHasRec)
		{
			x1x2->Fill(info.x1, info.x2, weight);
			visibleMass->Fill(info.visibleMass, weight);
			if(info.collinearMass > 0.0f) collinearMass->Fill(info.collinearMass, weight);
			if(info.svfitMass > 0.0f) svfitMass->Fill(info.svfitMass, weight);
			if(info.svfitMass > 0.0f && info.lepPt > 22. && info.tauPt > 30.) svfitMassHigh->Fill(info.svfitMass, weight);
			higgsPt->Fill(info.higgsPt, weight);
		}

		caloMet->Fill(info.caloMet, weight);
		pfMet->Fill(info.pfMet, weight);
		mvaPfMet->Fill(info.mvaPfMet, weight);
		nJets30->Fill(info.nJets30, weight);
		nBJets20->Fill(info.nBJets20, weight);
		if(info.leadJetPt > 0.0f) leadJetPt->Fill(info.leadJetPt, weight);
	}

	TH1D* nPV;
	TH1D* nPV_u;

	std::auto_ptr<LeptonHistograms<TLepton> > leptonHistograms;

	TH1D* diTauGenInitialMass;
	TH1D* diTauGenFinalMass;
	TH1D* diTauGenVisMass;
	TH1D* diTauGenPt;
	TH1D* diTauGenEta;

	TH1D* tauGenPt;
	TH1D* tauGenEta;
	TH1D* tauGenPhi;
	TH1D* tauPt;
	TH1D* tauEta;
	TH1D* tauPhi;
	TH1D* tauIso;

	TH1D* tauMomentumResolution;
	TProfile* tauMomentumResolutionVsPt;

	TH2D* x1x2;

	TH1D* visibleMass;
	TH1D* collinearMass;
	TH1D* svfitMass;
	TH1D* svfitMassHigh;
	TH1D* higgsPt;
	TH1D* caloMet;
	TH1D* pfMet;
	TH1D* mvaPfMet;
	TH1D* nJets30;
	TH1D* nBJets20;
	TH1D* leadJetPt;
};

template<typename TLepton>
struct Category
{
	Category(TDirectory* dir, const char* name, bool is_embedded):
		isEmbedded(is_embedded)
	{
		TDirectory* newdir = dir->mkdir(name);
		newdir->cd();

		if(is_embedded) orig.reset(new OrigHistograms<TLepton>(newdir, "orig"));
		genInitialMasscut.reset(new Histograms<TLepton>(newdir, "genInitialMasscut"));
		genFinalMasscut.reset(new Histograms<TLepton>(newdir, "genFinalMasscut"));
		genVisMasscut.reset(new Histograms<TLepton>(newdir, "genVisMasscut"));
		selected.reset(new Histograms<TLepton>(newdir, "selected"));

		lep_identified.reset(new Histograms<TLepton>(newdir, "lep_identified"));
		lep_isolated.reset(new Histograms<TLepton>(newdir, "lep_isolated"));
		tau_identified.reset(new Histograms<TLepton>(newdir, "tau_identified"));

		if(is_embedded) origRadiation.reset(new RadiationHistograms<TLepton>(newdir, "origRadiation"));
		if(is_embedded) recSelRadiation.reset(new RadiationHistograms<TLepton>(newdir, "recSelRadiation"));
		dir->cd();
	}

	void fill(const EventSelector<TLepton>& selector, const FillInfo& info, const TreeVars& vars, const TLepton& lepton, const Tau& tau, float weight)
	{
		if(selector.origAccepted && orig.get()) orig->fill(vars, lepton, tau, weight); // TODO: Add these variables to fillinfo

		// Require embedded events to be within acceptance for MC/Embedded comparison
		const bool origAccepted = (!isEmbedded || selector.origAccepted);
		if(origAccepted && selector.genInitialAccepted && selector.genInitialMasscut)
			genInitialMasscut->fill(info, weight);
		if(origAccepted && selector.genFinalAccepted && selector.genFinalMasscut)
			genFinalMasscut->fill(info, weight);

		if(origAccepted && selector.genFinalAccepted && selector.genFinalMasscut && selector.genVisAccepted)
		{
			genVisMasscut->fill(info, weight);
			if(selector.lep_matched && selector.lep_identified)
				lep_identified->fill(info, weight);
			if(selector.lep_matched && selector.lep_identified && selector.lep_isolated)
				lep_isolated->fill(info, weight);
			if(selector.tau_matched && selector.tau_identified)
				tau_identified->fill(info, weight);
			if(selector.lep_matched && selector.tau_matched && selector.selected)
				selected->fill(info, weight);
		}

		if(selector.origAccepted && origRadiation.get()) origRadiation->fill(vars, lepton, tau, weight); // TODO: add these variables to fillinfo
		if(selector.origAccepted && selector.genFinalAccepted && selector.genFinalMasscut && selector.genVisAccepted && selector.lep_matched && selector.tau_matched && selector.selected && recSelRadiation.get())
			recSelRadiation->fill(vars, lepton, tau, weight); // TODO: add these variables to fillinfo
	}

	const bool isEmbedded;

	std::auto_ptr<OrigHistograms<TLepton> > orig;
	std::auto_ptr<Histograms<TLepton> > genInitialMasscut;
	std::auto_ptr<Histograms<TLepton> > genFinalMasscut;
	std::auto_ptr<Histograms<TLepton> > genVisMasscut;
	std::auto_ptr<Histograms<TLepton> > selected;

	std::auto_ptr<Histograms<TLepton> > lep_identified;
	std::auto_ptr<Histograms<TLepton> > lep_isolated;
	std::auto_ptr<Histograms<TLepton> > tau_identified;

	std::auto_ptr<RadiationHistograms<TLepton> > origRadiation;
	std::auto_ptr<RadiationHistograms<TLepton> > recSelRadiation;
};

template<typename TLepton>
struct Out
{
	Out(TDirectory* dir, bool is_data, bool is_embedded):
		isEmbedded(is_embedded),
		inclusive(dir, "inclusive", is_embedded)
	{
		if(is_embedded)
		{
			inclusiveNoRadiation.reset(new Category<TLepton>(dir, "inclusiveNoRadiation", is_embedded));
			inclusiveRadiationFilter.reset(new Category<TLepton>(dir, "inclusiveRadiationFilter", is_embedded));
			inclusiveRadiationFilter2Sel1.reset(new Category<TLepton>(dir, "inclusiveRadiationFilter2Sel1", is_embedded));
			inclusiveRadiationFilter2Sel2.reset(new Category<TLepton>(dir, "inclusiveRadiationFilter2Sel2", is_embedded));
			inclusiveRadiationFilter2Sel3.reset(new Category<TLepton>(dir, "inclusiveRadiationFilter2Sel3", is_embedded));
		}

		nEvents = new TH1D("nEvents", "Total number of weighted events", 1, 0, 1);
		nPV = new TH1D("nPV", "Number of Primary Vertices", 100, -0.5, 99.5);
		if(!is_data) nTrueInteractions = new TH1D("nTrueInteractions", "Number of True Interactions", 100, -0.5, 99.5);
	}

	void fill(const TreeVars& vars, const TLepton& lepton, const Tau& tau, float weight)
	{
		const EventSelector<TLepton> selector(vars, lepton, tau);
		const FillInfo info(vars, lepton, tau);

		inclusive.fill(selector, info, vars, lepton, tau, weight);

		if(inclusiveNoRadiation.get())
		{
			assert(lepton.hasOrigGen && tau.hasOrigGen);
			if(lorentzEqual(lepton.p4origGenInitial, lepton.p4origGenFinal) && lorentzEqual(tau.p4origGenInitial, tau.p4origGenFinal))
				inclusiveNoRadiation->fill(selector, info ,vars, lepton, tau, weight);
		}

		if(inclusiveRadiationFilter.get())
		{
			if(vars.muonRadiationFilter)
				inclusiveRadiationFilter->fill(selector, info, vars, lepton, tau, weight);
		}

		if(inclusiveRadiationFilter2Sel1.get())
		{
			if(vars.muonRadiationFilter2Sel1)
				inclusiveRadiationFilter2Sel1->fill(selector, info, vars, lepton, tau, weight);
		}

		if(inclusiveRadiationFilter2Sel2.get())
		{
			if(vars.muonRadiationFilter2Sel2)
				inclusiveRadiationFilter2Sel2->fill(selector, info, vars, lepton, tau, weight);
		}

		if(inclusiveRadiationFilter2Sel3.get())
		{
			if(vars.muonRadiationFilter2Sel3)
				inclusiveRadiationFilter2Sel3->fill(selector, info, vars, lepton, tau, weight);
		}
	}

	void add_events(float n)
	{
		nEvents->SetBinContent(1, static_cast<double>(nEvents->GetBinContent(1) + n));
	}

	void add_pv(TH1F* nPV)
	{
		this->nPV->Add(nPV);
	}

	void add_interactions(TH1F* nTrueInteractions)
	{
		this->nTrueInteractions->Add(nTrueInteractions);
	}

	const bool isEmbedded;

	Category<TLepton> inclusive;
	std::auto_ptr<Category<TLepton> > inclusiveNoRadiation;
	std::auto_ptr<Category<TLepton> > inclusiveRadiationFilter;
	std::auto_ptr<Category<TLepton> > inclusiveRadiationFilter2Sel1;
	std::auto_ptr<Category<TLepton> > inclusiveRadiationFilter2Sel2;
	std::auto_ptr<Category<TLepton> > inclusiveRadiationFilter2Sel3;
	TH1D* nEvents;
	TH1D* nPV;
	TH1D* nTrueInteractions;
};

template<typename TLepton>
void register_lepton_branches(TTree* tree, TreeVars& vars);

template<>
void register_lepton_branches<Electron>(TTree* tree, TreeVars& vars)
{
	tree->SetBranchAddress("ElectronPt", &vars.electronPt);
	tree->SetBranchAddress("ElectronEta", &vars.electronEta);
	tree->SetBranchAddress("ElectronPhi", &vars.electronPhi);
	tree->SetBranchAddress("ElectronE", &vars.electronE);
	tree->SetBranchAddress("ElectronCharge", &vars.electronCharge);
	tree->SetBranchAddress("ElectronChargedParticlePtPfIso04", &vars.electronChargedParticlePfIso04);
	tree->SetBranchAddress("ElectronChargedHadronPtPfIso04", &vars.electronChargedHadronPfIso04);
	tree->SetBranchAddress("ElectronNeutralHadronEtPfIso04", &vars.electronNeutralHadronPfIso04);
	tree->SetBranchAddress("ElectronPhotonEtPfIso04", &vars.electronPhotonPfIso04);
	tree->SetBranchAddress("ElectronPUPtPfIso04", &vars.electronPUPfIso04);
	tree->SetBranchAddress("ElectronMVA", &vars.electronMvaId);
	tree->SetBranchAddress("ElectronQuality", &vars.electronQuality);
	tree->SetBranchAddress("ElectronHasConversion", &vars.electronHasConversion);
	tree->SetBranchAddress("ElectronGenPt", &vars.electronGenPt);
	tree->SetBranchAddress("ElectronGenEta", &vars.electronGenEta);
	tree->SetBranchAddress("ElectronGenPhi", &vars.electronGenPhi);
	tree->SetBranchAddress("ElectronGenE", &vars.electronGenE);
	tree->SetBranchAddress("ElectronGenFinalPt", &vars.electronGenFinalPt);
	tree->SetBranchAddress("ElectronGenFinalEta", &vars.electronGenFinalEta);
	tree->SetBranchAddress("ElectronGenFinalPhi", &vars.electronGenFinalPhi);
	tree->SetBranchAddress("ElectronGenFinalE", &vars.electronGenFinalE);
	tree->SetBranchAddress("ElectronGenCharge", &vars.electronGenCharge);
	tree->SetBranchAddress("ElectronGenPtVis", &vars.electronGenPtVis);
	tree->SetBranchAddress("ElectronGenEtaVis", &vars.electronGenEtaVis);
	tree->SetBranchAddress("ElectronGenPhiVis", &vars.electronGenPhiVis);
	tree->SetBranchAddress("ElectronGenEVis", &vars.electronGenEVis);
}

template<>
void register_lepton_branches<Muon>(TTree* tree, TreeVars& vars)
{
	tree->SetBranchAddress("MuonPt", &vars.muonPt);
	tree->SetBranchAddress("MuonEta", &vars.muonEta);
	tree->SetBranchAddress("MuonPhi", &vars.muonPhi);
	tree->SetBranchAddress("MuonE", &vars.muonE);
	tree->SetBranchAddress("MuonCharge", &vars.muonCharge);
	tree->SetBranchAddress("MuonChargedParticlePtPfIso04", &vars.muonChargedParticlePfIso04);
	tree->SetBranchAddress("MuonChargedHadronPtPfIso04", &vars.muonChargedHadronPfIso04);
	tree->SetBranchAddress("MuonNeutralHadronEtPfIso04", &vars.muonNeutralHadronPfIso04);
	tree->SetBranchAddress("MuonPhotonEtPfIso04", &vars.muonPhotonPfIso04);
	tree->SetBranchAddress("MuonPUPtPfIso04", &vars.muonPUPfIso04);
	tree->SetBranchAddress("MuonQuality", &vars.muonQuality);
	tree->SetBranchAddress("MuonGenPt", &vars.muonGenPt);
	tree->SetBranchAddress("MuonGenEta", &vars.muonGenEta);
	tree->SetBranchAddress("MuonGenPhi", &vars.muonGenPhi);
	tree->SetBranchAddress("MuonGenE", &vars.muonGenE);
	tree->SetBranchAddress("MuonGenFinalPt", &vars.muonGenFinalPt);
	tree->SetBranchAddress("MuonGenFinalEta", &vars.muonGenFinalEta);
	tree->SetBranchAddress("MuonGenFinalPhi", &vars.muonGenFinalPhi);
	tree->SetBranchAddress("MuonGenFinalE", &vars.muonGenFinalE);
	tree->SetBranchAddress("MuonGenCharge", &vars.muonGenCharge);
	tree->SetBranchAddress("MuonGenPtVis", &vars.muonGenPtVis);
	tree->SetBranchAddress("MuonGenEtaVis", &vars.muonGenEtaVis);
	tree->SetBranchAddress("MuonGenPhiVis", &vars.muonGenPhiVis);
	tree->SetBranchAddress("MuonGenEVis", &vars.muonGenEVis);
}

template<typename TLepton>
void process_file(const char* filename, Out<TLepton>& out, bool is_data, bool is_embedded, bool is_rh_embedded, const ZmumuEvtSelEffCorrWeight* zmumuWeight)
{
	TreeVars vars;

	std::string directory = "TauAna";
	TFile* file = new TFile(filename, "READ");
	TTree* tree = dynamic_cast<TTree*>(file->Get((directory + "/MuonTree").c_str()));

	tree->SetBranchAddress("Run", &vars.run);
	tree->SetBranchAddress("Lumi", &vars.lumi);
	tree->SetBranchAddress("Event", &vars.event);

	tree->SetBranchAddress("MinVisPtFilterWeight", &vars.minVisPtFilterWeight);
	tree->SetBranchAddress("TauSpinWeight", &vars.tauSpinWeight);
	tree->SetBranchAddress("MuonRadiationFilter", &vars.muonRadiationFilter);
	tree->SetBranchAddress("MuonRadiationFilter2Sel1", &vars.muonRadiationFilter2Sel1);
	tree->SetBranchAddress("MuonRadiationFilter2Sel2", &vars.muonRadiationFilter2Sel2);
	tree->SetBranchAddress("MuonRadiationFilter2Sel3", &vars.muonRadiationFilter2Sel3);

	tree->SetBranchAddress("NPV", &vars.nPV);
	if(!is_data) tree->SetBranchAddress("NTrueInteractions", &vars.nTrueInteractions);

	if(is_embedded)
	{
		tree->SetBranchAddress("OrigHltMu17Mu8", &vars.origHltMu17Mu8);

		tree->SetBranchAddress("OrigPosMuonPt", &vars.origPosMuonPt);
		tree->SetBranchAddress("OrigPosMuonEta", &vars.origPosMuonEta);
		tree->SetBranchAddress("OrigPosMuonPhi", &vars.origPosMuonPhi);
		tree->SetBranchAddress("OrigPosMuonE", &vars.origPosMuonE);
		tree->SetBranchAddress("OrigGenPosChargedLepPt", &vars.origGenPosMuonPt);
		tree->SetBranchAddress("OrigGenPosChargedLepEta", &vars.origGenPosMuonEta);
		tree->SetBranchAddress("OrigGenPosChargedLepPhi", &vars.origGenPosMuonPhi);
		tree->SetBranchAddress("OrigGenPosChargedLepE", &vars.origGenPosMuonE);
		tree->SetBranchAddress("OrigGenPosChargedLepFinalPt", &vars.origGenPosMuonFinalPt);
		tree->SetBranchAddress("OrigGenPosChargedLepFinalEta", &vars.origGenPosMuonFinalEta);
		tree->SetBranchAddress("OrigGenPosChargedLepFinalPhi", &vars.origGenPosMuonFinalPhi);
		tree->SetBranchAddress("OrigGenPosChargedLepFinalE", &vars.origGenPosMuonFinalE);
		tree->SetBranchAddress("OrigGenPosChargedLepPDG", &vars.origGenPosMuonPDG);

		tree->SetBranchAddress("OrigNegMuonPt", &vars.origNegMuonPt);
		tree->SetBranchAddress("OrigNegMuonEta", &vars.origNegMuonEta);
		tree->SetBranchAddress("OrigNegMuonPhi", &vars.origNegMuonPhi);
		tree->SetBranchAddress("OrigNegMuonE", &vars.origNegMuonE);
		tree->SetBranchAddress("OrigGenNegChargedLepPt", &vars.origGenNegMuonPt);
		tree->SetBranchAddress("OrigGenNegChargedLepEta", &vars.origGenNegMuonEta);
		tree->SetBranchAddress("OrigGenNegChargedLepPhi", &vars.origGenNegMuonPhi);
		tree->SetBranchAddress("OrigGenNegChargedLepE", &vars.origGenNegMuonE);
		tree->SetBranchAddress("OrigGenNegChargedLepFinalPt", &vars.origGenNegMuonFinalPt);
		tree->SetBranchAddress("OrigGenNegChargedLepFinalEta", &vars.origGenNegMuonFinalEta);
		tree->SetBranchAddress("OrigGenNegChargedLepFinalPhi", &vars.origGenNegMuonFinalPhi);
		tree->SetBranchAddress("OrigGenNegChargedLepFinalE", &vars.origGenNegMuonFinalE);
		tree->SetBranchAddress("OrigGenNegChargedLepPDG", &vars.origGenNegMuonPDG);
	}

	tree->SetBranchAddress("CaloMet", &vars.caloMet);
	tree->SetBranchAddress("CaloMetPhi", &vars.caloMetPhi);
	tree->SetBranchAddress("PfMet", &vars.pfMet);
	tree->SetBranchAddress("PfMetPhi", &vars.pfMetPhi);
	tree->SetBranchAddress("MvaPfMet", &vars.mvaPfMet);
	tree->SetBranchAddress("MvaPfMetPhi", &vars.mvaPfMetPhi);

	register_lepton_branches<TLepton>(tree, vars);

	tree->SetBranchAddress("TauPt", &vars.tauPt);
	tree->SetBranchAddress("TauEta", &vars.tauEta);
	tree->SetBranchAddress("TauPhi", &vars.tauPhi);
	tree->SetBranchAddress("TauE", &vars.tauE);
	tree->SetBranchAddress("TauCharge", &vars.tauCharge);
	tree->SetBranchAddress("TauDecayModeFinding", &vars.tauDecayModeFinding);
	tree->SetBranchAddress("TauQuality", &vars.tauQuality);
	tree->SetBranchAddress("TauIsolation3Hits", &vars.tauIsolation3Hits);
	tree->SetBranchAddress("TauLooseIsolation3Hits", &vars.tauLooseIsolation3Hits);
	tree->SetBranchAddress("TauMediumIsolation3Hits", &vars.tauMediumIsolation3Hits);
	tree->SetBranchAddress("TauTightIsolation3Hits", &vars.tauTightIsolation3Hits);
	tree->SetBranchAddress("TauAgainstMuonLoose", &vars.tauAgainstMuonLoose);
	tree->SetBranchAddress("TauAgainstMuonTight", &vars.tauAgainstMuonTight);
	tree->SetBranchAddress("TauAgainstElectronLoose", &vars.tauAgainstElectronLoose);
	tree->SetBranchAddress("TauAgainstElectronLooseMVA3", &vars.tauAgainstElectronLooseMVA3);
	tree->SetBranchAddress("TauAgainstElectronMediumMVA3", &vars.tauAgainstElectronMediumMVA3);
	tree->SetBranchAddress("TauAgainstElectronTightMVA3", &vars.tauAgainstElectronTightMVA3);
	tree->SetBranchAddress("TauGenPt", &vars.tauGenPt);
	tree->SetBranchAddress("TauGenEta", &vars.tauGenEta);
	tree->SetBranchAddress("TauGenPhi", &vars.tauGenPhi);
	tree->SetBranchAddress("TauGenE", &vars.tauGenE);
	tree->SetBranchAddress("TauGenFinalPt", &vars.tauGenFinalPt);
	tree->SetBranchAddress("TauGenFinalEta", &vars.tauGenFinalEta);
	tree->SetBranchAddress("TauGenFinalPhi", &vars.tauGenFinalPhi);
	tree->SetBranchAddress("TauGenFinalE", &vars.tauGenFinalE);
	tree->SetBranchAddress("TauGenCharge", &vars.tauGenCharge);
	tree->SetBranchAddress("TauGenPtVis", &vars.tauGenPtVis);
	tree->SetBranchAddress("TauGenEtaVis", &vars.tauGenEtaVis);
	tree->SetBranchAddress("TauGenPhiVis", &vars.tauGenPhiVis);
	tree->SetBranchAddress("TauGenEVis", &vars.tauGenEVis);

	tree->SetBranchAddress("NJets", &vars.nJets);
	tree->SetBranchAddress("JetPt", &vars.jetPt);
	tree->SetBranchAddress("JetEta", &vars.jetEta);
	tree->SetBranchAddress("JetPhi", &vars.jetPhi);
	tree->SetBranchAddress("JetEnergy", &vars.jetEnergy);
	tree->SetBranchAddress("JetBTag", &vars.jetBTag);

	if(tree->GetBranch("SVfitMass"))
		tree->SetBranchAddress("SVfitMass", &vars.svfitMass);

	TH1F* h_nPV = dynamic_cast<TH1F*>(file->Get((directory + "/h_nPV").c_str()));
	out.add_events(static_cast<unsigned int>(h_nPV->Integral()));
	out.add_pv(h_nPV);

	if(!is_data)
	{
		TH1F* h_nTrueInteractions = dynamic_cast<TH1F*>(file->Get((directory + "/h_nTrueInteractions").c_str()));
		out.add_interactions(h_nTrueInteractions);
	}

	unsigned int Filter = 0;
	const unsigned int CountTotal = tree->GetEntries();
	for(int i = 0; i < CountTotal; ++i)
	{
		tree->GetEntry(i);

		// Objects
		const TLepton lepton(vars);
		const Tau tau(vars);

		// PDG filter for embedded samples; only take muons from DY -> mumu
		if(!is_data && is_embedded)
			if(abs(vars.origGenPosMuonPDG) != 13 || abs(vars.origGenNegMuonPDG) != 13)
				continue;

		// There are some rare events for which the generator-level particles
		// do not match the original muons.
		if(is_embedded)
		{
			assert(lepton.hasOrigRec && tau.hasOrigRec);
			assert(lepton.hasGen && tau.hasGen);

			TLorentzVector orig = lepton.p4orig + tau.p4orig;
			TLorentzVector gen = lepton.p4genInitial + tau.p4genInitial;
			TLorentzVector diff = orig - gen;
			if(fabs(diff.Px()) > 1e-3 || fabs(diff.Py()) > 1e-3 || fabs(diff.Pz()) > 1e-3 || fabs(diff.E()) > 1e-3)
			{
				++Filter;
				continue;
			}
		}

		// Weighting
		const double zmumuEvtSelEff = ((zmumuWeight == NULL) ? 1. : zmumuWeight->getWeight(lepton, tau));
		double weight = vars.minVisPtFilterWeight * vars.tauSpinWeight * zmumuEvtSelEff;

		// Output
		out.fill(vars, lepton, tau, weight);
	}

	std::cout << "Filtered: " << Filter << "/" << CountTotal << std::endl;
	delete file;
}

template<typename TLepton>
void process_files(const std::vector<std::string>& files, Out<TLepton>& out, bool is_data, bool is_embedded, bool is_rh_embedded, const ZmumuEvtSelEffCorrWeight* zmumuWeight)
{
	for(std::vector<std::string>::const_iterator iter = files.begin(); iter != files.end(); ++iter)
	{
		std::cout << "\t[" << (iter - files.begin() + 1) << "/" << files.size() << "] " << *iter << "..." << std::endl;
		process_file(iter->c_str(), out, is_data, is_embedded, is_rh_embedded, zmumuWeight);
	}
}

void process_dataset(const char* dataset, const std::vector<std::string>& files)
{
	const bool is_data = (strstr(dataset, "Data") != NULL);
	const bool is_embedded = (strstr(dataset, "Embedded") != NULL);
	const bool is_rh_embedded = (strstr(dataset, "RHEmbedded") != NULL);

	std::string channel;
	if(strstr(dataset, "etau") != NULL) channel = "etau";
	else if(strstr(dataset, "mutau") != NULL) channel = "mutau";
	else throw std::runtime_error("Unknown channel for dataset \"" + std::string(dataset) + "\"");

	std::stringstream out_filename_stream;
	out_filename_stream << "out-";
	if(strncmp(dataset, "embedtest", 9) != 0) out_filename_stream << "embedtest_";
	out_filename_stream << dataset << "-" << channel;
	out_filename_stream << ".root";
	const std::string out_filename = out_filename_stream.str();

	std::auto_ptr<ZmumuEvtSelEffCorrWeight> zmumuWeight;
	if(is_embedded) zmumuWeight.reset(new ZmumuEvtSelEffCorrWeight("../data/ZmumuEvtSelEffCorrWeights.root"));

	std::cout << dataset << "... => " << out_filename << std::endl;
	std::cout << "\tData: " << is_data << std::endl;
	std::cout << "\tEmbedded: " << is_embedded << std::endl;
	std::cout << "\tRHEmbedded: " << is_rh_embedded << std::endl;

	TFile* outfile = new TFile(out_filename.c_str(), "RECREATE");

	if(channel == "etau")
	{
		Out<Electron> out(outfile, is_data, is_embedded);
		process_files(files, out, is_data, is_embedded, is_rh_embedded, zmumuWeight.get());
	}
	else if(channel == "mutau")
	{
		Out<Muon> out(outfile, is_data, is_embedded);
		process_files(files, out, is_data, is_embedded, is_rh_embedded, zmumuWeight.get());
	}
	else
	{
		assert(false);
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
