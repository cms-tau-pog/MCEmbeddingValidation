#include <TTree.h>
#include <TFile.h>
#include <TMath.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <Math/VectorUtil.h>
#include <sstream>
#include <memory>

#include "lumiFilter.h"

struct TreeVars {
	TreeVars() {
		run = lumi = event = ~0u;
		weight = -1.0f;
		nPV = 0; nTrueInteractions = -1.0f;
		pfMet = pfMetPhi = mvaPfMet = -1.0f;
		hltMu17Mu8 = false;
		PosMuon.Pt = PosMuon.Eta = PosMuon.Phi = PosMuon.E = -1.0f;
		PosMuon.HltMu17Mu8Leg8 = PosMuon.HltMu17Mu8Leg17 = false;
		PosMuon.TrackerIso03 = PosMuon.EcalIso03 = PosMuon.HcalIso03 = PosMuon.ChargedHadronPfIso04 = -1.0f;
		PosMuon.Quality = false;
		PosMuon.GenChargedLepPDG = 0;

		NegMuon.Pt = NegMuon.Eta = NegMuon.Phi = NegMuon.E = -1.0f;
		NegMuon.HltMu17Mu8Leg8 = NegMuon.HltMu17Mu8Leg17 = false;
		NegMuon.TrackerIso03 = NegMuon.EcalIso03 = NegMuon.HcalIso03 = NegMuon.ChargedHadronPfIso04 = -1.0f;
		NegMuon.Quality = false;
		NegMuon.GenChargedLepPDG = 0;

		nJets = 0;
		nTracks5 = nTracks10 = nTracks20 = nTracks30 = nTracks40 = 0;
		nGlobalMuons = nStandaloneMuons = nGlobalMuons = 0;
	}

	unsigned int run;
	unsigned int lumi;
	unsigned int event;
	float weight;

	unsigned int nPV;
	float nTrueInteractions;

	bool hltMu17Mu8;

	float pfMet;
	float pfMetType1;
	float pfMetPhi;
	float mvaPfMet;

	struct Muon {
		float Pt;
		float Eta;
		float Phi;
		float E;

		bool HltMu17Mu8Leg8;
		bool HltMu17Mu8Leg17;

		float TrackerIso03;
		float EcalIso03;
		float HcalIso03;

		float ChargedHadronPfIso04;
		bool Quality;

		int GenChargedLepPDG;
	};

	Muon PosMuon;
	Muon NegMuon;

	unsigned int nJets;
	float jetPt[20];
	float jetEta[20];
	float jetPhi[20];
	float jetEnergy[20];
	float jetBTag[20];

	unsigned int nTracks5;
	unsigned int nTracks10;
	unsigned int nTracks20;
	unsigned int nTracks30;
	unsigned int nTracks40;

	unsigned int nGlobalMuons;
	unsigned int nStandaloneMuons;
	unsigned int nPFMuons;
};

class Muon: public TreeVars::Muon
{
public:
	Muon(const TreeVars::Muon& base):
		TreeVars::Muon(base)
	{
		p4.SetPtEtaPhiE(Pt, Eta, Phi, E);
	}

	float RelCombIso() const
	{
		return (TrackerIso03 + EcalIso03 + HcalIso03) / Pt;
	}

	float RelChargedHadronPfIso() const
	{
		return (ChargedHadronPfIso04) / Pt;
	}

	TLorentzVector p4;
};

class Jet
{
public:
	Jet(const TreeVars& base, unsigned int n)
	{
		p4.SetPtEtaPhiE(base.jetPt[n], base.jetEta[n], base.jetPhi[n], base.jetEnergy[n]);
		btag = base.jetBTag[n];
	}

	TLorentzVector p4;
	float btag;
};

class IDIsoWeights
{
public:
	IDIsoWeights(const char* id_filename, const char* iso_filename)
	{
		f_id = new TFile(id_filename, "READ");
		f_iso = new TFile(iso_filename, "READ");

		idAbsEta0_0To0_9 = dynamic_cast<TGraph*>(f_id->Get("DATA_over_MC_Tight_pt_abseta<0.9"));
		idAbsEta0_9To1_2 = dynamic_cast<TGraph*>(f_id->Get("DATA_over_MC_Tight_pt_abseta0.9-1.2"));
		idAbsEta1_2To2_1 = dynamic_cast<TGraph*>(f_id->Get("DATA_over_MC_Tight_pt_abseta1.2-2.1"));
		idAbsEta2_1To2_4 = dynamic_cast<TGraph*>(f_id->Get("DATA_over_MC_Tight_pt_abseta2.1-2.4"));

		isoAbsEta0_0To0_9 = dynamic_cast<TGraph*>(f_iso->Get("DATA_over_MC_combRelIsoPF04dBeta<02_Tight_pt_abseta<0.9"));
		isoAbsEta0_9To1_2 = dynamic_cast<TGraph*>(f_iso->Get("DATA_over_MC_combRelIsoPF04dBeta<02_Tight_pt_abseta0.9-1.2"));
		isoAbsEta1_2To2_1 = dynamic_cast<TGraph*>(f_iso->Get("DATA_over_MC_combRelIsoPF04dBeta<02_Tight_pt_abseta1.2-2.1"));
		isoAbsEta2_1To2_4 = dynamic_cast<TGraph*>(f_iso->Get("DATA_over_MC_combRelIsoPF04dBeta<02_Tight_pt_abseta2.1-2.4"));

		assert(idAbsEta0_0To0_9 != NULL);
		assert(idAbsEta0_9To1_2 != NULL);
		assert(idAbsEta1_2To2_1 != NULL);
		assert(idAbsEta2_1To2_4 != NULL);

		assert(isoAbsEta0_0To0_9 != NULL);
		assert(isoAbsEta0_9To1_2 != NULL);
		assert(isoAbsEta1_2To2_1 != NULL);
		assert(isoAbsEta2_1To2_4 != NULL);
	}

	~IDIsoWeights()
	{
		delete f_id;
		delete f_iso;
	}

	double GetIDWeight(const Muon& muon)
	{
		if(fabs(muon.Eta) < 0.9) return idAbsEta0_0To0_9->Eval(muon.Pt);
		else if(fabs(muon.Eta) < 1.2) return idAbsEta0_9To1_2->Eval(muon.Pt);
		else if(fabs(muon.Eta) < 2.1) return idAbsEta1_2To2_1->Eval(muon.Pt);
		else return idAbsEta2_1To2_4->Eval(muon.Pt);
	}

	double GetIsoWeight(const Muon& muon)
	{
		if(fabs(muon.Eta) < 0.9) return isoAbsEta0_0To0_9->Eval(muon.Pt);
		else if(fabs(muon.Eta) < 1.2) return isoAbsEta0_9To1_2->Eval(muon.Pt);
		else if(fabs(muon.Eta) < 2.1) return isoAbsEta1_2To2_1->Eval(muon.Pt);
		else return isoAbsEta2_1To2_4->Eval(muon.Pt);
	}

	double GetWeight(const Muon& muon)
	{
		return GetIDWeight(muon) * GetIsoWeight(muon);
	}

private:
	TFile* f_id;
	TFile* f_iso;

	TGraph* idAbsEta0_0To0_9;
	TGraph* idAbsEta0_9To1_2;
	TGraph* idAbsEta1_2To2_1;
	TGraph* idAbsEta2_1To2_4;

	TGraph* isoAbsEta0_0To0_9;
	TGraph* isoAbsEta0_9To1_2;
	TGraph* isoAbsEta1_2To2_1;
	TGraph* isoAbsEta2_1To2_4;
};

class TriggerWeights
{
public:
	TriggerWeights(const char* filename)
	{
		f = new TFile(filename, "READ");

		hltMu17Mu8_Mu1_10To20_Mu2_20ToInfty = dynamic_cast<TH2F*>(f->Get("DATA_over_MC_Mu17Mu8_Tight_Mu1_10To20_&_Mu2_20ToInfty_with_STAT_uncrt"));
		hltMu17Mu8_Mu1_20ToInfty_Mu2_20ToInfty = dynamic_cast<TH2F*>(f->Get("DATA_over_MC_Mu17Mu8_Tight_Mu1_20ToInfty_&_Mu2_20ToInfty_with_STAT_uncrt"));

		assert(hltMu17Mu8_Mu1_10To20_Mu2_20ToInfty != NULL);
		assert(hltMu17Mu8_Mu1_20ToInfty_Mu2_20ToInfty != NULL);
	}

	~TriggerWeights() { delete f; }

	double GetWeight(const Muon& mu8, const Muon& mu17)
	{
		if(mu8.Pt > 20. && mu17.Pt > 20.)
		{
			const double minEta = std::min(std::min(fabs(mu8.Eta), fabs(mu17.Eta)), 2.3999);
			const double maxEta = std::min(std::max(fabs(mu8.Eta), fabs(mu17.Eta)), 2.3999);
			const int xbin = hltMu17Mu8_Mu1_20ToInfty_Mu2_20ToInfty->GetXaxis()->FindBin(maxEta);
			const int ybin = hltMu17Mu8_Mu1_20ToInfty_Mu2_20ToInfty->GetYaxis()->FindBin(minEta);
			const double value = hltMu17Mu8_Mu1_20ToInfty_Mu2_20ToInfty->GetBinContent(xbin, ybin);
			assert(value >= 0.5);
			return value;
		}
		else
		{
			const double eta1 = std::min(fabs(mu8.Eta), 2.3999);
			const double eta2 = std::min(fabs(mu17.Eta), 2.3999);
			const int xbin = hltMu17Mu8_Mu1_10To20_Mu2_20ToInfty->GetXaxis()->FindBin(eta1);
			const int ybin = hltMu17Mu8_Mu1_10To20_Mu2_20ToInfty->GetYaxis()->FindBin(eta2);
			const double value = hltMu17Mu8_Mu1_10To20_Mu2_20ToInfty->GetBinContent(xbin, ybin);
			assert(value >= 0.5);
			return value;
		}
	}
protected:
	TFile* f;
	TH2F* hltMu17Mu8_Mu1_10To20_Mu2_20ToInfty;
	TH2F* hltMu17Mu8_Mu1_20ToInfty_Mu2_20ToInfty;
};

struct DiMuonAcceptance
{
	DiMuonAcceptance(const Muon& first, const Muon& second)
	{
		if(first.p4.Pt() < 8 || second.p4.Pt() < 8) accepted = false;
		else if(first.p4.Pt() < 17 && second.p4.Pt() < 17) accepted = false;
		else if(fabs(first.p4.Eta()) > 2.5) accepted = false;
		else if(fabs(second.p4.Eta()) > 2.5) accepted = false;
		else if((first.p4 + second.p4).M() < 50.0) accepted = false; // There is some really low energy stuff that we don't want...
		else accepted = true;

		ided = first.Quality && second.Quality;
		isolated = first.RelChargedHadronPfIso() < 0.10 && second.RelChargedHadronPfIso() < 0.10;
		weakly_isolated = first.RelChargedHadronPfIso() > 0.30 && first.RelChargedHadronPfIso() < 1.00 && second.RelChargedHadronPfIso() > 0.30 && second.RelChargedHadronPfIso() < 1.00; // For QCD estimation
		masscut = (first.p4 + second.p4).M() > 50.0;
	}

	bool accepted;
	bool ided;
	bool isolated;
	bool weakly_isolated;
	bool masscut;
};

struct JetCategorization
{
	JetCategorization(const std::vector<Jet>& jets)
	{
		if(isVBF(jets)) category = VBF;
		else if(isBTag(jets)) category = BTAG;
		else if(is1Jet(jets)) category = _1JET;
		else category = _0JET;
	}

	bool isVBF(const std::vector<Jet>& jets) const
	{
		for(std::vector<Jet>::const_iterator iter = jets.begin(); iter != jets.end(); ++iter)
		{
			if(iter->p4.Pt() > 30.0)
			{
				for(std::vector<Jet>::const_iterator iter2 = iter + 1; iter2 != jets.end(); ++iter2)
				{
					if(iter2->p4.Pt() > 30.0 && fabs(iter->p4.Eta() - iter2->p4.Eta()) > 3.5 && (iter->p4 + iter2->p4).M() > 350 && isCentralJetVeto(jets, iter, iter2))
						return true;
				}
			}
		}

		return false;
	}

	bool isCentralJetVeto(const std::vector<Jet>& jets, const std::vector<Jet>::const_iterator& first_jet, const std::vector<Jet>::const_iterator& second_jet) const
	{
		const double min_eta = std::min(first_jet->p4.Eta(), second_jet->p4.Eta());
		const double max_eta = std::max(first_jet->p4.Eta(), second_jet->p4.Eta());

		for(std::vector<Jet>::const_iterator iter = jets.begin(); iter != jets.end(); ++iter)
		{
			if(iter == first_jet || iter == second_jet) continue;
			if(iter->p4.Pt() > 20.0 && iter->p4.Eta() > min_eta && iter->p4.Eta() < max_eta)
				return false;
		}

		return true;
	}

	bool isBTag(const std::vector<Jet>& jets) const
	{
		for(std::vector<Jet>::const_iterator iter = jets.begin(); iter != jets.end(); ++iter)
			if(iter->p4.Pt() > 20.0 && fabs(iter->p4.Eta()) < 2.4 && iter->btag > 0.679)
				return true;
		return false;
	}

	bool is1Jet(const std::vector<Jet>& jets) const
	{
		for(std::vector<Jet>::const_iterator iter = jets.begin(); iter != jets.end(); ++iter)
			if(iter->p4.Pt() > 30.0)
				return true;
		return false;
	}

	enum { VBF, BTAG, _1JET, _0JET } category;
};

struct Histograms {
	Histograms(TDirectory* dir, const char* name) {
		TDirectory* newdir = dir->mkdir(name);
		newdir->cd();

		nPV = new TH1D("nPV", "Number of primary vertices;N_{PV};Entries", 50, -0.5, 49.5);
		nPV_u = new TH1D("nPV_u", "Number of primary vertices, unweighted;N_{PV};Entries", 50, -0.5, 49.5);

		posMuonPt = new TH1D("posMuonPt", "Positive muon transverse momentum;p_{T};Entries/0.1 GeV", 1000, 0.0, 100.0);
		posMuonEta = new TH1D("posMuonEta", "Positive muon pseudo-rapidity;#eta;Entries/0.1", 180, -9.0, 9.0);
		posMuonPhi = new TH1D("posMuonPhi", "Positive muon azimuthal angle;#phi;Entries", 100, -M_PI, M_PI);
		posMuonIso = new TH1D("posMuonIso", "Positive muon relative PF isolation;I_{rel}^{PF}", 100, 0.0, 1.0);

		negMuonPt = new TH1D("negMuonPt", "Positive muon transverse momentum;p_{T};Entries/0.1 GeV", 1000, 0.0, 100.0);
		negMuonEta = new TH1D("negMuonEta", "Positive muon pseudo-rapidity;#eta;Entries/0.1", 180, -9.0, 9.0);
		negMuonPhi = new TH1D("negMuonPhi", "Positive muon azimuthal angle;#phi;Entries", 100, -M_PI, M_PI);
		negMuonIso = new TH1D("negMuonIso", "Positive muon relative PF isolation;I_{rel}^{PF}", 100, 0.0, 1.0);

		leadMuonPt = new TH1D("leadMuonPt", "Leading muon transverse momentum;p_{T};Entries/0.1 GeV", 1000, 0.0, 100.0);
		leadMuonEta = new TH1D("leadMuonEta", "Leading muon pseudo-rapidity;#eta;Entries/0.1", 180, -9.0, 9.0);
		leadMuonPhi = new TH1D("leadMuonPhi", "Leading muon azimuthal angle;#phi;Entries", 100, -M_PI, M_PI);
		leadMuonIso = new TH1D("leadMuonIso", "Leading muon relative PF isolation;I_{rel}^{PF}", 100, 0.0, 1.0);

		subMuonPt = new TH1D("subMuonPt", "Subleading muon transverse momentum;p_{T};Entries/0.1 GeV", 1000, 0.0, 100.0);
		subMuonEta = new TH1D("subMuonEta", "Subleading muon pseudo-rapidity;#eta;Entries/0.1", 180, -9.0, 9.0);
		subMuonPhi = new TH1D("subMuonPhi", "Subleading muon azimuthal angle;#phi;Entries", 100, -M_PI, M_PI);
		subMuonIso = new TH1D("subMuonIso", "Subleading muon relative PF isolation;I_{rel}^{PF}", 100, 0.0, 1.0);

		diMuonPt = new TH1D("diMuonPt", "Di-muon transverse momentum;p_{T};Entries/0.1 GeV", 1000, 0.0, 100.0);
		diMuonEta = new TH1D("diMuonEta", "Di-muon pseudo-rapidity;#eta;Entries/0.1", 180, -9.0, 9.0);
		diMuonPhi = new TH1D("diMuonPhi", "Di-muon azimuthal angle;#phi;Entries", 100, -M_PI, M_PI);
		diMuonDeltaPhi = new TH1D("diMuonDeltaPhi", "Di-muon #Delta #phi;#Delta Phi;Entries", 100, 0, M_PI);
		diMuonDeltaR = new TH1D("diMuonDeltaR", "Di-muon #Delta R;#Delta R;Entries/0.05", 200, 0.0, 10.0);
		diMuonMass = new TH1D("diMuonMass", "Di-muon invariant mass;m_{#mu#mu};Entries/0.1 GeV", 2000, 0.0, 200.0);

		nJets30 = new TH1D("nJets30", "Number of jets with p_{T} > 30 GeV;N_{jets};Entries", 20, -0.5, 19.5);
		bJetPt = new TH1D("bJetPt", "Transverse momentum of leading b-jet;p_{T};Entries/0.1 GeV", 1800, 20.0, 200.0);

		nTracks5 = new TH1D("nTracks5", "Number of tracks with p_{T} > 5 GeV;N_{tracks};Entries", 100, -0.5, 99.5);
		nTracks10 = new TH1D("nTracks10", "Number of tracks with p_{T} > 10 GeV;N_{tracks};Entries", 50, -0.5, 49.5);
		nTracks20 = new TH1D("nTracks20", "Number of tracks with p_{T} > 20 GeV;N_{tracks};Entries", 50, -0.5, 49.5);
		nTracks30 = new TH1D("nTracks30", "Number of tracks with p_{T} > 30 GeV;N_{tracks};Entries", 40, -0.5, 39.5);
		nTracks40 = new TH1D("nTracks40", "Number of tracks with p_{T} > 40 GeV;N_{tracks};Entries", 20, -0.5, 19.5);

		nGlobalMuons = new TH1D("nGlobalMuons", "Number of global muons;N_{muons};Entries", 10, -0.5, 9.5);
		nStandaloneMuons = new TH1D("nStandaloneMuons", "Number of standalone muons;N_{muons};Entries", 10, -0.5, 9.5);
		nPFMuons = new TH1D("nPFMuons", "Number of PF muons;N_{muons};Entries", 10, -0.5, 9.5);

		pfMet = new TH1D("pfMet", "Particle Flow missing transverse energy;E_{T}^{miss};Entries/0.1 GeV", 2000, 0.0, 200.0);
		pfMetType1 = new TH1D("pfMetType1", "Particle Flow missing transverse energy, type1 corrected;E_{T}^{miss};Entries/0.1 GeV", 2000, 0.0, 200.0);
		mvaPfMet = new TH1D("mvaPfMet", "MVA Particle Flow missing transverse energy;E_{T}^{miss};Entries/0.1 GeV", 2000, 0.0, 200.0);

		dir->cd();
	}

	void fill(const Muon& posMuon, const Muon& negMuon, unsigned int nJets, const TreeVars& vars, float weight)
	{
		const Muon& leadMuon = posMuon.p4.Pt() > negMuon.p4.Pt() ? posMuon : negMuon;
		const Muon& subMuon = posMuon.p4.Pt() <= negMuon.p4.Pt() ? posMuon : negMuon;

		nPV->Fill(vars.nPV, weight);
		nPV_u->Fill(vars.nPV);

		posMuonPt->Fill(posMuon.p4.Pt(), weight);
		posMuonEta->Fill(posMuon.p4.Eta(), weight);
		posMuonPhi->Fill(posMuon.p4.Phi(), weight);
		posMuonIso->Fill(posMuon.RelChargedHadronPfIso(), weight);

		negMuonPt->Fill(negMuon.p4.Pt(), weight);
		negMuonEta->Fill(negMuon.p4.Eta(), weight);
		negMuonPhi->Fill(negMuon.p4.Phi(), weight);
		negMuonIso->Fill(negMuon.RelChargedHadronPfIso(), weight);

		leadMuonPt->Fill(leadMuon.p4.Pt(), weight);
		leadMuonEta->Fill(leadMuon.p4.Eta(), weight);
		leadMuonPhi->Fill(leadMuon.p4.Phi(), weight);
		leadMuonIso->Fill(leadMuon.RelChargedHadronPfIso(), weight);

		subMuonPt->Fill(subMuon.p4.Pt(), weight);
		subMuonEta->Fill(subMuon.p4.Eta(), weight);
		subMuonPhi->Fill(subMuon.p4.Phi(), weight);
		subMuonIso->Fill(subMuon.RelChargedHadronPfIso(), weight);

		diMuonPt->Fill( (posMuon.p4 + negMuon.p4).Pt(), weight);
		diMuonEta->Fill( (posMuon.p4 + negMuon.p4).Eta(), weight);
		diMuonPhi->Fill( (posMuon.p4 + negMuon.p4).Phi(), weight);
		diMuonDeltaPhi->Fill(ROOT::Math::VectorUtil::DeltaPhi(posMuon.p4, negMuon.p4), weight);
		diMuonDeltaR->Fill(ROOT::Math::VectorUtil::DeltaR(posMuon.p4, negMuon.p4), weight);
		diMuonMass->Fill( (posMuon.p4 + negMuon.p4).M(), weight);

		nJets30->Fill(nJets, weight);
		float bJetPtVal = -1.;
		for(unsigned int i = 0; i < vars.nJets; ++i)
			if(vars.jetPt[i] > 20. && fabs(vars.jetEta[i]) < 2.4 && vars.jetBTag[i] > 0.679)
				{ bJetPtVal = vars.jetPt[i]; break; }
		if(bJetPtVal > 0.) bJetPt->Fill(bJetPtVal, weight);

		nTracks5->Fill(vars.nTracks5, weight);
		nTracks10->Fill(vars.nTracks10, weight);
		nTracks20->Fill(vars.nTracks20, weight);
		nTracks30->Fill(vars.nTracks30, weight);
		nTracks40->Fill(vars.nTracks40, weight);

		nGlobalMuons->Fill(vars.nGlobalMuons, weight);
		nStandaloneMuons->Fill(vars.nStandaloneMuons, weight);
		nPFMuons->Fill(vars.nPFMuons, weight);

		pfMet->Fill(vars.pfMet, weight);
		pfMetType1->Fill(vars.pfMetType1, weight);
		mvaPfMet->Fill(vars.mvaPfMet, weight);
	}

	TH1D* nPV;
	TH1D* nPV_u;

	TH1D* posMuonPt;
	TH1D* posMuonEta;
	TH1D* posMuonPhi;
	TH1D* posMuonIso;

	TH1D* negMuonPt;
	TH1D* negMuonEta;
	TH1D* negMuonPhi;
	TH1D* negMuonIso;

	TH1D* leadMuonPt;
	TH1D* leadMuonEta;
	TH1D* leadMuonPhi;
	TH1D* leadMuonIso;

	TH1D* subMuonPt;
	TH1D* subMuonEta;
	TH1D* subMuonPhi;
	TH1D* subMuonIso;

	TH1D* diMuonPt;
	TH1D* diMuonEta;
	TH1D* diMuonPhi;
	TH1D* diMuonDeltaPhi;
	TH1D* diMuonDeltaR;
	TH1D* diMuonMass;

	TH1D* nJets30;
	TH1D* bJetPt;

	TH1D* nTracks5;
	TH1D* nTracks10;
	TH1D* nTracks20;
	TH1D* nTracks30;
	TH1D* nTracks40;

	TH1D* nGlobalMuons;
	TH1D* nStandaloneMuons;
	TH1D* nPFMuons;

	TH1D* pfMet;
	TH1D* pfMetType1;
	TH1D* mvaPfMet;
};

struct Category
{
	Category(TDirectory* dir, const char* name)
	{
		TDirectory* newdir = dir->mkdir(name);
		newdir->cd();

		initial.reset(new Histograms(newdir, "initial"));
		accepted.reset(new Histograms(newdir, "accepted"));
		ided.reset(new Histograms(newdir, "ided"));
		isolated.reset(new Histograms(newdir, "isolated"));
		masscut.reset(new Histograms(newdir, "masscut"));
		isolated_masscut.reset(new Histograms(newdir, "isolated_masscut"));

		weakly_isolated.reset(new Histograms(newdir, "weakly_isolated"));
		weakly_isolated_masscut.reset(new Histograms(newdir, "weakly_isolated_masscut"));

		dir->cd();
	}

	void fill(const Muon& muon1, const Muon& muon2, unsigned int nJets, const TreeVars& vars, float weight)
	{
		initial->fill(muon1, muon2, nJets, vars, weight);

		DiMuonAcceptance acceptance(muon1, muon2);
		if(acceptance.accepted) accepted->fill(muon1, muon2, nJets, vars, weight);
		if(acceptance.accepted && acceptance.ided) ided->fill(muon1, muon2, nJets, vars, weight);
		if(acceptance.accepted && acceptance.ided && acceptance.isolated) isolated->fill(muon1, muon2, nJets, vars, weight);
		if(acceptance.accepted && acceptance.ided && acceptance.masscut) masscut->fill(muon1, muon2, nJets, vars, weight);
		if(acceptance.accepted && acceptance.ided && acceptance.isolated && acceptance.masscut) isolated_masscut->fill(muon1, muon2, nJets, vars, weight);
		if(acceptance.accepted && acceptance.ided && acceptance.weakly_isolated) weakly_isolated->fill(muon1, muon2, nJets, vars, weight);
		if(acceptance.accepted && acceptance.ided && acceptance.weakly_isolated && acceptance.masscut) weakly_isolated_masscut->fill(muon1, muon2, nJets, vars, weight);
	}

	std::auto_ptr<Histograms> initial;
	std::auto_ptr<Histograms> accepted;
	std::auto_ptr<Histograms> ided;
	std::auto_ptr<Histograms> isolated;
	std::auto_ptr<Histograms> masscut;
	std::auto_ptr<Histograms> isolated_masscut;
	std::auto_ptr<Histograms> weakly_isolated;
	std::auto_ptr<Histograms> weakly_isolated_masscut;
};

struct Out
{
	Out(TDirectory* dir, bool is_data):
		inclusive(dir, "inclusive"),
		inclusive_unweighted(dir, "inclusive_unweighted"),
		vbf(dir, "vbf"),
		btag(dir, "btag"),
		_1jet(dir, "1jet"),
		_0jet(dir, "0jet")
	{
		nEvents = new TH1D("nEvents", "Total number of events", 1, 0, 1);
		nPV = new TH1D("nPV", "Number of Primary Vertices", 50, -0.5, 49.5);
		if(!is_data) nTrueInteractions = new TH1D("nTrueInteractions", "Number of True Interactions", 100, -0.5, 99.5);
	}

	void fill(const Muon& muon1, const Muon& muon2, const std::vector<Jet>& jets, const TreeVars& vars, float weight)
	{
		const JetCategorization jetCategorization(jets);

		unsigned int nJets = 0;
		for(unsigned int i = 0; i < jets.size(); ++i)
			if(jets[i].p4.Pt() > 30)
				++nJets;

		inclusive.fill(muon1, muon2, nJets, vars, weight);
		inclusive_unweighted.fill(muon1, muon2, nJets, vars, 1.0);
		switch(jetCategorization.category)
		{
		case JetCategorization::VBF: vbf.fill(muon1, muon2, nJets, vars, weight); break;
		case JetCategorization::BTAG: btag.fill(muon1, muon2, nJets, vars, weight); break;
		case JetCategorization::_1JET: _1jet.fill(muon1, muon2, nJets, vars, weight); break;
		case JetCategorization::_0JET: _0jet.fill(muon1, muon2, nJets, vars, weight); break;
		default: assert(false); break;
		}
	}

	void add_events(unsigned int n)
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

	Category inclusive;
	Category inclusive_unweighted;
	Category vbf;
	Category btag;
	Category _1jet;
	Category _0jet;
	TH1D* nEvents;
	TH1D* nPV;
	TH1D* nTrueInteractions;
};

void process_file(const char* filename, Out& out, bool is_data, bool is_ss, int pdg_filter, TH1D* puWeightHisto, TriggerWeights* trigWeights, IDIsoWeights* idIsoWeights)
{
	TreeVars vars;

	std::string directory = is_ss ? "diMuonsSS" : "diMuons";
	TFile* file = new TFile(filename, "READ");
	TTree* tree = dynamic_cast<TTree*>(file->Get((directory + "/MuonTree").c_str()));

	tree->SetBranchAddress("Run", &vars.run);
	tree->SetBranchAddress("Lumi", &vars.lumi);
	tree->SetBranchAddress("Event", &vars.event);

	tree->SetBranchAddress("Weight", &vars.weight);
	tree->SetBranchAddress("NPV", &vars.nPV);
	if(!is_data) tree->SetBranchAddress("NTrueInteractions", &vars.nTrueInteractions);
	tree->SetBranchAddress("HltMu17Mu8", &vars.hltMu17Mu8);

	tree->SetBranchAddress("PfMet", &vars.pfMet);
	tree->SetBranchAddress("PfType1CorrectedMet", &vars.pfMetType1);
	tree->SetBranchAddress("PfMetPhi", &vars.pfMetPhi);
	//tree->SetBranchAddress("MvaPfMet", &vars.mvaPfMet);

	tree->SetBranchAddress("PosMuonPt", &vars.PosMuon.Pt);
	tree->SetBranchAddress("PosMuonEta", &vars.PosMuon.Eta);
	tree->SetBranchAddress("PosMuonPhi", &vars.PosMuon.Phi);
	tree->SetBranchAddress("PosMuonE", &vars.PosMuon.E);
	tree->SetBranchAddress("PosMuonTrackerIso03", &vars.PosMuon.TrackerIso03);
	tree->SetBranchAddress("PosMuonEcalIso03", &vars.PosMuon.EcalIso03);
	tree->SetBranchAddress("PosMuonHcalIso03", &vars.PosMuon.HcalIso03);
	tree->SetBranchAddress("PosMuonChargedHadronPfIso04", &vars.PosMuon.ChargedHadronPfIso04);
	tree->SetBranchAddress("PosMuonQuality", &vars.PosMuon.Quality);
	tree->SetBranchAddress("PosMuonHltMu17Mu8Leg8", &vars.PosMuon.HltMu17Mu8Leg8);
	tree->SetBranchAddress("PosMuonHltMu17Mu8Leg17", &vars.PosMuon.HltMu17Mu8Leg17);
	if(!is_data) tree->SetBranchAddress("GenPosChargedLepPDG", &vars.PosMuon.GenChargedLepPDG);

	tree->SetBranchAddress("NegMuonPt", &vars.NegMuon.Pt);
	tree->SetBranchAddress("NegMuonEta", &vars.NegMuon.Eta);
	tree->SetBranchAddress("NegMuonPhi", &vars.NegMuon.Phi);
	tree->SetBranchAddress("NegMuonE", &vars.NegMuon.E);
	tree->SetBranchAddress("NegMuonTrackerIso03", &vars.NegMuon.TrackerIso03);
	tree->SetBranchAddress("NegMuonEcalIso03", &vars.NegMuon.EcalIso03);
	tree->SetBranchAddress("NegMuonHcalIso03", &vars.NegMuon.HcalIso03);
	tree->SetBranchAddress("NegMuonChargedHadronPfIso04", &vars.NegMuon.ChargedHadronPfIso04);
	tree->SetBranchAddress("NegMuonQuality", &vars.NegMuon.Quality);
	tree->SetBranchAddress("NegMuonHltMu17Mu8Leg8", &vars.NegMuon.HltMu17Mu8Leg8);
	tree->SetBranchAddress("NegMuonHltMu17Mu8Leg17", &vars.NegMuon.HltMu17Mu8Leg17);
	if(!is_data) tree->SetBranchAddress("GenNegChargedLepPDG", &vars.NegMuon.GenChargedLepPDG);

	tree->SetBranchAddress("nJets", &vars.nJets);
	tree->SetBranchAddress("jetPt", &vars.jetPt);
	tree->SetBranchAddress("jetEta", &vars.jetEta);
	tree->SetBranchAddress("jetPhi", &vars.jetPhi);
	tree->SetBranchAddress("jetEnergy", &vars.jetEnergy);
	tree->SetBranchAddress("jetBTag", &vars.jetBTag);

	tree->SetBranchAddress("nTracks5", &vars.nTracks5);
	tree->SetBranchAddress("nTracks10", &vars.nTracks10);
	tree->SetBranchAddress("nTracks20", &vars.nTracks20);
	tree->SetBranchAddress("nTracks30", &vars.nTracks30);
	tree->SetBranchAddress("nTracks40", &vars.nTracks40);

	tree->SetBranchAddress("nGlobalMuons", &vars.nGlobalMuons);
	tree->SetBranchAddress("nStandaloneMuons", &vars.nStandaloneMuons);
	tree->SetBranchAddress("nPFMuons", &vars.nPFMuons);

	TH1F* h_nPV = dynamic_cast<TH1F*>(file->Get((directory + "/h_nPV").c_str()));
	out.add_events(static_cast<unsigned int>(h_nPV->GetEntries()));
	out.add_pv(h_nPV);

	if(!is_data)
	{
		TH1F* h_nTrueInteractions = dynamic_cast<TH1F*>(file->Get((directory + "/h_nTrueInteractions").c_str()));
		out.add_interactions(h_nTrueInteractions);
	}

	const unsigned int CountTotal = tree->GetEntries();
	for(int i = 0; i < CountTotal; ++i)
	{
		tree->GetEntry(i);

		// PDG filter
		if(!is_data && pdg_filter != 0 && abs(vars.NegMuon.GenChargedLepPDG) != pdg_filter && abs(vars.PosMuon.GenChargedLepPDG) != pdg_filter)
			continue;

		// Trigger
		if(!vars.hltMu17Mu8) continue;

		// Lumi Filter
		if(is_data && !LumiFilter::pass(vars.run, vars.lumi)) continue;

		// Muons
		Muon posMuon(vars.PosMuon);
		Muon negMuon(vars.NegMuon);

		// Trigger matching
		if(!(posMuon.HltMu17Mu8Leg8 && negMuon.HltMu17Mu8Leg17) && !(posMuon.HltMu17Mu8Leg17 && negMuon.HltMu17Mu8Leg8)) continue;

		// Jets
		std::vector<Jet> jets;
		for(unsigned int i = 0; i < vars.nJets; ++i)
			jets.push_back(Jet(vars, i));

		// PU weighting
		double weight = vars.weight;
		if(puWeightHisto) weight *= puWeightHisto->GetBinContent(puWeightHisto->FindBin(vars.nTrueInteractions));

		// Trigger weighting
		if(trigWeights)
		{
			if(posMuon.HltMu17Mu8Leg8 && negMuon.HltMu17Mu8Leg17)
				weight *= trigWeights->GetWeight(posMuon, negMuon);
			else
				weight *= trigWeights->GetWeight(negMuon, posMuon);
		}

		// ID/ISO Weighting
		if(idIsoWeights)
		{
			weight *= idIsoWeights->GetWeight(posMuon);
			weight *= idIsoWeights->GetWeight(negMuon);
		}

		// Output
		out.fill(posMuon, negMuon, jets, vars, weight);
	}

	delete file;
}

void process_dataset(const char* dataset, const std::vector<std::string>& files, bool is_ss, unsigned int pdg_filter)
{
	const bool is_data = (strstr(dataset, "Data") != NULL);

	std::stringstream out_filename_stream;
	out_filename_stream << "out-" << dataset;
	if(pdg_filter != 0) out_filename_stream << "-filter" << pdg_filter;
	if(is_ss) out_filename_stream << "-ss";
	out_filename_stream << ".root";
	const std::string out_filename = out_filename_stream.str();

	std::cout << dataset << "... => " << out_filename << std::endl;
	std::cout << "\tData: " << is_data << std::endl;

	TFile* outfile = new TFile(out_filename.c_str(), "RECREATE");
	Out out(outfile, is_data);

	TFile* puWeights = new TFile("../data/Pileup_Summer12_to_Run2012ABCD.root", "READ");
	TH1D* puWeightsHisto = NULL;
	if(puWeights->IsZombie())
		std::cout << "No PU reweighting file available: No PU reweighting performed!" << std::endl;
	else
		puWeightsHisto = dynamic_cast<TH1D*>(puWeights->Get("puWeight"));
	TriggerWeights trigWeights("../data/MuHLTEfficiencies_Run_2012ABCD_53X_DR03-2.root");
	IDIsoWeights idisoWeights("../data/MuonEfficiencies_Run2012ReReco_53X.root", "../data/MuonEfficiencies_ISO_Run_2012ReReco_53X.root");

	for(std::vector<std::string>::const_iterator iter = files.begin(); iter != files.end(); ++iter)
	{
		std::cout << "\t[" << (iter - files.begin() + 1) << "/" << files.size() << "] " << *iter << "..." << std::endl;
		process_file(iter->c_str(), out, is_data, is_ss, pdg_filter, is_data ? NULL : puWeightsHisto, is_data ? NULL : &trigWeights, is_data ? NULL : &idisoWeights);
	}

	outfile->Write();
	delete outfile;
}

int main(int argc, char* argv[])
{
	TH1::SetDefaultSumw2(true);

	bool is_ss = false;
	unsigned int pdg = 0;
	std::vector<std::string> files;
	for(int i = 1; i < argc; ++i)
	{
		char* endptr;
		unsigned long npdg = strtol(argv[i], &endptr, 10);
		if(*endptr == '\0')
		{
			pdg = npdg;
			std::cout << "PDG filter: " << pdg << std::endl;
		}
		else if(strcmp(argv[i], "SS") == 0)
		{
			is_ss = true;
		}
		else
		{
			files.push_back(argv[i]);
		}
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
		process_dataset(iter->first.c_str(), iter->second, is_ss, pdg);

	return 0;
}
