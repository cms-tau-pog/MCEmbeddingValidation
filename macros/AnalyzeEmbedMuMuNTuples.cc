#include <TTree.h>
#include <TFile.h>
#include <TMath.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <Math/VectorUtil.h>
#include <sstream>
#include <memory>

#include <TMVA/Reader.h>

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > RotLorentzVector;
typedef ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag> RotVector3;

RotLorentzVector rotate(const RotLorentzVector& p4, const RotLorentzVector& axis, double angle)
{
  TVector3 p3(p4.Px(), p4.Py(), p4.Pz());
  p3.Rotate(angle, TVector3(axis.x(), axis.y(), axis.z()).Unit());
  RotLorentzVector p4_rotated(p3.Px(), p3.Py(), p3.Pz(), p4.energy());
  assert(TMath::Abs(p3.Mag() - p4.P()) < (1.e-3*p4.P()));
  assert(TMath::Abs(p4.P() - p4_rotated.P()) < (1.e-3*p4.P()));
  assert(TMath::Abs(p4.M() - p4_rotated.M()) < (1.e-3*p4.P()));
  assert(TMath::Abs(ROOT::Math::VectorUtil::Angle(p4, axis) - ROOT::Math::VectorUtil::Angle(p4_rotated, axis)) < 1e-3);
  return p4_rotated;
}

std::pair<double, double> phiAngleInRF(const RotLorentzVector& first, const RotLorentzVector& second)
{
  RotLorentzVector pplus_lab(0., 0., 4000., sqrt(4000.*4000. + 0.939*0.939));
  RotLorentzVector pminus_lab(0., 0., -4000., sqrt(4000.*4000. + 0.939*0.939));

  RotLorentzVector muon1P4_lab = first;
  RotLorentzVector muon2P4_lab = second;
  RotLorentzVector zP4_lab = muon1P4_lab + muon2P4_lab;

  ROOT::Math::Boost boost_to_rf(zP4_lab.BoostToCM());
  RotLorentzVector muon1P4_rf = boost_to_rf(muon1P4_lab);
  RotLorentzVector muon2P4_rf = boost_to_rf(muon2P4_lab);
  RotLorentzVector pplus_rf = boost_to_rf(pplus_lab);
  RotLorentzVector pminus_rf = boost_to_rf(pminus_lab);

  // The proton momentum vector and the Z vector define a plane. Transform the proton momentum vector such that the two define the same plane, but the two vectors
  // are orthogonal to another. We do this by projecting the proton momentum vector on the Z momentum vector, and subtract the projected part from the initial vector.
  const RotVector3 pplus_rf_proj = pplus_rf.Vect() - pplus_rf.Vect().Dot(zP4_lab.Vect().Unit()) * zP4_lab.Vect().Unit();
  assert(fabs(pplus_rf_proj.Dot(zP4_lab.Vect())) < 1e-3);

  // Now, define a cartesian coordinate system such that the X axis is in the projected proton momentum direction, the Z axis is in the Z boson direction
  // and the Y axis is orthogonal on the two:
  const RotVector3 x = pplus_rf_proj.Unit();
  const RotVector3 z = zP4_lab.Vect().Unit();
  const RotVector3 y = x.Cross(z).Unit();

  // Get the x and y components which we use to define the phi angle:
  const double xcomp = muon1P4_rf.Vect().Dot(x);
  const double ycomp = muon1P4_rf.Vect().Dot(y);

  const double phi = atan2(ycomp,xcomp);
  return std::make_pair(phi, ROOT::Math::VectorUtil::Angle(muon1P4_rf, zP4_lab));
}

std::pair<double, double> phiAngleInRF(const TLorentzVector& first, const TLorentzVector& second)
{
  RotLorentzVector rFirst(first.Px(), first.Py(), first.Pz(), first.E());
  RotLorentzVector rSecond(second.Px(), second.Py(), second.Pz(), second.E());
  return phiAngleInRF(rFirst, rSecond);
}

std::pair<RotLorentzVector, RotLorentzVector> mirrorInRF(const RotLorentzVector& first, const RotLorentzVector& second)
{
  const RotLorentzVector pplus_lab(0., 0., 4000., sqrt(4000.*4000. + 0.939*0.939));
  const RotLorentzVector pminus_lab(0., 0., -4000., sqrt(4000.*4000. + 0.939*0.939));

  const RotLorentzVector muon1P4_lab = first;
  const RotLorentzVector muon2P4_lab = second;
  const RotLorentzVector zP4_lab = muon1P4_lab + muon2P4_lab;

  ROOT::Math::Boost boost_to_rf(zP4_lab.BoostToCM());
  const RotLorentzVector muon1P4_rf = boost_to_rf(muon1P4_lab);
  const RotLorentzVector muon2P4_rf = boost_to_rf(muon2P4_lab);
  const RotLorentzVector pplus_rf = boost_to_rf(pplus_lab);
  const RotLorentzVector pminus_rf = boost_to_rf(pminus_lab);

  // The proton momentum vector and the Z vector define a plane. Transform the proton momentum vector such that the two define the same plane, but the two vectors
  // are orthogonal to another. We do this by projecting the proton momentum vector on the Z momentum vector, and subtract the projected part from the initial vector.
  const RotVector3 pplus_rf_proj = pplus_rf.Vect() - pplus_rf.Vect().Dot(zP4_lab.Vect().Unit()) * zP4_lab.Vect().Unit();
  assert(fabs(pplus_rf_proj.Dot(zP4_lab.Vect())) < 1e-3);

  // Define a cartesian coordinate system such that the X axis is in the projected proton momentum direction, the Z axis is in the Z boson direction
  // and the Y axis is orthogonal on the two:
  const RotVector3 x = pplus_rf_proj.Unit();
  const RotVector3 z = zP4_lab.Vect().Unit();
  const RotVector3 y = x.Cross(z).Unit();

  // Mirror both vectors on the x-z plane
  const double ycomp1 = muon1P4_rf.Vect().Dot(y);
  const RotVector3 y1 = muon1P4_rf.Vect() - 2 * ycomp1 * y;
  assert(fabs(ycomp1 + y1.Dot(y)) < 1e-3);

  const double ycomp2 = muon2P4_rf.Vect().Dot(y);
  const RotVector3 y2 = muon2P4_rf.Vect() - 2 * ycomp2 * y;
  assert(fabs(ycomp2 + y2.Dot(y)) < 1e-3);

  // Boost back into rest-frame
  RotLorentzVector mirror1P4_rf, mirror2P4_rf;
  mirror1P4_rf.SetPxPyPzE(y1.x(), y1.y(), y1.z(), muon1P4_rf.E());
  mirror2P4_rf.SetPxPyPzE(y2.x(), y2.y(), y2.z(), muon2P4_rf.E());

  ROOT::Math::Boost boost_to_lab(boost_to_rf.Inverse());
  const RotLorentzVector mirror1P4_lab = boost_to_lab(mirror1P4_rf);
  const RotLorentzVector mirror2P4_lab = boost_to_lab(mirror2P4_rf);

  assert( fabs(zP4_lab.mass() - (mirror1P4_lab + mirror2P4_lab).mass()) <= (1.e-3*zP4_lab.mass()) &&
          fabs(zP4_lab.pt()   - (mirror1P4_lab + mirror2P4_lab).pt())   <= (1.e-3*zP4_lab.pt()));

  return std::make_pair(mirror1P4_lab, mirror2P4_lab);
}

std::pair<RotLorentzVector, RotLorentzVector> rotateInRF(const RotLorentzVector& first, const RotLorentzVector& second, double angle)
{
  RotLorentzVector muon1P4_lab = first;
  RotLorentzVector muon2P4_lab = second;
  RotLorentzVector zP4_lab = muon1P4_lab + muon2P4_lab;

  ROOT::Math::Boost boost_to_rf(zP4_lab.BoostToCM());
  ROOT::Math::Boost boost_to_lab(boost_to_rf.Inverse());

  //RotLorentzVector zP4_rf = boost_to_rf(zP4_lab);
  RotLorentzVector muon1P4_rf = boost_to_rf(muon1P4_lab);
  RotLorentzVector muon2P4_rf = boost_to_rf(muon2P4_lab);

  muon1P4_rf = rotate(muon1P4_rf, zP4_lab, angle * M_PI / 180.);
  muon2P4_rf = rotate(muon2P4_rf, zP4_lab, angle * M_PI / 180.);

  RotLorentzVector lep1P4_lab = boost_to_lab(muon1P4_rf);
  RotLorentzVector lep2P4_lab = boost_to_lab(muon2P4_rf);

  assert( fabs(zP4_lab.mass() - (lep1P4_lab + lep2P4_lab).mass()) <= (1.e-3*zP4_lab.mass()) &&
          fabs(zP4_lab.pt()   - (lep1P4_lab + lep2P4_lab).pt())   <= (1.e-3*zP4_lab.pt()));

  return std::make_pair(lep1P4_lab, lep2P4_lab);
}

std::pair<TLorentzVector, TLorentzVector> rotateInRF(const TLorentzVector& first, const TLorentzVector& second, double angle)
{
  RotLorentzVector rFirst(first.Px(), first.Py(), first.Pz(), first.E());
  RotLorentzVector rSecond(second.Px(), second.Py(), second.Pz(), second.E());
  std::pair<RotLorentzVector, RotLorentzVector> rot = rotateInRF(rFirst, rSecond, angle);
  TLorentzVector rrFirst, rrSecond;
  rrFirst.SetPxPyPzE(rot.first.px(), rot.first.py(), rot.first.pz(), rot.first.energy());
  rrSecond.SetPxPyPzE(rot.second.px(), rot.second.py(), rot.second.pz(), rot.second.energy());
  return std::make_pair(rrFirst, rrSecond);
}

std::pair<TLorentzVector, TLorentzVector> mirrorInRF(const TLorentzVector& first, const TLorentzVector& second)
{
  RotLorentzVector rFirst(first.Px(), first.Py(), first.Pz(), first.E());
  RotLorentzVector rSecond(second.Px(), second.Py(), second.Pz(), second.E());
  std::pair<RotLorentzVector, RotLorentzVector> rot = mirrorInRF(rFirst, rSecond);
  TLorentzVector rrFirst, rrSecond;
  rrFirst.SetPxPyPzE(rot.first.px(), rot.first.py(), rot.first.pz(), rot.first.energy());
  rrSecond.SetPxPyPzE(rot.second.px(), rot.second.py(), rot.second.pz(), rot.second.energy());
  return std::make_pair(rrFirst, rrSecond);
}

const bool lorentzEqual(const TLorentzVector& first, const TLorentzVector& second)
{
	return fabs(first.Px() - second.Px()) < 1e-3 && fabs(first.Py() - second.Py()) < 1e-3 && fabs(first.Pz() - second.Pz()) < 1e-3 && fabs(first.E() - second.E()) < 1e-3;
}

static TLorentzVector TLorentzVector_PtEtaPhiE(double pt, double eta, double phi, double e)
{
	TLorentzVector v;
	v.SetPtEtaPhiE(pt, eta, phi, e);
	return v;
}

struct TreeVars {
	TreeVars() {
		run = lumi = event = ~0u;
		muonRadiationFilter = muonRadiationFilter2Sel1 = muonRadiationFilter2Sel2 = muonRadiationFilter2Sel3 = true;
		nPV = 0; nTrueInteractions = -1.0f; nGenMEJets = 0;
		caloMet = caloMetPhi = pfMet = pfMetPhi = -1.0f;
		origCaloMet = origCaloMetPhi = origPfMet = origPfMetPhi = -1.0f;
		hltMu17Mu8 = origHltMu17Mu8 = false;

		PosMuon.Pt = PosMuon.Eta = PosMuon.Phi = PosMuon.E = -1.0f;
		PosMuon.GenPt = PosMuon.GenEta = PosMuon.GenPhi = PosMuon.GenE = -1.0f;
		PosMuon.GenFinalPt = PosMuon.GenFinalEta = PosMuon.GenFinalPhi = PosMuon.GenFinalE = -1.0f;
		PosMuon.OrigPt = PosMuon.OrigEta = PosMuon.OrigPhi = PosMuon.OrigE = -1.0f;
		PosMuon.OrigGenPt = PosMuon.OrigGenEta = PosMuon.OrigGenPhi = PosMuon.OrigGenE = -1.0f;
		PosMuon.OrigGenFinalPt = PosMuon.OrigGenFinalEta = PosMuon.OrigGenFinalPhi = PosMuon.OrigGenFinalE = -1.0f;
		PosMuon.HltMu17Mu8Leg8 = PosMuon.HltMu17Mu8Leg17 = false;
		PosMuon.ChargedHadronPtPfIso04 = PosMuon.ChargedParticlePtPfIso04 = PosMuon.NeutralHadronEtPfIso04 = PosMuon.PhotonEtPfIso04 = PosMuon.PUPtPfIso04 = -1.0f;
		PosMuon.Quality = false;
		PosMuon.GenChargedLepPDG = PosMuon.OrigGenChargedLepPDG = 0;

		NegMuon.Pt = NegMuon.Eta = NegMuon.Phi = NegMuon.E = -1.0f;
		NegMuon.GenPt = NegMuon.GenEta = NegMuon.GenPhi = NegMuon.GenE = -1.0f;
		NegMuon.GenFinalPt = NegMuon.GenFinalEta = NegMuon.GenFinalPhi = NegMuon.GenFinalE = -1.0f;
		NegMuon.OrigPt = NegMuon.OrigEta = NegMuon.OrigPhi = NegMuon.OrigE = -1.0f;
		NegMuon.OrigGenPt = NegMuon.OrigGenEta = NegMuon.OrigGenPhi = NegMuon.OrigGenE = -1.0f;
		NegMuon.OrigGenFinalPt = NegMuon.OrigGenFinalEta = NegMuon.OrigGenFinalPhi = NegMuon.OrigGenFinalE = -1.0f;
		NegMuon.HltMu17Mu8Leg8 = NegMuon.HltMu17Mu8Leg17 = false;
		NegMuon.ChargedHadronPtPfIso04 = NegMuon.ChargedParticlePtPfIso04 = NegMuon.NeutralHadronEtPfIso04 = NegMuon.PhotonEtPfIso04 = NegMuon.PUPtPfIso04 = -1.0f;
		NegMuon.Quality = false;
		NegMuon.GenChargedLepPDG = NegMuon.OrigGenChargedLepPDG = 0;

		nJets = 0;
		nTracks5 = nTracks10 = nTracks20 = nTracks30 = nTracks40 = 0;
	}

	unsigned int run;
	unsigned int lumi;
	unsigned int event;

	bool muonRadiationFilter;
	bool muonRadiationFilter2Sel1;
	bool muonRadiationFilter2Sel2;
	bool muonRadiationFilter2Sel3;

	unsigned int nPV;
	float nTrueInteractions;
	int nGenMEJets;

	bool hltMu17Mu8;
	bool origHltMu17Mu8;

	float caloMet;
	float caloMetPhi;
	float pfMet;
	float pfMetPhi;

	float origCaloMet;
	float origCaloMetPhi;
	float origPfMet;
	float origPfMetPhi;

	struct Muon {
		float Pt;
		float Eta;
		float Phi;
		float E;

		float GenPt;
		float GenEta;
		float GenPhi;
		float GenE;

		float GenFinalPt;
		float GenFinalEta;
		float GenFinalPhi;
		float GenFinalE;

		float OrigPt;
		float OrigEta;
		float OrigPhi;
		float OrigE;

		float OrigGenPt;
		float OrigGenEta;
		float OrigGenPhi;
		float OrigGenE;

		float OrigGenFinalPt;
		float OrigGenFinalEta;
		float OrigGenFinalPhi;
		float OrigGenFinalE;

		bool HltMu17Mu8Leg8;
		bool HltMu17Mu8Leg17;

		float ChargedHadronPtPfIso04;
		float ChargedParticlePtPfIso04;
		float NeutralHadronEtPfIso04;
		float PhotonEtPfIso04;
		float PUPtPfIso04;
		bool Quality;

		int GenChargedLepPDG;
		int OrigGenChargedLepPDG;
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
};

class Muon
{
public:
	Muon(const TreeVars::Muon& muon, int charge):
		hasGen(muon.GenPt > 0.0f),
		hasRec(muon.Pt > 0.0f),
		hasOrigGen(muon.OrigGenPt > 0.0f),
		hasOrigRec(muon.OrigPt > 0.0f),
		p4(TLorentzVector_PtEtaPhiE(muon.Pt, muon.Eta, muon.Phi, muon.E)),
		p4genInitial(TLorentzVector_PtEtaPhiE(muon.GenPt, muon.GenEta, muon.GenPhi, muon.GenE)),
		p4genFinal(TLorentzVector_PtEtaPhiE(muon.GenFinalPt, muon.GenFinalEta, muon.GenFinalPhi, muon.GenFinalE)),
		p4orig(TLorentzVector_PtEtaPhiE(muon.OrigPt, muon.OrigEta, muon.OrigPhi, muon.OrigE)),
		p4origGenInitial(TLorentzVector_PtEtaPhiE(muon.OrigGenPt, muon.OrigGenEta, muon.OrigGenPhi, muon.OrigGenE)),
		p4origGenFinal(TLorentzVector_PtEtaPhiE(muon.OrigGenFinalPt, muon.OrigGenFinalEta, muon.OrigGenFinalPhi, muon.OrigGenFinalE)),
		genCharge(muon.GenChargedLepPDG > 0 ? -1 : +1),
		charge(charge),
		quality(muon.Quality),
		chargedParticleIsolation(1.0f/muon.Pt * (muon.ChargedParticlePtPfIso04)),
		chargedHadronIsolation(1.0f/muon.Pt * (muon.ChargedHadronPtPfIso04)),
		neutralHadronIsolation(1.0f/muon.Pt * (muon.NeutralHadronEtPfIso04)),
		gammaIsolation(1.0f/muon.Pt * (muon.PhotonEtPfIso04)),
		puIsolation(1.0f/muon.Pt * (muon.PUPtPfIso04)),
		isolation(1.0f/muon.Pt * (muon.ChargedParticlePtPfIso04 + std::max(0.0f, muon.NeutralHadronEtPfIso04 + muon.PhotonEtPfIso04 - 0.5f * muon.PUPtPfIso04)))
	{
	}

	const bool hasGen;
	const bool hasRec;
	const bool hasOrigGen;
	const bool hasOrigRec;

	const TLorentzVector p4;
	const TLorentzVector p4genInitial;
	const TLorentzVector p4genFinal;
	const TLorentzVector p4orig;
	const TLorentzVector p4origGenInitial;
	const TLorentzVector p4origGenFinal;

	const int genCharge;
	const int charge;

	const bool quality;
	const float chargedParticleIsolation;
	const float chargedHadronIsolation;
	const float neutralHadronIsolation;
	const float gammaIsolation;
	const float puIsolation;
	const float isolation;
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

class ZmumuEvtSelEffCorrWeight
{
public:
	ZmumuEvtSelEffCorrWeight(const char* filename, const char* prefix)
	{
		muMinusPt_vs_muPlusPt = muMinusEta_vs_muPlusEta = NULL;
		file = new TFile(filename, "READ");
		if(file->IsZombie()) return;

		muMinusPt_vs_muPlusPt = dynamic_cast<TH2D*>(file->Get((std::string(prefix) + "ZmumuEvtSelEff_muMinusPt_vs_muPlusPt").c_str()));
		muMinusEta_vs_muPlusEta = dynamic_cast<TH2D*>(file->Get((std::string(prefix) + "ZmumuEvtSelEffCorr_muMinusEta_vs_muPlusEta").c_str()));
	}

	bool hasLevel1() const { return muMinusPt_vs_muPlusPt != NULL; }
	bool hasLevel2() const { return muMinusEta_vs_muPlusEta != NULL; }

	static int findBin(TAxis* axis, double x)
	{
		int bin = axis->FindFixBin(x);
		if(bin < 1) bin = 1;
		if(bin > axis->GetNbins()) bin = axis->GetNbins();
		return bin;
	}

	double GetPtWeight(const Muon& posMuon, const Muon& negMuon, bool orig) const
	{
		assert(hasLevel1());
		assert(!orig || (posMuon.hasOrigRec && negMuon.hasOrigRec));
		assert(orig || (posMuon.hasRec && negMuon.hasRec));

		const TLorentzVector& plus = orig ? posMuon.p4orig : posMuon.p4;
		const TLorentzVector& minus = orig ? negMuon.p4orig : negMuon.p4;

		const int binPtX = findBin(muMinusPt_vs_muPlusPt->GetXaxis(), plus.Pt());
		const int binPtY = findBin(muMinusPt_vs_muPlusPt->GetYaxis(), minus.Pt());
		double effPt = muMinusPt_vs_muPlusPt->GetBinContent(binPtX, binPtY);

		if(effPt < 1e-1) 1.; //std::cout << "ptplus=" << plus.Pt() << ", ptminus=" << minus.Pt() << std::endl;
		return 1./effPt;
	}

	double GetPtEtaWeight(const Muon& posMuon, const Muon& negMuon, bool orig) const
	{
		assert(hasLevel2());
		const double ptWeight = GetPtWeight(posMuon, negMuon, orig);
		assert(!orig || (posMuon.hasOrigRec && negMuon.hasOrigRec));
		assert(orig || (posMuon.hasRec && negMuon.hasRec));

		const TLorentzVector& plus = orig ? posMuon.p4orig : posMuon.p4;
		const TLorentzVector& minus = orig ? negMuon.p4orig : negMuon.p4;

		const int binEtaX = findBin(muMinusEta_vs_muPlusEta->GetXaxis(), plus.Eta());
		const int binEtaY = findBin(muMinusEta_vs_muPlusEta->GetYaxis(), minus.Eta());
		double effEta = muMinusEta_vs_muPlusEta->GetBinContent(binEtaX, binEtaY);

		if(effEta < 1e-1) return ptWeight; //std::cout << "etaplus=" << plus.Eta() << ", etaminus=" << minus.Eta() << std::endl;
		return ptWeight * 1./(effEta);
	}

private:
	TFile* file;
	TH2D* muMinusPt_vs_muPlusPt;
	TH2D* muMinusEta_vs_muPlusEta;
};

struct DiMuonAcceptance
{
	bool acceptZmumu(const TLorentzVector& first, const TLorentzVector& second)
	{
		// Embedding cuts (very loose), for Zmumu efficiency measurement
		if(first.Pt() < 8 || second.Pt() < 8) return false;
		else if(first.Pt() < 17 && second.Pt() < 17) return false;
		else if(fabs(first.Eta()) > 2.5) return false;
		else if(fabs(second.Eta()) > 2.5) return false;
		else if((first + second).M() < 20.0) return false;
		return true;
	}

	bool acceptEmb(const TLorentzVector& first, const TLorentzVector& second)
	{
		// Embedding cuts (very loose). Use 10/20 for pT cuts to be consistent with PF embedding.
		if(first.Pt() < 10 || second.Pt() < 10) return false;
		else if(first.Pt() < 20 && second.Pt() < 20) return false;
		else if(fabs(first.Eta()) > 2.5) return false;
		else if(fabs(second.Eta()) > 2.5) return false;
		else if((first + second).M() < 20.0) return false;
		return true;
	}

	bool acceptRec(const TLorentzVector& first, const TLorentzVector& second)
	{
		// Standard cuts.
		if(first.Pt() < 10 || second.Pt() < 10) return false;
		else if(first.Pt() < 20 && second.Pt() < 20) return false;
		else if(fabs(first.Eta()) > 2.1) return false;
		else if(fabs(second.Eta()) > 2.1) return false;
		else if((first + second).M() < 20.0) return false;
		return true;
	}

	DiMuonAcceptance(const Muon& first, const Muon& second, const TreeVars& vars):
		zmumuGen(first.hasGen && second.hasGen && acceptZmumu(first.p4genInitial, second.p4genInitial) && (first.p4genInitial + second.p4genInitial).M() > 50.),
		zmumuRec(first.hasRec && second.hasRec && vars.hltMu17Mu8 && acceptZmumu(first.p4, second.p4) && first.quality && second.quality && first.chargedHadronIsolation < 0.1 && second.chargedHadronIsolation < 0.1 && (first.p4 + second.p4).M() > 50.),
		zmumuGenEmb(first.hasGen && second.hasGen && acceptEmb(first.p4genInitial, second.p4genInitial) && (first.p4genInitial + second.p4genInitial).M() > 50. && fabs(first.p4genInitial.Eta()) < 2.4 && fabs(second.p4genInitial.Eta()) < 2.4),
		zmumuRecEmb(first.hasRec && second.hasRec && vars.hltMu17Mu8 && acceptEmb(first.p4, second.p4) && first.quality && second.quality && first.chargedHadronIsolation < 0.1 && second.chargedHadronIsolation < 0.1 && (first.p4 + second.p4).M() > 50. && fabs(first.p4.Eta()) < 2.4 && fabs(second.p4.Eta()) < 2.4),
		origGenAccepted(first.hasOrigGen && second.hasOrigGen && acceptEmb(first.p4origGenFinal, second.p4origGenFinal)),
		origAccepted(vars.origHltMu17Mu8 && first.hasOrigRec && second.hasOrigRec && acceptEmb(first.p4orig, second.p4orig)),
		genInitialAccepted(first.hasGen && second.hasGen && acceptEmb(first.p4genInitial, second.p4genInitial) && fabs(first.p4genInitial.Eta()) < 2.4 && fabs(second.p4genInitial.Eta()) < 2.4), // PF embedded and RH embedded have a huge different in generator level eta at eta=2.5, maybe it comes from the different muon ID? So for now, go only until 2.4
		genInitialMasscut( (first.p4genInitial + second.p4genInitial).M() > 50.),
		genFinalAccepted(first.hasGen && second.hasGen && acceptEmb(first.p4genFinal, second.p4genFinal) && fabs(first.p4genFinal.Eta()) < 2.4 && fabs(second.p4genFinal.Eta()) < 2.4), // PF embedded and RH embedded have a huge different in generator level eta at eta=2.5, maybe it comes from the different muon ID? So for now, go only until 2.4
		genFinalMasscut( (first.p4genFinal + second.p4genFinal).M() > 50.),
		recAccepted(first.hasRec && second.hasRec && acceptRec(first.p4, second.p4)),
		recIdentified(first.quality && second.quality),
		recIsolated(first.isolation < 0.1 && second.isolation < 0.1),
		recMasscut( (first.p4 + second.p4).M() > 60.),
		recLowmass( (first.p4 + second.p4).M() < 80.)
	{
	}

	const bool zmumuGen; // pt 17/8, eta 2.5
	const bool zmumuRec; // pt 17/8, eta 2.5
	const bool zmumuGenEmb; // pt 20/10, eta 2.4
	const bool zmumuRecEmb; // pt 20/10, eta 2.4

	const bool origGenAccepted;    // original generator muons within acceptance (pt 20/10, eta 2.4)
	const bool origAccepted;       // original reconstructed muons within acceptance (pt 20/10, eta 2.4)
	const bool genInitialAccepted; // pre-FSR generator muons within acceptance (pt 20/10, eta 2.4)
	const bool genInitialMasscut;  // pre-FSR generator muons > 50 GeV
	const bool genFinalAccepted;   // post-FSR generator muons within acceptance (pt 20/10, eta 2.4)
	const bool genFinalMasscut;    // post-FSR generator muons > 50 GeV
	const bool recAccepted;        // reconstructed muons within acceptance (pt 20/10, eta 2.1)
	const bool recIdentified;      // reconstructed muons survive muon ID criteria
	const bool recIsolated;        // reconstructed muons survive relative combined PF isolation
	const bool recMasscut;         // reconstructed muons > 60 GeV
	const bool recLowmass;         // reconstructed muons < 80 GeV
};

// Histograms for studying the FSR filter: profile plots of filter results vs. FSR on generator level
struct RadiationHistograms {
	RadiationHistograms(TDirectory* dir, const char* name)
	{
		TDirectory* newdir = dir->mkdir(name);
		newdir->cd();

		const double relBinning[] = { 0., 0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.004, 0.006, 0.01, 0.015, 0.02, 0.03, 0.05, 0.1, 1.0 };
		const double absBinning[] = { 0., 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 100.0 };

		const int relNBins = sizeof(relBinning)/sizeof(relBinning[0]) - 1;
		const int absNBins = sizeof(absBinning)/sizeof(absBinning[0]) - 1;

		absFilter = new TProfile("absFilter", "Filter Result;E_{loss}^{rad} [GeV];Fiter Result", absNBins, absBinning);
		absFilter2Sel1 = new TProfile("absFilter2Sel1", "Filter Result;E_{loss}^{rad} [GeV];Fiter Result", absNBins, absBinning);
		absFilter2Sel2 = new TProfile("absFilter2Sel2", "Filter Result;E_{loss}^{rad} [GeV];Fiter Result", absNBins, absBinning);
		absFilter2Sel3 = new TProfile("absFilter2Sel3", "Filter Result;E_{loss}^{rad} [GeV];Fiter Result", absNBins, absBinning);

		dir->cd();
	}

	void fill(const Muon& posMuon, const Muon& negMuon, const TreeVars& vars, float weight)
	{
		assert(posMuon.hasOrigGen && posMuon.hasOrigGen);

		const double Einitial = posMuon.p4origGenInitial.E() + negMuon.p4origGenInitial.E();
		const double Efinal = posMuon.p4origGenFinal.E() + negMuon.p4origGenFinal.E();

		const double Eloss = Einitial - Efinal;
		assert(Eloss >= -1e-3);
		assert(Eloss <= Einitial+1e-3);
		const double ElossRel = Eloss / Einitial;

		absFilter->Fill(Eloss, vars.muonRadiationFilter ? 1.0 : 0.0, weight);
		absFilter2Sel1->Fill(Eloss, vars.muonRadiationFilter2Sel1 ? 1.0 : 0.0, weight);
		absFilter2Sel2->Fill(Eloss, vars.muonRadiationFilter2Sel2 ? 1.0 : 0.0, weight);
		absFilter2Sel3->Fill(Eloss, vars.muonRadiationFilter2Sel3 ? 1.0 : 0.0, weight);
	}

	TProfile* absFilter;
	TProfile* absFilter2Sel1;
	TProfile* absFilter2Sel2;
	TProfile* absFilter2Sel3;
};

// Histograms for creating the Zmumu event selection efficiency corrections
struct ZmumuEffHistograms {
	ZmumuEffHistograms(TDirectory* dir, const char* name)
	{
		TDirectory* newdir = dir->mkdir(name);
		newdir->cd();

		double pt_binning[] = { 0., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 24., 25., 26., 27., 28., 29., 30., 31., 32., 33., 34., 35., 36., 37., 38., 39., 40., 42., 44., 46., 48., 50., 55., 60., 70., 80., 100., 150., 250. };
		unsigned int nptbins = sizeof(pt_binning)/sizeof(pt_binning[0]) - 1;

		genPt = new TH2D("genPt", "genPt;posMuon p_{T};negMuon p_{T}", nptbins, pt_binning, nptbins, pt_binning);
		genEta = new TH2D("genEta", "genEta;posMuon #eta;negMuon #eta", 50, -2.5, 2.5, 50, -2.5, 2.5);
		recPt = new TH2D("recPt", "recPt;posMuon p_{T};negMuon p_{T}", nptbins, pt_binning, nptbins, pt_binning);
		recEta = new TH2D("recEta", "recEta;posMuon #eta;negMuon #eta", 50, -2.5, 2.5, 50, -2.5, 2.5);

		dir->cd();
	}

	void fill(const Muon& posMuon, const Muon& negMuon, float weight)
	{
		if(posMuon.hasGen && negMuon.hasGen)
		{
			genPt->Fill(posMuon.p4genInitial.Pt(), negMuon.p4genInitial.Pt(), weight);
			genEta->Fill(posMuon.p4genInitial.Eta(), negMuon.p4genInitial.Eta(), weight);
		}

		if(posMuon.hasRec && negMuon.hasRec)
		{
			recPt->Fill(posMuon.p4.Pt(), negMuon.p4.Pt(), weight);
			recEta->Fill(posMuon.p4.Eta(), negMuon.p4.Eta(), weight);
		}
	}

	TH2D* genPt;
	TH2D* genEta;
	TH2D* recPt;
	TH2D* recEta;
};

// Histograms which show what happens when applying a four-vector transformation on
// the input muons.
struct TransformHistograms
{
	TransformHistograms(TDirectory* dir, const char* name)
	{
		TDirectory* newdir = dir->mkdir(name);
		newdir->cd();

		diMuonMass = new TH1D("diMuonMass", "diMuonMass;M_{#mu#mu};Entries", 200, 0., 200.);
		diMuonPt = new TH1D("diMuonPt", "diMuonPt;p_{T};Entries", 200, 0., 200.);

		posMuonPt = new TH1D("posMuonPt", "posMuonPt;p_{T};Entries", 200, 0., 200.);
		posMuonEta = new TH1D("posMuonEta", "posMuonEta;#eta;Entries", 180, -9., 9.);
		posMuonPhi = new TH1D("posMuonPhi", "posMuonPhi;#phi;Entries", 100, -M_PI, M_PI);
		negMuonPt = new TH1D("negMuonPt", "negMuonPt;p_{T};Entries", 200, 0., 200.);
		negMuonEta = new TH1D("negMuonEta", "negMuonEta;#eta;Entries", 180, -9., 9.);
		negMuonPhi = new TH1D("negMuonPhi", "negMuonPhi;#phi;Entries", 100, -M_PI, M_PI);

		posMuonRotPhi = new TH1D("posMuonRotPhi", "posMuonRotPhi;#phi;Entries", 100, -M_PI, M_PI);
		posMuonRotTheta = new TH1D("posMuonRotTheta", "posMuonRotTheta;#phi;Entries", 100, 0, M_PI);
		deltaR = new TH1D("deltaR", "deltaR;#Delta R;Entries", 100, 0., 5.);

		dir->cd();
	}

	void fill(const TLorentzVector& posMuon, const TLorentzVector& negMuon,
	          const TLorentzVector& origPosMuon, const TLorentzVector& origNegMuon,
		  float weight)
	{
		const TLorentzVector diMuon = posMuon + negMuon;

		diMuonMass->Fill(diMuon.M(), weight);
		diMuonPt->Fill(diMuon.Pt(), weight);

		posMuonPt->Fill(posMuon.Pt(), weight);
		posMuonEta->Fill(posMuon.Eta(), weight);
		posMuonPhi->Fill(posMuon.Phi(), weight);
		negMuonPt->Fill(negMuon.Pt(), weight);
		negMuonEta->Fill(negMuon.Eta(), weight);
		negMuonPhi->Fill(negMuon.Phi(), weight);

                const std::pair<double, double> rot = phiAngleInRF(posMuon, negMuon);
		posMuonRotPhi->Fill(rot.first, weight);
		posMuonRotTheta->Fill(rot.second, weight);

		const double posMuonDeltaR = std::min(ROOT::Math::VectorUtil::DeltaR(posMuon, origPosMuon), ROOT::Math::VectorUtil::DeltaR(posMuon, origNegMuon));
		const double negMuonDeltaR = std::min(ROOT::Math::VectorUtil::DeltaR(negMuon, origPosMuon), ROOT::Math::VectorUtil::DeltaR(negMuon, origNegMuon));
		deltaR->Fill(posMuonDeltaR, weight);
		deltaR->Fill(negMuonDeltaR, weight);
	}

	TH1D* diMuonMass;
	TH1D* diMuonPt;

	TH1D* posMuonPt;
	TH1D* posMuonEta;
	TH1D* posMuonPhi;

	TH1D* negMuonPt;
	TH1D* negMuonEta;
	TH1D* negMuonPhi;

	TH1D* posMuonRotPhi;
	TH1D* posMuonRotTheta;

	TH1D* deltaR;
};

struct Histograms {
	Histograms(TDirectory* dir, const char* name, bool is_data, bool is_embedded, bool is_gen_embedded) {
		TDirectory* newdir = dir->mkdir(name);
		newdir->cd();

		nPV = new TH1D("nPV", "Number of primary vertices;N_{PV};Entries", 50, -0.5, 49.5);
		nPV_u = new TH1D("nPV_u", "Number of primary vertices, unweighted;N_{PV};Entries", 50, -0.5, 49.5);

		if(!is_data)
		{
			posGenMuonPt = new TH1D("posGenMuonPt", "Positive generated muon transverse momentum;p_{T};Entries/GeV", 200, 0.0, 200.0);
			posGenMuonEta = new TH1D("posGenMuonEta", "Positive generated muon pseudo-rapidity;#eta;Entries/0.1", 180, -9.0, 9.0);
			posGenMuonPhi = new TH1D("posGenMuonPhi", "Positive generated muon azimuthal angle;#phi;Entries", 100, -M_PI, M_PI);

			negGenMuonPt = new TH1D("negGenMuonPt", "Positive generated muon transverse momentum;p_{T};Entries/GeV", 200, 0.0, 200.0);
			negGenMuonEta = new TH1D("negGenMuonEta", "Positive generated muon pseudo-rapidity;#eta;Entries/0.1", 180, -9.0, 9.0);
			negGenMuonPhi = new TH1D("negGenMuonPhi", "Positive generated muon azimuthal angle;#phi;Entries", 100, -M_PI, M_PI);

			genDiMuonPt = new TH1D("genDiMuonPt", "Generated Di-muon transverse momentum;p_{T};Entries/GeV", 200, 0.0, 200.0);
			genDiMuonEta = new TH1D("genDiMuonEta", "Generated Di-muon pseudo-rapidity;#eta;Entries/0.1", 180, -9.0, 9.0);
			genInitialDiMuonMass = new TH1D("genInitialDiMuonMass", "Generated Di-muon invariant mass;m_{#mu#mu};Entries/GeV", 400, 0.0, 200.0);
			genFinalDiMuonMass = new TH1D("genFinalDiMuonMass", "Generated Di-muon invariant mass;m_{#mu#mu};Entries/GeV", 400, 0.0, 200.0);

			genDiMuonDeltaPt = new TH1D("genDiMuonDeltaPt", "Generated DeltaPt of the two muons;#Delta p_{T};Entries/GeV", 200, 0., 200.);
			posGenDeltaRZ = new TH1D("posGenDeltaRZ", "Generated DeltaR between positive muon and Z;#Delta R;Entries", 100, 0., 5.);
		}

		posMuonPt = new TH1D("posMuonPt", "Positive muon transverse momentum;p_{T};Entries/GeV", 200, 0.0, 200.0);
		posMuonEta = new TH1D("posMuonEta", "Positive muon pseudo-rapidity;#eta;Entries/0.1", 180, -9.0, 9.0);
		posMuonPhi = new TH1D("posMuonPhi", "Positive muon azimuthal angle;#phi;Entries", 100, -M_PI, M_PI);
		posMuonIso = new TH1D("posMuonIso", "Positive muon relative PF isolation;I_{rel}^{PF}", 100, 0.0, 1.0);

		posMuonChargedParticleIso = new TH1D("posMuonChargedParticleIso", "Positive muon relative charged particle isolation;I_{rel}^{charged}", 100, 0.0, 1.0);
		posMuonChargedHadronIso = new TH1D("posMuonChargedHadronIso", "Positive muon relative charged hadron isolation;I_{rel}^{charged}", 100, 0.0, 1.0);
		posMuonNeutralHadronIso = new TH1D("posMuonNeutralHadronIso", "Positive muon relative neutral hadron isolation;I_{rel}^{charged}", 100, 0.0, 1.0);
		posMuonGammaIso = new TH1D("posMuonGammaIso", "Positive muon relative gamma isolation;I_{rel}^{charged}", 100, 0.0, 1.0);
		posMuonPUIso = new TH1D("posMuonPUIso", "Positive muon relative PU isolation;I_{rel}^{charged}", 100, 0.0, 1.0);

		negMuonPt = new TH1D("negMuonPt", "Positive muon transverse momentum;p_{T};Entries/GeV", 200, 0.0, 200.0);
		negMuonEta = new TH1D("negMuonEta", "Positive muon pseudo-rapidity;#eta;Entries/0.1", 180, -9.0, 9.0);
		negMuonPhi = new TH1D("negMuonPhi", "Positive muon azimuthal angle;#phi;Entries", 100, -M_PI, M_PI);
		negMuonIso = new TH1D("negMuonIso", "Positive muon relative PF isolation;I_{rel}^{PF}", 100, 0.0, 1.0);

		leadMuonPt = new TH1D("leadMuonPt", "Leading muon transverse momentum;p_{T};Entries/GeV", 200, 0.0, 200.0);
		leadMuonEta = new TH1D("leadMuonEta", "Leading muon pseudo-rapidity;#eta;Entries/0.1", 180, -9.0, 9.0);
		leadMuonPhi = new TH1D("leadMuonPhi", "Leading muon azimuthal angle;#phi;Entries", 100, -M_PI, M_PI);
		leadMuonIso = new TH1D("leadMuonIso", "Leading muon relative PF isolation;I_{rel}^{PF}", 100, 0.0, 1.0);

		subMuonPt = new TH1D("subMuonPt", "Subleading muon transverse momentum;p_{T};Entries/GeV", 200, 0.0, 200.0);
		subMuonEta = new TH1D("subMuonEta", "Subleading muon pseudo-rapidity;#eta;Entries/0.1", 180, -9.0, 9.0);
		subMuonPhi = new TH1D("subMuonPhi", "Subleading muon azimuthal angle;#phi;Entries", 100, -M_PI, M_PI);
		subMuonIso = new TH1D("subMuonIso", "Subleading muon relative PF isolation;I_{rel}^{PF}", 100, 0.0, 1.0);

		diMuonPt = new TH1D("diMuonPt", "Di-muon transverse momentum;p_{T};Entries/GeV", 200, 0.0, 200.0);
		diMuonEta = new TH1D("diMuonEta", "Di-muon pseudo-rapidity;#eta;Entries/0.1", 180, -9.0, 9.0);
		diMuonPhi = new TH1D("diMuonPhi", "Di-muon azimuthal angle;#phi;Entries", 100, -M_PI, M_PI);
		diMuonDeltaPhi = new TH1D("diMuonDeltaPhi", "Di-muon #Delta #phi;#Delta Phi;Entries", 100, 0, M_PI);
		diMuonDeltaR = new TH1D("diMuonDeltaR", "Di-muon #Delta R;#Delta R;Entries/0.05", 200, 0.0, 10.0);
		diMuonMass = new TH1D("diMuonMass", "Di-muon invariant mass;m_{#mu#mu};Entries/GeV", 400, 0.0, 200.0);

		diMuonDeltaPt = new TH1D("diMuonDeltaPt", "DeltaPt of the two muons;#Delta p_{T};Entries/GeV", 200, 0., 200.);
		posDeltaRZ = new TH1D("posDeltaRZ", "DeltaR between positive muon and Z;#Delta R;Entries", 100, 0., 5.);

		nJets30 = new TH1D("nJets30", "Number of jets with p_{T} > 30 GeV;N_{jets};Entries", 20, -0.5, 19.5);
		leadJetPt = new TH1D("leadJetPt", "Leading jet transverse momentum;p_{T} [GeV];Entries/GeV", 270, 30., 300.);

		nTracks5 = new TH1D("nTracks5", "Number of tracks with p_{T} > 5 GeV;N_{tracks};Entries", 100, -0.5, 99.5);
		nTracks10 = new TH1D("nTracks10", "Number of tracks with p_{T} > 10 GeV;N_{tracks};Entries", 50, -0.5, 49.5);
		nTracks20 = new TH1D("nTracks20", "Number of tracks with p_{T} > 20 GeV;N_{tracks};Entries", 50, -0.5, 49.5);
		nTracks30 = new TH1D("nTracks30", "Number of tracks with p_{T} > 30 GeV;N_{tracks};Entries", 40, -0.5, 39.5);
		nTracks40 = new TH1D("nTracks40", "Number of tracks with p_{T} > 40 GeV;N_{tracks};Entries", 20, -0.5, 19.5);

		caloMet = new TH1D("caloMet", "Calorimeter missing transverse energy;E_{T}^{miss};Entries/GeV", 200, 0.0, 200.0);
		pfMet = new TH1D("pfMet", "Particle Flow missing transverse energy;E_{T}^{miss};Entries/GeV", 200, 0.0, 200.0);

		deltaCaloMetOrigCaloMet = deltaPfMetOrigPfMet = NULL;
		if(is_embedded) deltaCaloMetOrigCaloMet = new TH1D("deltaCaloMetOrigCaloMet", "Delta CaloMet;#DeltaE_{T}^{miss};Entries/0.1 GeV", 200, 0.0, 20.0);
		if(is_embedded) deltaPfMetOrigPfMet = new TH1D("deltaPfMetOrigPfMet", "Delta PfMet;#DeltaE_{T}^{miss};Entries/0.1 GeV", 200, 0.0, 20.0);

		dir->cd();
	}

	void fill(const Muon& posMuon, const Muon& negMuon, const std::vector<Jet>& jets30, const TreeVars& vars, float weight)
	{
		nPV->Fill(vars.nPV, weight);
		nPV_u->Fill(vars.nPV);

		if(posMuon.hasGen && negMuon.hasGen)
		{
			posGenMuonPt->Fill(posMuon.p4genFinal.Pt(), weight);
			posGenMuonEta->Fill(posMuon.p4genFinal.Eta(), weight);
			posGenMuonPhi->Fill(posMuon.p4genFinal.Phi(), weight);

			negGenMuonPt->Fill(negMuon.p4genFinal.Pt(), weight);
			negGenMuonEta->Fill(negMuon.p4genFinal.Eta(), weight);
			negGenMuonPhi->Fill(negMuon.p4genFinal.Phi(), weight);

			const TLorentzVector diMuonInitial = posMuon.p4genInitial + negMuon.p4genInitial;
			const TLorentzVector diMuonFinal = posMuon.p4genFinal + negMuon.p4genFinal;
			genDiMuonPt->Fill(diMuonFinal.Pt(), weight);
			genDiMuonEta->Fill(diMuonFinal.Eta(), weight);
			genInitialDiMuonMass->Fill(diMuonInitial.M(), weight);
			genFinalDiMuonMass->Fill(diMuonFinal.M(), weight);

			genDiMuonDeltaPt->Fill(fabs(posMuon.p4genFinal.Pt() - negMuon.p4genFinal.Pt()), weight);
			posGenDeltaRZ->Fill(ROOT::Math::VectorUtil::DeltaR(posMuon.p4genFinal, diMuonFinal), weight);
		}

		if(posMuon.hasRec && negMuon.hasRec)
		{
			const Muon& leadMuon = posMuon.p4.Pt() > negMuon.p4.Pt() ? posMuon : negMuon;
			const Muon& subMuon = posMuon.p4.Pt() <= negMuon.p4.Pt() ? posMuon : negMuon;

			posMuonPt->Fill(posMuon.p4.Pt(), weight);
			posMuonEta->Fill(posMuon.p4.Eta(), weight);
			posMuonPhi->Fill(posMuon.p4.Phi(), weight);
			posMuonIso->Fill(posMuon.isolation, weight);

			posMuonChargedParticleIso->Fill(posMuon.chargedParticleIsolation);
			posMuonChargedHadronIso->Fill(posMuon.chargedHadronIsolation);
			posMuonNeutralHadronIso->Fill(posMuon.neutralHadronIsolation);
			posMuonGammaIso->Fill(posMuon.gammaIsolation);
			posMuonPUIso->Fill(posMuon.puIsolation);

			negMuonPt->Fill(negMuon.p4.Pt(), weight);
			negMuonEta->Fill(negMuon.p4.Eta(), weight);
			negMuonPhi->Fill(negMuon.p4.Phi(), weight);
			negMuonIso->Fill(negMuon.isolation, weight);

			leadMuonPt->Fill(leadMuon.p4.Pt(), weight);
			leadMuonEta->Fill(leadMuon.p4.Eta(), weight);
			leadMuonPhi->Fill(leadMuon.p4.Phi(), weight);
			leadMuonIso->Fill(leadMuon.isolation, weight);

			subMuonPt->Fill(subMuon.p4.Pt(), weight);
			subMuonEta->Fill(subMuon.p4.Eta(), weight);
			subMuonPhi->Fill(subMuon.p4.Phi(), weight);
			subMuonIso->Fill(subMuon.isolation, weight);

			diMuonPt->Fill( (posMuon.p4 + negMuon.p4).Pt(), weight);
			diMuonEta->Fill( (posMuon.p4 + negMuon.p4).Eta(), weight);
			diMuonPhi->Fill( (posMuon.p4 + negMuon.p4).Phi(), weight);
			diMuonDeltaPhi->Fill(ROOT::Math::VectorUtil::DeltaPhi(posMuon.p4, negMuon.p4), weight);
			diMuonDeltaR->Fill(ROOT::Math::VectorUtil::DeltaR(posMuon.p4, negMuon.p4), weight);
			diMuonMass->Fill( (posMuon.p4 + negMuon.p4).M(), weight);

			diMuonDeltaPt->Fill(fabs(posMuon.p4.Pt() - negMuon.p4.Pt()), weight);
			posDeltaRZ->Fill(ROOT::Math::VectorUtil::DeltaR(posMuon.p4, (posMuon.p4 + negMuon.p4)), weight);
		}

		nJets30->Fill(jets30.size(), weight);
		if(!jets30.empty())
			leadJetPt->Fill(jets30[0].p4.Pt(), weight);

		nTracks5->Fill(vars.nTracks5, weight);
		nTracks10->Fill(vars.nTracks10, weight);
		nTracks20->Fill(vars.nTracks20, weight);
		nTracks30->Fill(vars.nTracks30, weight);
		nTracks40->Fill(vars.nTracks40, weight);

		TVector2 caloMetVec, pfMetVec;
		caloMetVec.SetMagPhi(vars.caloMet, vars.caloMetPhi);
		pfMetVec.SetMagPhi(vars.pfMet, vars.pfMetPhi);
		caloMet->Fill(caloMetVec.Mod(), weight);
		pfMet->Fill(pfMetVec.Mod(), weight);

		if(posMuon.hasOrigRec && posMuon.hasRec && negMuon.hasOrigRec && negMuon.hasRec)
		{
			TVector2 origCaloMetVec, origPfMetVec;
			origCaloMetVec.SetMagPhi(vars.origCaloMet, vars.origCaloMetPhi);
			origPfMetVec.SetMagPhi(vars.origPfMet, vars.origPfMetPhi);

			if(deltaCaloMetOrigCaloMet) deltaCaloMetOrigCaloMet->Fill( (caloMetVec - origCaloMetVec).Mod(), weight);
			if(deltaPfMetOrigPfMet) deltaPfMetOrigPfMet->Fill( (pfMetVec - origPfMetVec).Mod(), weight);
		}
	}

	TH1D* nPV;
	TH1D* nPV_u;

	TH1D* posMuonFSR;
	TH1D* negMuonFSR;
	TH1D* diMuonFSR;

	TH1D* posGenMuonPt;
	TH1D* posGenMuonEta;
	TH1D* posGenMuonPhi;

	TH1D* negGenMuonPt;
	TH1D* negGenMuonEta;
	TH1D* negGenMuonPhi;

	TH1D* genDiMuonPt;
	TH1D* genDiMuonEta;
	TH1D* genInitialDiMuonMass;
	TH1D* genFinalDiMuonMass;

	TH1D* genDiMuonDeltaPt;
	TH1D* posGenDeltaRZ;

	TH1D* posMuonPt;
	TH1D* posMuonEta;
	TH1D* posMuonPhi;
	TH1D* posMuonIso;

	TH1D* posMuonChargedParticleIso;
	TH1D* posMuonChargedHadronIso;
	TH1D* posMuonNeutralHadronIso;
	TH1D* posMuonGammaIso;
	TH1D* posMuonPUIso;

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

	TH1D* diMuonDeltaPt;
	TH1D* posDeltaRZ;

	TH1D* nJets30;
	TH1D* leadJetPt;

	TH1D* nTracks5;
	TH1D* nTracks10;
	TH1D* nTracks20;
	TH1D* nTracks30;
	TH1D* nTracks40;

	TH1D* caloMet;
	TH1D* pfMet;

	TH1D* deltaCaloMetOrigCaloMet;
	TH1D* deltaPfMetOrigPfMet;
};

// These are the definitions of various four-vector transformations. These should only be studied
// with samples produced with the "writeAllEvents" flag set to TRUE in the mumu_ntuple_producer_cfg.py,
// otherwise events are lost due to acceptance effects.
struct TransformCategory
{
	TransformCategory(TDirectory* dir, const char* name)
	{
		TDirectory* newdir = dir->mkdir(name);
		newdir->cd();

		h_id.reset(new TransformHistograms(newdir, "id"));
		h_rot0.reset(new TransformHistograms(newdir, "rot0"));
		h_rot45.reset(new TransformHistograms(newdir, "rot45"));
		h_rotM45.reset(new TransformHistograms(newdir, "rotM45"));
		h_rot90.reset(new TransformHistograms(newdir, "rot90"));
		h_rotM90.reset(new TransformHistograms(newdir, "rotM90"));
		h_rot180.reset(new TransformHistograms(newdir, "rot180"));
		h_mirror.reset(new TransformHistograms(newdir, "mirror"));

		dir->cd();
	}

	void fill(const TLorentzVector& posMuon, const TLorentzVector& negMuon, float weight)
	{
		const std::pair<TLorentzVector, TLorentzVector> rot0 = rotateInRF(posMuon, negMuon, 0.);
		const std::pair<TLorentzVector, TLorentzVector> rot45 = rotateInRF(posMuon, negMuon, 45.);
		const std::pair<TLorentzVector, TLorentzVector> rot90 = rotateInRF(posMuon, negMuon, 90.);
		const std::pair<TLorentzVector, TLorentzVector> rot180 = rotateInRF(posMuon, negMuon, 180.);
		const std::pair<TLorentzVector, TLorentzVector> rot45neg = rotateInRF(posMuon, negMuon, -45.);
		const std::pair<TLorentzVector, TLorentzVector> rot90neg = rotateInRF(posMuon, negMuon, -90.);
		const std::pair<TLorentzVector, TLorentzVector> rot180neg = rotateInRF(posMuon, negMuon, -180.);
		const std::pair<TLorentzVector, TLorentzVector> mirror = mirrorInRF(posMuon, negMuon);

		h_id->fill(posMuon, negMuon, posMuon, negMuon, weight);
		h_rot0->fill(rot0.first, rot0.second, posMuon, negMuon, weight);
		h_rot45->fill(rot45.first, rot45.second, posMuon, negMuon, weight);
		h_rotM45->fill(rot45neg.first, rot45neg.second, posMuon, negMuon, weight);
		h_rot90->fill(rot90.first, rot90.second, posMuon, negMuon, weight);
		h_rotM90->fill(rot90neg.first, rot90neg.second, posMuon, negMuon, weight);
		h_rot180->fill(rot180.first, rot180.second, posMuon, negMuon, weight);
		h_mirror->fill(mirror.first, mirror.second, posMuon, negMuon, weight);
	}

	std::auto_ptr<TransformHistograms> h_id;
	std::auto_ptr<TransformHistograms> h_rot0;
	std::auto_ptr<TransformHistograms> h_rot45;
	std::auto_ptr<TransformHistograms> h_rotM45;
	std::auto_ptr<TransformHistograms> h_rot90;
	std::auto_ptr<TransformHistograms> h_rotM90;
	std::auto_ptr<TransformHistograms> h_rot180;
	std::auto_ptr<TransformHistograms> h_mirror;
};

struct Category
{
	Category(TDirectory* dir, const char* name, bool is_data, bool is_embedded, bool is_gen_embedded, bool allowOrigAccepted):
		isEmbedded(is_embedded),
		isGenEmbedded(is_gen_embedded)
	{
		TDirectory* newdir = dir->mkdir(name);
		newdir->cd();

		if( (!is_data && !is_embedded) || (is_embedded && is_gen_embedded)) genFinalAccepted.reset(new Histograms(newdir, "genFinalAccepted", is_data, is_embedded, is_gen_embedded));

		if(is_embedded && allowOrigAccepted) origAccepted.reset(new Histograms(newdir, "origAccepted", is_data, is_embedded, is_gen_embedded));
		if(is_embedded) origAndGenAccepted.reset(new Histograms(newdir, "origAndGenAccepted", is_data, is_embedded, is_gen_embedded));
		genInitialMasscut.reset(new Histograms(newdir, "genInitialMasscut", is_data, is_embedded, is_gen_embedded));
		genFinalMasscut.reset(new Histograms(newdir, "genFinalMasscut", is_data, is_embedded, is_gen_embedded));
		recAccepted.reset(new Histograms(newdir, "recAccepted", is_data, is_embedded, is_gen_embedded));
		recMasscut.reset(new Histograms(newdir, "recMasscut", is_data, is_embedded, is_gen_embedded));
		recLowmass.reset(new Histograms(newdir, "recLowmass", is_data, is_embedded, is_gen_embedded));
		recMasscutNoIsolation.reset(new Histograms(newdir, "recMasscutNoIsolation", is_data, is_embedded, is_gen_embedded));

		dir->cd();
	}

	void fill(const DiMuonAcceptance& acceptance, const Muon& muon1, const Muon& muon2, const std::vector<Jet>& jets30, const TreeVars& vars, float weight)
	{
		const bool isOrigAccepted = (!isEmbedded || (isEmbedded && !isGenEmbedded && acceptance.origAccepted) || (isEmbedded && isGenEmbedded && acceptance.origGenAccepted));
		const bool recTriggered = (isEmbedded || (!isEmbedded && vars.hltMu17Mu8));

		if(acceptance.genFinalAccepted) if(genFinalAccepted.get()) genFinalAccepted->fill(muon1, muon2, jets30, vars, weight);
		if(isOrigAccepted) if(origAccepted.get()) origAccepted->fill(muon1, muon2, jets30, vars, weight);
		if(isOrigAccepted && acceptance.genInitialAccepted) if(origAndGenAccepted.get()) origAndGenAccepted->fill(muon1, muon2, jets30, vars, weight);
		if(isOrigAccepted && acceptance.genInitialAccepted && acceptance.genInitialMasscut) if(genInitialMasscut.get()) genInitialMasscut->fill(muon1, muon2, jets30, vars, weight);
		if(isOrigAccepted && acceptance.genFinalAccepted && acceptance.genFinalMasscut) if(genFinalMasscut.get()) genFinalMasscut->fill(muon1, muon2, jets30, vars, weight);
		if(isOrigAccepted && acceptance.genFinalAccepted && acceptance.genFinalMasscut && acceptance.recAccepted && recTriggered) recAccepted->fill(muon1, muon2, jets30, vars, weight);
		if(isOrigAccepted && acceptance.genFinalAccepted && acceptance.genFinalMasscut && acceptance.recAccepted && recTriggered && acceptance.recIdentified && acceptance.recIsolated && acceptance.recMasscut) recMasscut->fill(muon1, muon2, jets30, vars, weight);
		if(isOrigAccepted && acceptance.genFinalAccepted && acceptance.genFinalMasscut && acceptance.recAccepted && recTriggered && acceptance.recIdentified && acceptance.recIsolated && acceptance.recLowmass) recLowmass->fill(muon1, muon2, jets30, vars, weight);
		if(isOrigAccepted && acceptance.genFinalAccepted && acceptance.genFinalMasscut && acceptance.recAccepted && recTriggered && acceptance.recIdentified && acceptance.recMasscut) recMasscutNoIsolation->fill(muon1, muon2, jets30, vars, weight);
	}

	const bool isEmbedded;
	const bool isGenEmbedded;

	std::auto_ptr<Histograms> genFinalAccepted;
	std::auto_ptr<Histograms> origAccepted;
	std::auto_ptr<Histograms> origAndGenAccepted;
	std::auto_ptr<Histograms> genInitialMasscut;
	std::auto_ptr<Histograms> genFinalMasscut;
	std::auto_ptr<Histograms> recAccepted;
	std::auto_ptr<Histograms> recMasscut;
	std::auto_ptr<Histograms> recLowmass;
	std::auto_ptr<Histograms> recMasscutNoIsolation;
};

struct Out
{
	struct JetPtSorter
	{
	public:
		const bool operator()(const Jet& jet1, const Jet& jet2)
		{
			return jet1.p4.Pt() > jet2.p4.Pt();
		}
	};

	Out(TDirectory* dir, bool is_data, bool is_embedded, bool is_gen_embedded, const ZmumuEvtSelEffCorrWeight *zmumuEvtSelWeight, const ZmumuEvtSelEffCorrWeight *zmumuEvtSelWeightEMB):
		isEmbedded(is_embedded),
		isGenEmbedded(is_gen_embedded),
		inclusive(dir, "inclusive", is_data, is_embedded, is_gen_embedded, true)
	{
		// For four-vector transformation study
		if(!is_data && !is_embedded)
		{
			transformInitial.reset(new TransformCategory(dir, "transformInitial"));
			transformFinal.reset(new TransformCategory(dir, "transformFinal"));

			transformInitial_Eta_25_15.reset(new TransformCategory(dir, "transformInitial_Eta_25_15"));
			transformInitial_Eta_15_5.reset(new TransformCategory(dir, "transformInitial_Eta_15_5"));
			transformInitial_Eta_5_M5.reset(new TransformCategory(dir, "transformInitial_Eta_5_M5"));
			transformInitial_Eta_M5_M15.reset(new TransformCategory(dir, "transformInitial_Eta_M5_M15"));
			transformInitial_Eta_M15_M25.reset(new TransformCategory(dir, "transformInitial_Eta_M15_M25"));

			transformInitial_Pt_0_20.reset(new TransformCategory(dir, "transformInitial_Pt_0_20"));
			transformInitial_Pt_20_40.reset(new TransformCategory(dir, "transformInitial_Pt_20_40"));
			transformInitial_Pt_40_inf.reset(new TransformCategory(dir, "transformInitial_Pt_40_inf"));

			transformInitial_ZPt_0_2.reset(new TransformCategory(dir, "transformInitial_ZPt_0_2"));
			transformInitial_ZPt_2_10.reset(new TransformCategory(dir, "transformInitial_ZPt_2_10"));
			transformInitial_ZPt_10_20.reset(new TransformCategory(dir, "transformInitial_ZPt_10_20"));
			transformInitial_ZPt_20_40.reset(new TransformCategory(dir, "transformInitial_ZPt_20_40"));
			transformInitial_ZPt_40_inf.reset(new TransformCategory(dir, "transformInitial_ZPt_40_inf"));
		}

		// Filter FSR on generator level:
		if(!is_data)
			inclusiveNoRadiation.reset(new Category(dir, "inclusiveNoRadiation", is_data, is_embedded, is_gen_embedded, false));

		// For reco-level FSR filter study:
		if(is_embedded && !is_data)
			radiation.reset(new RadiationHistograms(dir, "radiation"));

		// For determining the Zmumu event selection efficiency corrections:
		if(!is_embedded && !is_data)
		{
			zmumuGen.reset(new ZmumuEffHistograms(dir, "zmumuGen"));
			zmumuRec.reset(new ZmumuEffHistograms(dir, "zmumuRec"));
			if(zmumuEvtSelWeight && zmumuEvtSelWeight->hasLevel1()) zmumuRec_PtCorr.reset(new ZmumuEffHistograms(dir, "zmumuRec_PtCorr"));
			if(zmumuEvtSelWeight && zmumuEvtSelWeight->hasLevel2()) zmumuRec_PtEtaCorr.reset(new ZmumuEffHistograms(dir, "zmumuRec_PtEtaCorr"));

			zmumuGenEmb.reset(new ZmumuEffHistograms(dir, "zmumuGenEmb"));
			zmumuRecEmb.reset(new ZmumuEffHistograms(dir, "zmumuRecEmb"));
			if(zmumuEvtSelWeightEMB && zmumuEvtSelWeightEMB->hasLevel1()) zmumuRecEmb_PtCorr.reset(new ZmumuEffHistograms(dir, "zmumuRecEmb_PtCorr"));
			if(zmumuEvtSelWeightEMB && zmumuEvtSelWeightEMB->hasLevel2()) zmumuRecEmb_PtEtaCorr.reset(new ZmumuEffHistograms(dir, "zmumuRecEmb_PtEtaCorr"));
		}

		// Standard monitoring distributions:
		nEvents = new TH1D("nEvents", "Total number of events", 1, 0, 1);
		nPV = new TH1D("nPV", "Number of Primary Vertices", 50, -0.5, 49.5);
		if(!is_data) nTrueInteractions = new TH1D("nTrueInteractions", "Number of True Interactions", 100, -0.5, 99.5);
	}

	void fill(const Muon& posMuon, const Muon& negMuon, const std::vector<Jet>& jets, const TreeVars& vars, const ZmumuEvtSelEffCorrWeight *zmumuEvtSelWeight, const ZmumuEvtSelEffCorrWeight* zmumuEvtSelWeightEMB, float weight)
	{
		if(!isEmbedded)
		{
			const double diMuonPt = (posMuon.p4genInitial + negMuon.p4genInitial).Pt();
			transformInitial->fill(posMuon.p4genInitial, negMuon.p4genInitial, weight);
			transformFinal->fill(posMuon.p4genFinal, negMuon.p4genFinal, weight);

			if(posMuon.p4genInitial.Eta() < 2.5 && posMuon.p4genInitial.Eta() > 1.5) transformInitial_Eta_25_15->fill(posMuon.p4genInitial, negMuon.p4genInitial, weight);
			if(posMuon.p4genInitial.Eta() < 1.5 && posMuon.p4genInitial.Eta() > 0.5) transformInitial_Eta_15_5->fill(posMuon.p4genInitial, negMuon.p4genInitial, weight);
			if(posMuon.p4genInitial.Eta() < 0.5 && posMuon.p4genInitial.Eta() > -0.5) transformInitial_Eta_5_M5->fill(posMuon.p4genInitial, negMuon.p4genInitial, weight);
			if(posMuon.p4genInitial.Eta() < -0.5 && posMuon.p4genInitial.Eta() > -1.5) transformInitial_Eta_M5_M15->fill(posMuon.p4genInitial, negMuon.p4genInitial, weight);
			if(posMuon.p4genInitial.Eta() < -1.5 && posMuon.p4genInitial.Eta() > -2.5) transformInitial_Eta_M15_M25->fill(posMuon.p4genInitial, negMuon.p4genInitial, weight);

			if(posMuon.p4genInitial.Pt() < 20.) transformInitial_Pt_0_20->fill(posMuon.p4genInitial, negMuon.p4genInitial, weight);
			if(posMuon.p4genInitial.Pt() > 20. && posMuon.p4genInitial.Pt() < 40.) transformInitial_Pt_20_40->fill(posMuon.p4genInitial, negMuon.p4genInitial, weight);
			if(posMuon.p4genInitial.Pt() > 40.) transformInitial_Pt_40_inf->fill(posMuon.p4genInitial, negMuon.p4genInitial, weight);

			if(diMuonPt < 2.) transformInitial_ZPt_0_2->fill(posMuon.p4genInitial, negMuon.p4genInitial, weight);
			if(diMuonPt >= 2. && diMuonPt < 10.) transformInitial_ZPt_2_10->fill(posMuon.p4genInitial, negMuon.p4genInitial, weight);
			if(diMuonPt >= 10. && diMuonPt < 20.) transformInitial_ZPt_10_20->fill(posMuon.p4genInitial, negMuon.p4genInitial, weight);
			if(diMuonPt >= 20. && diMuonPt < 40.) transformInitial_ZPt_20_40->fill(posMuon.p4genInitial, negMuon.p4genInitial, weight);
			if(diMuonPt >= 40.) transformInitial_ZPt_40_inf->fill(posMuon.p4genInitial, negMuon.p4genInitial, weight);
		}

		const DiMuonAcceptance acceptance(posMuon, negMuon, vars);

		std::vector<Jet> jets30;
		for(unsigned int i = 0; i < jets.size(); ++i)
			if(jets[i].p4.Pt() > 30)
				jets30.push_back(jets[i]);
		std::sort(jets30.begin(), jets30.end(), JetPtSorter());

		// If we have a reco-embedded sample, then we need to have Zmumu efficiency corrections
		assert(!isEmbedded || isGenEmbedded || (zmumuEvtSelWeight && zmumuEvtSelWeight->hasLevel2()));
		if(isEmbedded && !isGenEmbedded)
			weight *= zmumuEvtSelWeightEMB->GetPtEtaWeight(posMuon, negMuon, isEmbedded);

		inclusive.fill(acceptance, posMuon, negMuon, jets30, vars, weight);

		// No-radiation category: Filter FSR on generator level; for embedded samples filter it
		// based on the original input muons.
		if(inclusiveNoRadiation.get())
		{
			if(isEmbedded)
			{
				assert(posMuon.hasOrigGen && negMuon.hasOrigGen);
				if(lorentzEqual(posMuon.p4origGenInitial, posMuon.p4origGenFinal) && lorentzEqual(negMuon.p4origGenInitial, negMuon.p4origGenFinal))
					inclusiveNoRadiation->fill(acceptance, posMuon, negMuon, jets30, vars, weight);
			}
			else
			{
				assert(posMuon.hasGen && negMuon.hasGen);
				if(lorentzEqual(posMuon.p4genInitial, posMuon.p4genFinal) && lorentzEqual(negMuon.p4genInitial, negMuon.p4genFinal))
					inclusiveNoRadiation->fill(acceptance, posMuon, negMuon, jets30, vars, weight);
			}
		}

		// reco-level radiation filter study
		if(radiation.get())
			if(acceptance.origAccepted)
				radiation->fill(posMuon, negMuon, vars, weight);

		assert(!isEmbedded || isGenEmbedded || (posMuon.hasOrigRec && negMuon.hasOrigRec));
		assert(!isEmbedded || !isGenEmbedded || (posMuon.hasOrigGen && negMuon.hasOrigGen));

		// For obtaining Zmumu selection efficiency correction:
		if(acceptance.zmumuGen && zmumuGen.get()) zmumuGen->fill(posMuon, negMuon, 1.0f);
		if(acceptance.zmumuRec && zmumuRec.get()) zmumuRec->fill(posMuon, negMuon, 1.0f);
		if(acceptance.zmumuRec && zmumuRec_PtCorr.get()) zmumuRec_PtCorr->fill(posMuon, negMuon, zmumuEvtSelWeight->GetPtWeight(posMuon, negMuon, isEmbedded));
		if(acceptance.zmumuRec && zmumuRec_PtEtaCorr.get()) zmumuRec_PtEtaCorr->fill(posMuon, negMuon, zmumuEvtSelWeight->GetPtEtaWeight(posMuon, negMuon, isEmbedded));

		if(acceptance.zmumuGenEmb && zmumuGenEmb.get()) zmumuGenEmb->fill(posMuon, negMuon, 1.0f);
		if(acceptance.zmumuRecEmb && zmumuRecEmb.get()) zmumuRecEmb->fill(posMuon, negMuon, 1.0f);
		if(acceptance.zmumuRecEmb && zmumuRecEmb_PtCorr.get()) zmumuRecEmb_PtCorr->fill(posMuon, negMuon, zmumuEvtSelWeightEMB->GetPtWeight(posMuon, negMuon, isEmbedded));
		if(acceptance.zmumuRecEmb && zmumuRecEmb_PtEtaCorr.get()) zmumuRecEmb_PtEtaCorr->fill(posMuon, negMuon, zmumuEvtSelWeightEMB->GetPtEtaWeight(posMuon, negMuon, isEmbedded));
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

	const bool isEmbedded;
	const bool isGenEmbedded;

	std::auto_ptr<TransformCategory> transformInitial;
	std::auto_ptr<TransformCategory> transformFinal;

	std::auto_ptr<TransformCategory> transformInitial_Eta_25_15;
	std::auto_ptr<TransformCategory> transformInitial_Eta_15_5;
	std::auto_ptr<TransformCategory> transformInitial_Eta_5_M5;
	std::auto_ptr<TransformCategory> transformInitial_Eta_M5_M15;
	std::auto_ptr<TransformCategory> transformInitial_Eta_M15_M25;

	std::auto_ptr<TransformCategory> transformInitial_Pt_0_20;
	std::auto_ptr<TransformCategory> transformInitial_Pt_20_40;
	std::auto_ptr<TransformCategory> transformInitial_Pt_40_inf;

	std::auto_ptr<TransformCategory> transformInitial_ZPt_0_2;
	std::auto_ptr<TransformCategory> transformInitial_ZPt_2_10;
	std::auto_ptr<TransformCategory> transformInitial_ZPt_10_20;
	std::auto_ptr<TransformCategory> transformInitial_ZPt_20_40;
	std::auto_ptr<TransformCategory> transformInitial_ZPt_40_inf;

	Category inclusive;
	std::auto_ptr<Category> inclusiveNoRadiation;
	std::auto_ptr<RadiationHistograms> radiation;

	std::auto_ptr<ZmumuEffHistograms> zmumuGen;
	std::auto_ptr<ZmumuEffHistograms> zmumuRec;
	std::auto_ptr<ZmumuEffHistograms> zmumuRec_PtCorr;
	std::auto_ptr<ZmumuEffHistograms> zmumuRec_PtEtaCorr;

	std::auto_ptr<ZmumuEffHistograms> zmumuGenEmb;
	std::auto_ptr<ZmumuEffHistograms> zmumuRecEmb;
	std::auto_ptr<ZmumuEffHistograms> zmumuRecEmb_PtCorr;
	std::auto_ptr<ZmumuEffHistograms> zmumuRecEmb_PtEtaCorr;

	TH1D* nEvents;
	TH1D* nPV;
	TH1D* nTrueInteractions;
};

void process_file(const char* filename, Out& out, bool is_data, bool is_embedded, bool is_gen_embedded, const ZmumuEvtSelEffCorrWeight *zmumuEvtSelWeight, const ZmumuEvtSelEffCorrWeight *zmumuEvtSelWeightEMB, int pdg_filter)
{
	TreeVars vars;

	std::string directory = "diMuons";
	TFile* file = new TFile(filename, "READ");
	TTree* tree = dynamic_cast<TTree*>(file->Get((directory + "/MuonTree").c_str()));

	tree->SetBranchAddress("Run", &vars.run);
	tree->SetBranchAddress("Lumi", &vars.lumi);
	tree->SetBranchAddress("Event", &vars.event);

	if(is_embedded)
	{
		tree->SetBranchAddress("MuonRadiationFilter", &vars.muonRadiationFilter);
		tree->SetBranchAddress("MuonRadiationFilter2Sel1", &vars.muonRadiationFilter2Sel1);
		tree->SetBranchAddress("MuonRadiationFilter2Sel2", &vars.muonRadiationFilter2Sel2);
		tree->SetBranchAddress("MuonRadiationFilter2Sel3", &vars.muonRadiationFilter2Sel3);
	}

	tree->SetBranchAddress("NPV", &vars.nPV);
	if(!is_data) tree->SetBranchAddress("NTrueInteractions", &vars.nTrueInteractions);
	if(!is_data) tree->SetBranchAddress("NGenMEJets", &vars.nGenMEJets);
	tree->SetBranchAddress("HltMu17Mu8", &vars.hltMu17Mu8);
	if(is_embedded && !is_gen_embedded) tree->SetBranchAddress("OrigHltMu17Mu8", &vars.origHltMu17Mu8);

	tree->SetBranchAddress("CaloMet", &vars.caloMet);
	tree->SetBranchAddress("CaloMetPhi", &vars.caloMetPhi);
	tree->SetBranchAddress("PfMet", &vars.pfMet);
	tree->SetBranchAddress("PfMetPhi", &vars.pfMetPhi);
	if(is_embedded && !is_gen_embedded)
	{
		tree->SetBranchAddress("OrigCaloMet", &vars.origCaloMet);
		tree->SetBranchAddress("OrigCaloMetPhi", &vars.origCaloMetPhi);
		tree->SetBranchAddress("OrigPfMet", &vars.origPfMet);
		tree->SetBranchAddress("OrigPfMetPhi", &vars.origPfMetPhi);
	}

	tree->SetBranchAddress("PosMuonPt", &vars.PosMuon.Pt);
	tree->SetBranchAddress("PosMuonEta", &vars.PosMuon.Eta);
	tree->SetBranchAddress("PosMuonPhi", &vars.PosMuon.Phi);
	tree->SetBranchAddress("PosMuonE", &vars.PosMuon.E);
	tree->SetBranchAddress("PosMuonChargedParticlePtPfIso04", &vars.PosMuon.ChargedParticlePtPfIso04);
	tree->SetBranchAddress("PosMuonChargedHadronPfIso04", &vars.PosMuon.ChargedHadronPtPfIso04);
	tree->SetBranchAddress("PosMuonNeutralHadronEtPfIso04", &vars.PosMuon.NeutralHadronEtPfIso04);
	tree->SetBranchAddress("PosMuonPhotonEtPfIso04", &vars.PosMuon.PhotonEtPfIso04);
	tree->SetBranchAddress("PosMuonPUPtPfIso04", &vars.PosMuon.PUPtPfIso04);
	tree->SetBranchAddress("PosMuonQuality", &vars.PosMuon.Quality);
	tree->SetBranchAddress("PosMuonHltMu17Mu8Leg8", &vars.PosMuon.HltMu17Mu8Leg8);
	tree->SetBranchAddress("PosMuonHltMu17Mu8Leg17", &vars.PosMuon.HltMu17Mu8Leg17);

	if(!is_data || is_embedded)
	{
		tree->SetBranchAddress("GenPosChargedLepPt", &vars.PosMuon.GenPt);
		tree->SetBranchAddress("GenPosChargedLepEta", &vars.PosMuon.GenEta);
		tree->SetBranchAddress("GenPosChargedLepPhi", &vars.PosMuon.GenPhi);
		tree->SetBranchAddress("GenPosChargedLepE", &vars.PosMuon.GenE);
		tree->SetBranchAddress("GenPosChargedLepFinalPt", &vars.PosMuon.GenFinalPt);
		tree->SetBranchAddress("GenPosChargedLepFinalEta", &vars.PosMuon.GenFinalEta);
		tree->SetBranchAddress("GenPosChargedLepFinalPhi", &vars.PosMuon.GenFinalPhi);
		tree->SetBranchAddress("GenPosChargedLepFinalE", &vars.PosMuon.GenFinalE);
		tree->SetBranchAddress("GenPosChargedLepPDG", &vars.PosMuon.GenChargedLepPDG);
	}

	if(is_embedded)
	{
		if(!is_gen_embedded)
		{
			tree->SetBranchAddress("OrigPosMuonPt", &vars.PosMuon.OrigPt);
			tree->SetBranchAddress("OrigPosMuonEta", &vars.PosMuon.OrigEta);
			tree->SetBranchAddress("OrigPosMuonPhi", &vars.PosMuon.OrigPhi);
			tree->SetBranchAddress("OrigPosMuonE", &vars.PosMuon.OrigE);
		}

		if(!is_data)
		{
			tree->SetBranchAddress("OrigGenPosChargedLepPt", &vars.PosMuon.OrigGenPt);
			tree->SetBranchAddress("OrigGenPosChargedLepEta", &vars.PosMuon.OrigGenEta);
			tree->SetBranchAddress("OrigGenPosChargedLepPhi", &vars.PosMuon.OrigGenPhi);
			tree->SetBranchAddress("OrigGenPosChargedLepE", &vars.PosMuon.OrigGenE);
			tree->SetBranchAddress("OrigGenPosChargedLepFinalPt", &vars.PosMuon.OrigGenFinalPt);
			tree->SetBranchAddress("OrigGenPosChargedLepFinalEta", &vars.PosMuon.OrigGenFinalEta);
			tree->SetBranchAddress("OrigGenPosChargedLepFinalPhi", &vars.PosMuon.OrigGenFinalPhi);
			tree->SetBranchAddress("OrigGenPosChargedLepFinalE", &vars.PosMuon.OrigGenFinalE);
			tree->SetBranchAddress("OrigGenPosChargedLepPDG", &vars.PosMuon.OrigGenChargedLepPDG);
		}
	}

	tree->SetBranchAddress("NegMuonPt", &vars.NegMuon.Pt);
	tree->SetBranchAddress("NegMuonEta", &vars.NegMuon.Eta);
	tree->SetBranchAddress("NegMuonPhi", &vars.NegMuon.Phi);
	tree->SetBranchAddress("NegMuonE", &vars.NegMuon.E);
	tree->SetBranchAddress("NegMuonChargedParticlePtPfIso04", &vars.NegMuon.ChargedParticlePtPfIso04);
	tree->SetBranchAddress("NegMuonChargedHadronPfIso04", &vars.NegMuon.ChargedHadronPtPfIso04);
	tree->SetBranchAddress("NegMuonNeutralHadronEtPfIso04", &vars.NegMuon.NeutralHadronEtPfIso04);
	tree->SetBranchAddress("NegMuonPhotonEtPfIso04", &vars.NegMuon.PhotonEtPfIso04);
	tree->SetBranchAddress("NegMuonPUPtPfIso04", &vars.NegMuon.PUPtPfIso04);
	tree->SetBranchAddress("NegMuonQuality", &vars.NegMuon.Quality);
	tree->SetBranchAddress("NegMuonHltMu17Mu8Leg8", &vars.NegMuon.HltMu17Mu8Leg8);
	tree->SetBranchAddress("NegMuonHltMu17Mu8Leg17", &vars.NegMuon.HltMu17Mu8Leg17);

	if(!is_data || is_embedded)
	{
		tree->SetBranchAddress("GenNegChargedLepPt", &vars.NegMuon.GenPt);
		tree->SetBranchAddress("GenNegChargedLepEta", &vars.NegMuon.GenEta);
		tree->SetBranchAddress("GenNegChargedLepPhi", &vars.NegMuon.GenPhi);
		tree->SetBranchAddress("GenNegChargedLepE", &vars.NegMuon.GenE);
		tree->SetBranchAddress("GenNegChargedLepFinalPt", &vars.NegMuon.GenFinalPt);
		tree->SetBranchAddress("GenNegChargedLepFinalEta", &vars.NegMuon.GenFinalEta);
		tree->SetBranchAddress("GenNegChargedLepFinalPhi", &vars.NegMuon.GenFinalPhi);
		tree->SetBranchAddress("GenNegChargedLepFinalE", &vars.NegMuon.GenFinalE);
		tree->SetBranchAddress("GenNegChargedLepPDG", &vars.NegMuon.GenChargedLepPDG);
	}

	if(is_embedded)
	{
		if(!is_gen_embedded)
		{
			tree->SetBranchAddress("OrigNegMuonPt", &vars.NegMuon.OrigPt);
			tree->SetBranchAddress("OrigNegMuonEta", &vars.NegMuon.OrigEta);
			tree->SetBranchAddress("OrigNegMuonPhi", &vars.NegMuon.OrigPhi);
			tree->SetBranchAddress("OrigNegMuonE", &vars.NegMuon.OrigE);
		}

		if(!is_data)
		{
			tree->SetBranchAddress("OrigGenNegChargedLepPt", &vars.NegMuon.OrigGenPt);
			tree->SetBranchAddress("OrigGenNegChargedLepEta", &vars.NegMuon.OrigGenEta);
			tree->SetBranchAddress("OrigGenNegChargedLepPhi", &vars.NegMuon.OrigGenPhi);
			tree->SetBranchAddress("OrigGenNegChargedLepE", &vars.NegMuon.OrigGenE);
			tree->SetBranchAddress("OrigGenNegChargedLepFinalPt", &vars.NegMuon.OrigGenFinalPt);
			tree->SetBranchAddress("OrigGenNegChargedLepFinalEta", &vars.NegMuon.OrigGenFinalEta);
			tree->SetBranchAddress("OrigGenNegChargedLepFinalPhi", &vars.NegMuon.OrigGenFinalPhi);
			tree->SetBranchAddress("OrigGenNegChargedLepFinalE", &vars.NegMuon.OrigGenFinalE);
			tree->SetBranchAddress("OrigGenNegChargedLepPDG", &vars.NegMuon.OrigGenChargedLepPDG);
		}
	}

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

	TH1F* h_nPV = dynamic_cast<TH1F*>(file->Get((directory + "/h_nPV").c_str()));
	out.add_events(static_cast<unsigned int>(h_nPV->GetEntries()));
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

		// PDG filter
		if(!is_data && pdg_filter != 0)
		{
			if(!is_embedded && (abs(vars.NegMuon.GenChargedLepPDG) != pdg_filter || abs(vars.PosMuon.GenChargedLepPDG) != pdg_filter))
				continue;
			if(is_embedded && (abs(vars.NegMuon.OrigGenChargedLepPDG) != pdg_filter || abs(vars.PosMuon.OrigGenChargedLepPDG) != pdg_filter))
				continue;
		}

		// Trigger
		//if(!is_embedded && !vars.hltMu17Mu8) continue;

		// Muons
		Muon posMuon(vars.PosMuon, +1);
		Muon negMuon(vars.NegMuon, -1);

		// There are some rare events for which the Z four vector is not preserved.
		// For GenEmbedding, this fails in ~1-2% of the cases. I am not sure why.
		// Maybe the genMuonsFromZs in the embedding code use a slightly different
		// algorithm for finding the muons after FSR than we do.
		if(is_embedded)
		{
			assert(is_gen_embedded || (posMuon.hasOrigRec && negMuon.hasOrigRec));
			assert(!is_gen_embedded || (posMuon.hasOrigGen && negMuon.hasOrigGen));
			assert(posMuon.hasGen && negMuon.hasGen);

			TLorentzVector orig =
				is_gen_embedded ?
				posMuon.p4origGenFinal + negMuon.p4origGenFinal :
				posMuon.p4orig + negMuon.p4orig;
			TLorentzVector gen = posMuon.p4genInitial + negMuon.p4genInitial;
			TLorentzVector diff = orig - gen;
			if(fabs(diff.Px()) > 1e-3 || fabs(diff.Py()) > 1e-3 || fabs(diff.Pz()) > 1e-3 || fabs(diff.E()) > 1e-3)
			{
				++Filter;
				continue;
			}
		}

		// Trigger matching
		//if(!is_embedded && !(posMuon.HltMu17Mu8Leg8 && negMuon.HltMu17Mu8Leg17) && !(posMuon.HltMu17Mu8Leg17 && negMuon.HltMu17Mu8Leg8)) continue;

		// Jets
		std::vector<Jet> jets;
		for(unsigned int i = 0; i < vars.nJets; ++i)
			jets.push_back(Jet(vars, i));

		double weight = 1.0;

		// Output
		out.fill(posMuon, negMuon, jets, vars, zmumuEvtSelWeight, zmumuEvtSelWeightEMB, weight);
	}

	std::cout << "Filtered: " << Filter << "/" << CountTotal << std::endl;
	delete file;
}

void process_dataset(const char* dataset, const std::vector<std::string>& files, unsigned int pdg_filter)
{
	const bool is_data = (strstr(dataset, "Data") != NULL);
	const bool is_embedded = (strstr(dataset, "Embed") != NULL);
	const bool is_gen_embedded = (strstr(dataset, "GenEmbed") != NULL);

	std::stringstream out_filename_stream;
	out_filename_stream << "out-" << dataset;
	if(pdg_filter != 0) out_filename_stream << "-filter" << pdg_filter;
	out_filename_stream << ".root";
	const std::string out_filename = out_filename_stream.str();

	std::auto_ptr<ZmumuEvtSelEffCorrWeight> zmumuEvtSelWeight, zmumuEvtSelWeightEMB;
	zmumuEvtSelWeight.reset(new ZmumuEvtSelEffCorrWeight("../data/ZmumuEvtSelEffCorrWeights.root", "")); // also for normal MC, so we can calculate L2 (eta-eta) corrections
	zmumuEvtSelWeightEMB.reset(new ZmumuEvtSelEffCorrWeight("../data/ZmumuEvtSelEffCorrWeights.root", "EMB_")); // also for normal MC, so we can calculate L2 (eta-eta) corrections

	std::cout << dataset << "... => " << out_filename << std::endl;
	std::cout << "\tData: " << is_data << std::endl;
	std::cout << "\tEmbedded: " << is_embedded << std::endl;
	std::cout << "\tGenEmbedded: " << is_gen_embedded << std::endl;

	TFile* outfile = new TFile(out_filename.c_str(), "RECREATE");
	Out out(outfile, is_data, is_embedded, is_gen_embedded, zmumuEvtSelWeight.get(), zmumuEvtSelWeightEMB.get());

	for(std::vector<std::string>::const_iterator iter = files.begin(); iter != files.end(); ++iter)
	{
		std::cout << "\t[" << (iter - files.begin() + 1) << "/" << files.size() << "] " << *iter << "..." << std::endl;
		process_file(iter->c_str(), out, is_data, is_embedded, is_gen_embedded, zmumuEvtSelWeight.get(), zmumuEvtSelWeightEMB.get(), pdg_filter);
	}

	outfile->Write();
	delete outfile;
}

int main(int argc, char* argv[])
{
	TH1::SetDefaultSumw2(true);

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
		process_dataset(iter->first.c_str(), iter->second, pdg);

	return 0;
}
