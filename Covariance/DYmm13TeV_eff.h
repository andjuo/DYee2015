#ifndef DYmm13TeV_eff_h
#define DYmm13TeV_eff_h

#include "inputs.h"
#include "DYbinning.h"
#include <TLorentzVector.h>

namespace DYtools {
//const int EtaPtFIMax=25; // assummed maximum number of (eta,pt) bins
const int EtaPtFIMax=50; // assummed maximum number of (eta,pt) bins

// -------------------------------------------------------------

inline
  int FlatIndex(const TH2D* h2BinDefs, double eta, double pt, int allowOverflow) {
  int ibin=h2BinDefs->GetXaxis()->FindBin(eta);
  int jbin=h2BinDefs->GetYaxis()->FindBin(pt);
  if (allowOverflow && (jbin==h2BinDefs->GetNbinsY()+1)) jbin=h2BinDefs->GetNbinsY();
  int fi= (jbin-1) + (ibin-1)*h2BinDefs->GetNbinsY();
  if ((ibin==0) || (jbin==0)) fi=-1;
  else if ((ibin==h2BinDefs->GetNbinsX()+1) ||
	   (jbin==h2BinDefs->GetNbinsY()+1)) fi=-1;
  if (0) {
    std::cout << "eta=" << eta << ", pt=" << pt << ", ibin=" << ibin
	      << ", jbin=" << jbin << ", fi=" << fi << "\n";
  }
  return fi;
}

// -------------------------------------------------------------

inline
  int FlatIndex(const TH2D* h2BinDefs, int iEta, int iPt, int allowOverflow) {
  int ibin=iEta;
  int jbin=iPt;
  if (allowOverflow && (jbin==h2BinDefs->GetNbinsY()+1)) jbin=h2BinDefs->GetNbinsY();
  int fi= (jbin-1) + (ibin-1)*h2BinDefs->GetNbinsY();
  if ((ibin==0) || (jbin==0)) fi=-1;
  else if ((ibin==h2BinDefs->GetNbinsX()+1) ||
	   (jbin==h2BinDefs->GetNbinsY()+1)) fi=-1;
  //std::cout << "eta=" << eta << ", pt=" << pt << ", ibin=" << ibin
  //	    << ", jbin=" << jbin << ", fi=" << fi << "\n";

  if (0) {
    std::cout << " FlatIndex(iEta=" << iEta << ", iPt=" << iPt << ") "
	      << " eta=" << h2BinDefs->GetXaxis()->GetBinCenter(iEta)
	      << ", pT=" << h2BinDefs->GetYaxis()->GetBinCenter(iPt)
	      << ", fi=" << fi
	      << "\n";
  }
  return fi;
}

// -------------------------------------------------------------

inline
int FlatIndex(const TH2D* h2BinDefs, const TLorentzVector *v, int allowOverflow=0) {
  return FlatIndex(h2BinDefs,v->Eta(),v->Pt(), allowOverflow);
}

// -------------------------------------------------------------

inline
int GetValue(const TH2D* h2, double eta, double pt, double &v, double &err,
	     double defaultVal=0., double defaultErr=0.)
{
  int ibin=h2->GetXaxis()->FindBin(eta);
  int jbin=h2->GetYaxis()->FindBin(pt);
  int ok=1;
  if ((ibin==0) || (jbin==0)) ok=0;
  else if ((ibin==h2->GetNbinsX()+1) ||
	   (jbin==h2->GetNbinsY()+1)) ok=0;
  if (ok) {
    v=  h2->GetBinContent(ibin,jbin);
    err=h2->GetBinError  (ibin,jbin);
  }
  else {
    v=defaultVal;
    err=defaultErr;
  }
  //std::cout << "eta=" << eta << ", pt=" << pt << ", ibin=" << ibin
  //	    << ", jbin=" << jbin << ", ok=" << ok
  //	    << ": " << v << " +- " << err
  //	    << "\n";
  return ok;
}

// -------------------------------------------------------------

inline
int GetValueIdx(const TH2D* h2, int iEta, int iPt, double &v, double &err,
	     double defaultVal=0., double defaultErr=0.)
{
  int ibin=iEta;
  int jbin=iPt;
  int ok=1;
  if ((ibin==0) || (jbin==0)) ok=0;
  else if ((ibin==h2->GetNbinsX()+1) ||
	   (jbin==h2->GetNbinsY()+1)) ok=0;
  if (ok) {
    v=  h2->GetBinContent(ibin,jbin);
    err=h2->GetBinError  (ibin,jbin);
  }
  else {
    v=defaultVal;
    err=defaultErr;
  }
  //std::cout << "eta=" << eta << ", pt=" << pt << ", ibin=" << ibin
  //	    << ", jbin=" << jbin << ", ok=" << ok
  //	    << ": " << v << " +- " << err
  //	    << "\n";
  return ok;
}


}; // namespace DYtools

// -------------------------------------------------------------
// -------------------------------------------------------------

class DYTnPEff_t {
public:
  TH2D *h2Eff_RecoID_Data, *h2Eff_Iso_Data;
  TH2D *h2Eff_HLT4p2_Data, *h2Eff_HLT4p3_Data;
  TH2D *h2Eff_RecoID_MC, *h2Eff_Iso_MC;
  TH2D *h2Eff_HLT4p2_MC, *h2Eff_HLT4p3_MC;
protected:
  std::vector<TH2D*> h2VReco, h2VIso, h2VHLT;
  //double fDefVal, fDefErr;
  mutable double fEffReco, fEffIso, fEffHlt1, fEffHlt2;
  //double fEffRecoErr, fEffIsoErr, fEffHlt1Err, fEffHlt2Err;
  int fElChannel;

protected:
  int ptrsOk(const std::vector<TH2D*> &v) const
  {
    int ok=1;
    for (unsigned int i=0; ok && (i<v.size()); i++)
      if (!v[i]) ok=0;
    return ok;
  }

public:
  DYTnPEff_t(int setElChannel=0);
  DYTnPEff_t(const DYTnPEff_t &e, TString tag);

  int isElChannel() const { return fElChannel; }

  int ptrsOk() const {
    int ok= (
	     (h2VReco.size()==2) && (h2VIso.size()==2) &&
	     (h2VHLT.size()==4)
	     ) ? 1:0;
    if (ok) ok=(ptrsOk(h2VReco) && ptrsOk(h2VIso) && ptrsOk(h2VHLT)) ? 1:0;
    return ok;
  }

  double totalEffIdx(int ibin1, int jbin1, int ibin2, int jbin2,
		     int mc, int hlt4p3) const
  {
    fEffReco=
      h2VReco[mc]->GetBinContent(ibin1,jbin1) *
      h2VReco[mc]->GetBinContent(ibin2,jbin2);
    fEffIso=
      h2VIso[mc]->GetBinContent(ibin1,jbin1) *
      h2VIso[mc]->GetBinContent(ibin2,jbin2);
    fEffHlt1= h2VHLT[mc+2*hlt4p3]->GetBinContent(ibin1,jbin1);
    fEffHlt2= h2VHLT[mc+2*hlt4p3]->GetBinContent(ibin2,jbin2);
    double effHlt= 1 - (1-fEffHlt1)*(1-fEffHlt2);
    double eff= fEffReco *fEffIso * effHlt;
    if (0) {
      std::cout << "totalEffIdx(" << ibin1 << "," << jbin1 << "; "
		<< ibin2 << "," << jbin2 << "; "
		<< h2VReco[mc]->GetBinContent(ibin1,jbin1) << " x "
		<< h2VReco[mc]->GetBinContent(ibin2,jbin2) << "; "
		<< h2VIso[mc]->GetBinContent(ibin1,jbin1) << " x "
		<< h2VIso[mc]->GetBinContent(ibin2,jbin2) << "; "
		<< h2VHLT[mc+2*hlt4p3]->GetBinContent(ibin1,jbin1) << " x "
		<< h2VHLT[mc+2*hlt4p3]->GetBinContent(ibin2,jbin2) << "; "
		<< "tot=" << eff << "\n";
    }
    return eff;
  }

  double totalEff(double eta1, double pt1, double eta2, double pt2,
		  int mc, int hlt4p3) const
  {
    int ibin1=h2Eff_RecoID_Data->GetXaxis()->FindBin(eta1);
    int ibin2=h2Eff_RecoID_Data->GetXaxis()->FindBin(eta2);
    int jbin1=h2Eff_RecoID_Data->GetYaxis()->FindBin(pt1);
    int jbin2=h2Eff_RecoID_Data->GetYaxis()->FindBin(pt2);
    if ((ibin1==0) || (jbin1==0) || (ibin2==0) || (jbin2==0)) return 0.;
    if (ibin1==h2Eff_RecoID_Data->GetNbinsX()+1) ibin1--;
    if (ibin2==h2Eff_RecoID_Data->GetNbinsX()+1) ibin2--;
    if (jbin1==h2Eff_RecoID_Data->GetNbinsY()+1) jbin1--;
    if (jbin2==h2Eff_RecoID_Data->GetNbinsY()+1) jbin2--;
    return totalEffIdx(ibin1,jbin1,ibin2,jbin2, mc, hlt4p3);
  }

  double totalEff(const TLorentzVector *v1, const TLorentzVector *v2,
		  int mc, int hlt4p3) const
  { return totalEff(v1->Eta(),v1->Pt(),v2->Eta(),v2->Pt(),mc,hlt4p3); }


  double scaleFactor(double eta1, double pt1, double eta2, double pt2,
		     int hlt4p3) const
  {
    //std::cout << "scaleFactor( " << Form(" (%2.2lf,%2.0lf)", eta1,pt1)
    //	      << Form(" (%2.2lf,%2.0lf)", eta2,pt2)
    //	      << " hlt4p3=" << hlt4p3 << "\n";
    double effData= totalEff(eta1,pt1,eta2,pt2, 0,hlt4p3);
    //std::cout << " totalEffData=" << effData << "\n";
    double effMC= totalEff(eta1,pt1,eta2,pt2, 1,hlt4p3);
    //std::cout << " totalEffMC=" << effMC << "\n";
    double rho= (effMC==0.) ? 0. : effData/effMC;
    //std::cout << " rho=" << rho << "\n";
    return rho;
  }

  double scaleFactor(const TLorentzVector *v1, const TLorentzVector *v2,
		     int hlt4p3) const
  { return scaleFactor(v1->Eta(),v1->Pt(),v2->Eta(),v2->Pt(),hlt4p3); }

  double scaleFactorIdx(int ibin1, int jbin1, int ibin2, int jbin2, int hlt4p3) const
  {
    //std::cout << "scaleFactorIdx( " << Form(" (%2d,%2d)", ibin1,jbin1)
    //	      << Form(" (%2d,%2d)", ibin2,jbin2)
    //	      << " hlt4p3=" << hlt4p3 << "\n";
    double effData= totalEffIdx(ibin1,jbin1,ibin2,jbin2, 0,hlt4p3);
    //std::cout << " totalEffData=" << effData << "\n";
    double effMC= totalEffIdx(ibin1,jbin1,ibin2,jbin2, 1,hlt4p3);
    //std::cout << " totalEffMC=" << effMC << "\n";
    double rho= (effMC==0.) ? 0. : effData/effMC;
    //std::cout << " rho=" << rho << "\n";
    return rho;
  }

  int updateVectors(); // copy TH2D* pointers to arrays
  int assign(const DYTnPEff_t &e, TString tag);
  int randomize(const DYTnPEff_t &e, TString tag);
  int randomize(const DYTnPEff_t &e_errPos, const DYTnPEff_t &e_errNeg,
		TString tag);

  int checkBinning() const;
  void printNumbers() const;
  void printEffRatios(const DYTnPEff_t &e, int compareErrs=0) const;

  // internal save/load
  int save(TFile &fout, TString subdirTag="");
  int load(TFile &fin, TString subdir="DYTnPEff", TString tag="");

  // load from Kyeongpil's files
  int load_DYMM_TnP(TFile &fin, int version, TString tag="");
  int load_DYMM_TnP_asymmEff(TFile &fin, int version, int posErr, TString tag="");

};

// -------------------------------------------------------------
// -------------------------------------------------------------

// EventSpace_t contains (eta1,pt1)x(eta2,pt2) information for each mass bin

class EventSpace_t {
  TString fName;
  const TH2D* fh2EffBinDef;
  std::vector<TH2D*> fh2ESV;
public:
  EventSpace_t(TString setName="mainES", TH2D *set_effBinDef=NULL);
  EventSpace_t(TString setName, const EventSpace_t &es);

  int checkPtrs() const;
  TString name() const { return fName; }
  const TH2D* h2EffBinDef() const { return fh2EffBinDef; }
  unsigned int size() const { return fh2ESV.size(); }
  const std::vector<TH2D*> h2ESvec() const { return fh2ESV; }

  template<class idx_t>
    const TH2D* h2ES(idx_t iMass) const { return fh2ESV[iMass]; }

  const TH2D* h2ES(double mass) const {
    int iMass= DYtools::massIdx(mass);
    if (iMass==-1) return NULL;
    if (iMass>int(fh2ESV.size())) { std::cout << "not ready\n"; return NULL; }
    return fh2ESV[iMass];
  }

  int assign(const EventSpace_t &es);

  int fill(const TLorentzVector *v1, const TLorentzVector *v2, double weight) {
    const int allow_overflow=1;
    int fi1= DYtools::FlatIndex(fh2EffBinDef, v1, allow_overflow);
    int fi2= DYtools::FlatIndex(fh2EffBinDef, v2, allow_overflow);
    //std::cout << "fi1=" << fi1 << ", fi2=" << fi2 << "\n";
    if ((fi1<0) || (fi2<0)) {
      std::cout << "negative fi\n";
      return 0;
    }
    if ((fi1>DYtools::EtaPtFIMax) || (fi2>DYtools::EtaPtFIMax)) {
      std::cout << "too small EtaPtFIMax=" << DYtools::EtaPtFIMax
		<< ", fi1=" << fi1 << ", fi2=" << fi2 << "\n";
      return 0;
    }
    double m= (*v1 + *v2).M();
    int iMass= DYtools::massIdx(m);
    //std::cout << "fill m=" << m << ", iMass=" << iMass << ", weight=" << weight << "\n";
    if (iMass==-1) return 0;
    fh2ESV[iMass]->Fill(fi1,fi2, weight);
    return 1;
  }

  TH1D* calculateScaleFactor(const DYTnPEff_t &eff, int hlt4p3,
			     TString hName, TString hTitle) const;

  std::vector<TH1D*> avgAxisValues(std::vector<TH2D*> *h2ptSpace=NULL,
				   std::vector<TH2D*> *h2etaSpace=NULL) const;

  int save(TFile &fout);
  int load(TFile &fin, TString subdir);

 protected:
  int prepareVector();

};

// -------------------------------------------------------------
// -------------------------------------------------------------



#endif
