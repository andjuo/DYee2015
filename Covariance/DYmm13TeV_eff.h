#ifndef DYmm13TeV_eff_h
#define DYmm13TeV_eff_h

#include "inputs.h"
#include "DYbinning.h"
#include <TLorentzVector.h>

namespace DYtools {
//const int EtaPtFIMax=25; // assummed maximum number of (eta,pt) bins
#ifndef def_fewMassBins
  const int EtaPtFIMax=50; // assummed maximum number of (eta,pt) bins
#else
  const int EtaPtFIMax=10; // debugging with a small number of bins
#endif

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
  typedef enum { _kind_RecoID=0, _kind_Iso, _kind_HLT4p2, _kind_HLT4p3
  } TEffKind_t;
  typedef enum { _misc_empty=0, _misc_hlt4p2_plain=1, _misc_hlt4p3_plain,
		 _misc_hlt4p2_leadLep, _misc_hlt4p3_leadLep,
		 _misc_Last } TMiscStudyFlag_t;
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
  ~DYTnPEff_t();

  int isElChannel() const { return fElChannel; }

  std::vector<TH2D*>* h2vecPtr(int isrc) {
    if (isrc==0) return &h2VReco;
    if (isrc==1) return &h2VIso;
    if (isrc==2) return &h2VHLT;
    return NULL;
  }

  const std::vector<TH2D*>* h2vecPtr(int isrc) const {
    if (isrc==0) return &h2VReco;
    if (isrc==1) return &h2VIso;
    if (isrc==2) return &h2VHLT;
    return NULL;
  }

  int setHisto(int kind, int isMC, TH2D *h2);
  int copyHisto(int kind, int isMC, const TH2D *h2);
  const TH2D* getHisto(int ikind, int isMC) const;
  int fullList_setHisto(unsigned int i, TH2D *h2);

  unsigned int fullListSize() const
  { return (h2VReco.size()+h2VIso.size()+h2VHLT.size()); }

  TH2D* h2fullList(unsigned int i)
  {
    if (i>=h2VReco.size()) i-=h2VReco.size();
    else return h2VReco[i];
    if (i>=h2VIso.size()) i-=h2VIso.size();
    else return h2VIso[i];
    if (i>=h2VHLT.size()) {
      std::cout << "h2fullList: i=" << i << " is too large\n";
      return NULL;
    }
    return h2VHLT[i];
  }

  const TH2D* h2fullList(unsigned int i) const
  {
    if (i>=h2VReco.size()) i-=h2VReco.size();
    else return h2VReco[i];
    if (i>=h2VIso.size()) i-=h2VIso.size();
    else return h2VIso[i];
    if (i>=h2VHLT.size()) {
      std::cout << "h2fullList: i=" << i << " is too large\n";
      return NULL;
    }
    return h2VHLT[i];
  }

  TString fullListHistoName(unsigned int i) const
  { return TString(this->h2fullList(i)->GetName()); }
  TString fullListHistoTitle(unsigned int i) const
  { return TString(this->h2fullList(i)->GetTitle()); }

  void removeError();
  void setError(const DYTnPEff_t &e);
  void resetAll();
  int excludeGap();

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


  // scale factors
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

  // misc set
  double totalEffIdx_misc(int ibin1, int jbin1, int ibin2, int jbin2,
			  int mc, int miscFlag) const;
  double totalEff_misc(double eta1, double pt1, double eta2, double pt2,
		       int mc, int miscFlag) const;

  double totalEff_misc(const TLorentzVector *v1, const TLorentzVector *v2,
		       int mc, int miscFlag) const
  { return totalEff_misc(v1->Eta(),v1->Pt(),v2->Eta(),v2->Pt(),mc,miscFlag); }
  double scaleFactor_misc(double eta1, double pt1, double eta2, double pt2,
			  int miscFlag) const;
  double scaleFactor_misc(const TLorentzVector *v1,
			  const TLorentzVector *v2,
			  int miscFlag) const
  { return scaleFactor_misc(v1->Eta(),v1->Pt(),v2->Eta(),v2->Pt(),miscFlag); }
  double scaleFactorIdx_misc(int ibin1, int jbin1, int ibin2, int jbin2,
			     int miscFlag) const;


  int updateVectors(); // copy TH2D* pointers to arrays
  int assign(const DYTnPEff_t &e, TString tag);
  int assign_DiffAsUnc(const DYTnPEff_t &nominal, int relative); // histos have to be created
  int randomize(const DYTnPEff_t &e, TString tag, int systematic=0);
  int randomizeRelErr(const DYTnPEff_t &e0, const DYTnPEff_t &eForRelUnc,
		      TString tag, int systematic=0);
  int randomize(const DYTnPEff_t &e_errPos, const DYTnPEff_t &e_errNeg,
		TString tag, int systematic=0);

  int multiply(const DYTnPEff_t &e);
  int add(const DYTnPEff_t &e);
  int addShiftByUnc(const DYTnPEff_t &e, double nSigmas_data,double nSigmas_mc);

  // (x,dx) * (1,e0* dRelUnc)
  int includeRelUnc(const DYTnPEff_t &e0, const DYTnPEff_t &eRelUnc);
  //int includeUnc(const DYTnPEff_t &eUnc); // use add, instead

  int checkBinning() const;
  void printNumbers(TString fileTag="") const;
  void displayAll() const;
  void displayEffTot(int data_or_mc= 1+2, TString tag="", int hlt4p3=-1,
		     int misc=_misc_empty) const;
  void plotProfiles(int idxOnly=-1,
		    int skipGap=0, int plotVsPt=0, DYTnPEff_t *tnp2=NULL,
		    TString label1="", TString label2="",
		    DYTnPEff_t *tnp3=NULL, TString label3="");
  void plotSFProfiles(int idxOnly=-1,
		      int skipGap=0, int plotVsPt=0, DYTnPEff_t *tnp2=NULL,
		      TString label1="", TString label2="",
		      DYTnPEff_t *tnp3=NULL, TString label3="");
  void listNumbers() const;
  void listNumbersBrief(int nLinesX, int nLinesY) const;
  void printEffRatios(const DYTnPEff_t &e, int compareErrs=0,
		      double markIfDiffRelTol=0.) const;

  // internal save/load
  int save(TFile &fout, TString subdirTag="") const;
  int load(TFile &fin, TString subdir="DYTnPEff", TString tag="");
  int load(TString fname, TString subdir="DYTnPEff", TString tag="");

  // load from Kyeongpil's files
  int load_DYMM_TnP(TFile &fin, int version, TString tag="");
  int load_DYMM_TnP_asymmEff(TFile &fin, int version, int posErr, TString tag="");

};

// -------------------------------------------------------------

class DYTnPEffColl_t {
public:
  typedef enum { _rnd_ampl=0, _rnd_barrel_vs_endcap,
		 _rnd_low_vs_high_pt } TRandomizationKind_t;

protected:
  DYTnPEff_t fTnPEff; // efficiencies with statistical uncertainty
  std::vector<DYTnPEff_t*> fTnPEffSrcV;  // alternative efficiencies for systematics
  std::vector<DYTnPEff_t*> fTnPEffSystV; // derived systematic uncertainties
  int fElChannel;
public:
  DYTnPEffColl_t(int setlElChannel=0);
  // copy all or a subset of syst.uncertainties
  DYTnPEffColl_t(const DYTnPEffColl_t &e, TString tag, int iSrc=-1);

  const DYTnPEff_t& getTnPWithStatUnc() const { return fTnPEff; }
  DYTnPEff_t* getTnPWithSystUnc(int idx) const;
  DYTnPEff_t* getTnPWithTotUnc(TString tag="_totUnc") const;
  template<class idx_t>
  const DYTnPEff_t* getTnPSystUncSource(idx_t idx) const {return fTnPEffSrcV[idx];}
  template<class idx_t>
  const DYTnPEff_t* getTnPSystUnc(idx_t idx) const { return fTnPEffSystV[idx]; }
  DYTnPEff_t* getTnPSource(int srcIdx) const;
  DYTnPEff_t* getTnPShiftByUnc(TString tag, int srcIdx,
			       double nSigmas_data, double nSigmas_mc) const;

  int isElChannel() const { return fElChannel; }
  int ptrsOk() const;
  int assign(const DYTnPEffColl_t &coll, TString tag, int iSrc=-1);
  int excludeGap();

  double scaleFactor(const TLorentzVector *v1, const TLorentzVector *v2,
		     int hlt4p3) const
  { return fTnPEff.scaleFactor(v1,v2,hlt4p3); }

  double scaleFactorIdx(int ibin1, int jbin1, int ibin2, int jbin2, int hlt4p3) const
  { return fTnPEff.scaleFactorIdx(ibin1,jbin1,ibin2,jbin2,hlt4p3); }

  int assignStatErr(const DYTnPEff_t &e, TString tag);
  int addSystErrSource(const DYTnPEff_t &e, TString tag);

  // Produce randomized set of the scale factors
  // scrIdx=-1111 -- both stat and syst uncertainties
  // srcIdx= 0 -- only stat uncertainty
  // srcIdx= number -- syst uncertainty from src (number-1)
  // srcIdx= -number -- the central value of randomization would be 1, the uncertainty - relative
  // srcIdx= 9999 -- no randomization, just put it together
  DYTnPEff_t* randomize(int srcIdx, TString tag) const;

  // rndKind from TRandomizationKind_t
  DYTnPEff_t* randomizeByKind(int rndKind, int hlt4p3, TString tag,
			      int iSrcOnly=-1,
			      int maxSigmaData=0, int maxSigmaMC=0,
      TH2D **h2chk=NULL, unsigned int ihChk=0, unsigned int iSrcChk=0) const;

  //void printNumbers() const;
  void displayAll(int includeSrc=0) const;
  //void plotProfiles(int skipGap, int plotVsPt) const;
  void listNumbers(int includeSrc=0) const;

  int save(TFile &fout, TString subdirTag="") const;
  int save(TString fname, TString subdirTag="") const;
  int load(TFile &fin, TString subdir="DYTnPEffColl", TString tagList="");
  int load(TString fname, TString subdir="DYTnPEffColl", TString tagList="");
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
      if ((v1->Pt()>10) && (v2->Pt()>10)) {
	std::cout << "negative fi in " << v1->Pt() << ", " << v1->Eta()
		  << " : " << fi1 << "; " << v2->Pt() << ", "
		  << v2->Eta() << " : " << fi2 << "\n";
      }
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
  TH1D* calculateScaleFactor_misc(const DYTnPEff_t &eff, int miscFlag,
				  TString hName, TString hTitle,
				  int followFIBin1=0, int followFIBin2=0,
				  TH1D *h1followedBin=NULL) const;

  void displayAll() const;
  void listAll() const;
  std::vector<TH1D*> avgAxisValues(std::vector<TH2D*> *h2ptSpace=NULL,
				   std::vector<TH2D*> *h2etaSpace=NULL) const;

  int save(TFile &fout);
  int load(TFile &fin, TString subdir);
  int load(TString fname, TString subdir);

 protected:
  int prepareVector();

};

// -------------------------------------------------------------
// -------------------------------------------------------------

TString etaRangeStr(const TH2D *h2, int iEta, int useAbsEta, int *isGap=NULL);
TString ptRangeStr(const TH2D *h2, int iPt);

void plotEffs(const std::vector<TH2D*> &hV,
	      const std::vector<HistoStyle_t> &hsV,
	      const std::vector<TString> &labels,
	      TString plotNameBase, TString effName,
	      int allInOne,
	      int isEff,
	      int skipGap, int plotVsPt);
void setAutoRanges(const std::vector<TH1D*> &h1V,
		   double &range_min, double &range_max, int isEffRange,
		   int silent);

// -------------------------------------------------------------

#endif
