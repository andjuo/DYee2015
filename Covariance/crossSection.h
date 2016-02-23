#ifndef crossSection_H
#define crossSection_H

#include "inputs.h"
#include "../RooUnfold/src/RooUnfoldResponse.h"
#include "../RooUnfold/src/RooUnfoldBayes.h"

// -----------------------------------------------------------

typedef enum { _varNone=0, _varYield, _varBkg, _varSig, _varDetRes, _varFSRRes,
    _varEff, _varRho, _varAcc, _varEffAcc, _varLast } TVaried_t;

typedef enum { _bkgZZ=0, _bkgWZ, _bkgWW, _bkgTTbar, _bkgDYtautau, _bkgTW,
	       _bkgWJets, _bkgQCD, _bkgLast } TBkg_t;

typedef enum { _csPostFsrInAcc=0, _csPostFsrFullSp,
	       _csPreFsrInAcc, _csPreFsrFullSp } TCSType_t;

TString variedVarName(TVaried_t);
TString bkgName(TBkg_t);
TString csTypeName(TCSType_t);

inline
TBkg_t next(TBkg_t &b)
{ if (b<_bkgLast) b=TBkg_t(int(b)+1); return b; }

// -----------------------------------------------------------

RooUnfoldResponse* loadRooUnfoldResponse(TString fname, TString fieldName, TString name);
RooUnfoldResponse* loadRooUnfoldResponse(TFile &fin, TString fieldName, TString name);
RooUnfoldBayes* loadRooUnfoldBayes(TFile &fin, TString fieldName, TString name);
void plotHisto(RooUnfoldResponse &rs, TString cNameBase, int logx=0, int logy=0);

// -----------------------------------------------------------

class MuonCrossSection_t;

// -----------------------------------------------------------

class CrossSection_t {
  friend class MuonCrossSection_t;
 protected:
  TString fName, fTag;
  TH1D *fh1Yield, *fh1Bkg;
  RooUnfoldResponse *fDetRes, *fFSRRes;
  TH1D *fh1Eff, *fh1Rho, *fh1Acc;
  TH1D* fh1EffAcc;
  TH1D* fh1Theory;
  TVaried_t fVar;
  TCSType_t fCSType;
  int fNIters;
  double fLumi;
 protected:
  TH1D *fh1Signal, *fh1Unf, *fh1UnfRhoCorr, *fh1UnfRhoEffCorr;
  TH1D *fh1UnfRhoEffAccCorr, *fh1PreFsr, *fh1PreFsrCS;
  TH1D *fh1Varied;
  RooUnfoldBayes *fDetResBayes, *fFSRBayes;
 public:
  CrossSection_t(TString setName="xs", TString setTag="",
		 TCSType_t setCSType=_csPreFsrFullSp);
  CrossSection_t(const CrossSection_t &cs,
		 TString setName="xsNew", TString setTag="newTag");
  void clearAllPtrs();

  const TH1D *h1Yield() const { return fh1Yield; }
  void h1Yield(const TH1D* h1) { fh1Yield=copy(h1,"h1Yield",fTag); }
  const TH1D *h1Bkg() const { return fh1Bkg; }
  void h1Bkg(const TH1D* h1, double scale=1.) {
    fh1Bkg=copy(h1,"h1Bkg",fTag);
    if (scale!=1.) fh1Bkg->Scale(scale);
  }
  const TH1D *h1Eff() const { return fh1Eff; }
  void h1Eff(const TH1D *h1) { fh1Eff=copy(h1,"h1Eff",fTag); }
  const TH1D *h1Rho() const { return fh1Rho; }
  void h1Rho(const TH1D *h1) { fh1Rho=copy(h1,"h1Rho",fTag); }
  const TH1D *h1Acc() const { return fh1Acc; }
  void h1Acc(const TH1D *h1) { fh1Acc=copy(h1,"h1Acc",fTag); }
  const TH1D *h1EffAcc() const { return fh1EffAcc; }
  void h1EffAcc(const TH1D *h1) { fh1EffAcc=copy(h1,"h1EffAcc",fTag); }
  const RooUnfoldResponse* detRes() const { return fDetRes; }
  void detRes(const RooUnfoldResponse &rs) { fDetRes= new RooUnfoldResponse(rs); }
  const RooUnfoldResponse* fsrRes() const { return fFSRRes; }
  void fsrRes(const RooUnfoldResponse &rs) { fFSRRes= new RooUnfoldResponse(rs); }
  TString name() const { return fName; }
  void name(TString new_name) { fName=new_name; }
  TString tag() const { return fTag; }
  void tag(TString new_tag) { fTag=new_tag; }
  TVaried_t var() const { return fVar; }
  void var(TVaried_t new_var) { fVar=new_var; }
  TCSType_t csType() const { return fCSType; }
  void csType(TCSType_t setCS) { fCSType=setCS; }
  int nIters() const { return fNIters; }
  void nIters(int setIters) { fNIters=setIters; }
  double lumi() const { return fLumi; }
  void lumi(double set_lumi) { fLumi=set_lumi; }

  const TH1D* h1Theory() const { return fh1Theory; }
  void h1Theory(const TH1D *h1) { fh1Theory=copy(h1,"h1Theory",fTag); }

  const TH1D* h1Signal() const { return fh1Signal; }
  const TH1D* h1Unf() const { return fh1Unf; }
  const TH1D* h1UnfRhoCorr() const { return fh1UnfRhoCorr; }
  const TH1D* h1UnfRhoEffCorr() const { return fh1UnfRhoEffCorr; }
  const TH1D* h1UnfRhoEffAccCorr() const { return fh1UnfRhoEffAccCorr; }
  const TH1D* h1PreFsr() const { return fh1PreFsr; }
  const TH1D* h1PreFsrCS() const { return fh1PreFsrCS; }

  TH1D* getVariedHisto(TVaried_t new_var);

  int checkPtrs(int *failCode=NULL) const;
  int assign(const CrossSection_t &cs);

  TH1D* calcCrossSection();
  TH1D* calcCrossSection(TVaried_t new_var, int idx);
  TCanvas* plotCrossSection(TString canvName="cs");

  int sampleRndVec(TVaried_t new_var, int sampleSize,
		   std::vector<TH1D*> &rndCS);
  int sampleRndVec(TVaried_t new_var, const std::vector<TH1D*> &rndHistos,
		   std::vector<TH1D*> &rndCS);
  int deriveCov(const std::vector<TH1D*> &rndCS, TH1D **h1avgCS, TH2D **h2cov);

  int save(TString fname) const;
  int load(TString fname, TString setTag="");

  // ----------------------------------------------

  TH1D* copy(const TH1D *h1, TString newName, TString useTag) const;
  TH1D* copy(const TH1D *h1, TString useTag) const;
  TH1D* copy(const TH1D *h1) const { return copy(h1,fTag); }
  TH1D* copy(const TH1 *h1unf, TString newName, TString newTitle, const TH1D *h1protoType) const;

  // ----------------------------------------------

};

// -----------------------------------------------------------

class MuonCrossSection_t {
 protected:
  CrossSection_t fCSa, fCSb;
  TString fTag;
  std::vector<TH1D*> fh1BkgV;
  TVectorD fBkgWeights; //  xsec/nEvents
  TH1D *fh1Bkg;
  TH1D *fh1CS;
  TH1D *fh1Theory;
 protected:
  RooUnfoldBayes *fFSRBayes;
 public:
  MuonCrossSection_t(TString setName, TString setTag,
		     double lumiA=0, double lumiB=0);

  void clear() { fh1BkgV.clear(); }

  const CrossSection_t& csA() const { return fCSa; }
  CrossSection_t& editCSa() { return fCSa; }
  const CrossSection_t& csB() const { return fCSb; }
  CrossSection_t& editCSb() { return fCSb; }

  TString tag() const { return fTag; }
  void lumi(double lumiA, double lumiB) { fCSa.lumi(lumiA); fCSb.lumi(lumiB); }

  const TH1D* h1Bkg() const { return fh1Bkg; }
  void h1Bkg(const TH1D *h1, int clearVec=1);
  const std::vector<TH1D*>& bkgV() const { return fh1BkgV; }
  const TVectorD& bkgW() const { return fBkgWeights; }
  void setBkgV(const std::vector<TH1D*> &setBkg, const std::vector<double> &weights);
  int recalcBkg(const std::vector<double> &weights);

  const TH1D* h1CS() const { return fh1CS; }
  const TH1D* h1Theory() const { return fh1Theory; }
  void h1Theory(const TH1D* h1) { fh1Theory=fCSa.copy(h1,"h1Theory",fTag); }

  TH1D* calcCrossSection();
  TCanvas* plotCrossSection(TString canvName="cs", int recalculate=0);

  int sampleRndVec(TVaried_t new_var, int sampleSize,
		   std::vector<TH1D*> &rndCS,
		   std::vector<TH1D*> *rndCSa_out=NULL,
		   std::vector<TH1D*> *rndCSb_out=NULL,
		   std::vector<TH1D*> *rndVarVec_out=NULL);
  int deriveCov(const std::vector<TH1D*> &rndCS, TVaried_t var,
		TH1D **h1avgCS_out, TH2D **h2cov_out);

  int save(TString fnameBase);
  int load(TString fnameBase, TString setTag);

 protected:
  TH1D* calcPreFsrCS_sumAB(const TH1D *h1a, const TH1D *h1b,TString useTag);
};

// -----------------------------------------------------------


#endif
