#ifndef crossSection_H
#define crossSection_H

#include "inputs.h"
#include "../RooUnfold/src/RooUnfoldResponse.h"
#include "../RooUnfold/src/RooUnfoldBayes.h"
#include "../RooUnfold/src/RooUnfoldInvert.h"
#include "DYbinning.h"

// -----------------------------------------------------------

typedef enum { _varNone=0, _varYield, _varBkg, _varBkgXS, _varSig, _varDetRes,
	       _varFSRRes, _varFSRRes_Poisson,
	       _varEff, _varRho, _varRhoFile,
	       _varAcc, _varEffAcc, _varLast,
	       _varRhoSyst, _varTheory, _varYieldPoisson } TVaried_t;

typedef enum { _bkgZZ=0, _bkgWZ, _bkgWW, _bkgTTbar, _bkgDYtautau, _bkgTW,
	       _bkgWJets, _bkgQCD, _bkgLast } TBkg_t;

typedef enum { _csPostFsrInAcc=0, _csPostFsrFullSp,
	       _csPreFsrInAcc, _csPreFsrFullSp } TCSType_t;

typedef enum { _tnpSyst_none=0,
	       _tnpSyst_bkgPdf, _tnpSyst_sigPdf, _tnpSyst_NLOvsLO,
	       _tnpSyst_tagDef, _tnpSyst_last } TTnPSystType_t;

typedef std::map<TVaried_t,TString> finfomap_t;
typedef std::map<TVaried_t,double> fvarweightmap_t;

TString variedVarName(TVaried_t);
TString bkgName(TBkg_t);
TString csTypeName(TCSType_t);
TString tnpSystName(TTnPSystType_t, int shortName=0);

inline
TBkg_t next(TBkg_t &b)
{ if (b<_bkgLast) b=TBkg_t(int(b)+1); return b; }

inline
TTnPSystType_t next(TTnPSystType_t &s)
{ if (s<_tnpSyst_last) s=TTnPSystType_t(int(s)+1); return s; }

// -----------------------------------------------------------

template<class type_t>
int findVaried(const std::map<TVaried_t,type_t> &aMap, TVaried_t var,
	       int verbose=0)
{
  int ii=0, idx=-1;
  for (finfomap_t::const_iterator it= aMap.begin();
       it!=aMap.end(); it++, ii++) {
    //std::cout << "chk ii=" << ii << " " << variedVarName(it->first) << "\n";
    if (it->first == var) {
      idx=ii;
      break;
    }
  }
  if ((idx==-1) && verbose) {
    std::cout << "failed to find " << variedVarName(var) << " in the map\n";
  }
  return idx;
}

// -----------------------------------------------------------

RooUnfoldResponse* loadRooUnfoldResponse(TString fname, TString fieldName,
					 TString name, int warnIfAbsent=1);
RooUnfoldResponse* loadRooUnfoldResponse(TFile &fin, TString fieldName,
					 TString name, int warnIfAbsent=1);
RooUnfoldBayes* loadRooUnfoldBayes(TFile &fin, TString fieldName, TString name);
void plotHisto(RooUnfoldResponse &rs, TString cNameBase, int logx=0, int logy=0);
RooUnfoldResponse* randomizeWithinErr(const RooUnfoldResponse *R, TString name,
				      int poissonRnd=0);

// -----------------------------------------------------------

class MuonCrossSection_t;

// -----------------------------------------------------------

struct RndVecInfo_t {
  TString fFName, fHistoNameBase, fHistoNameBaseV2;
  int fMaxSamples;
public:
  RndVecInfo_t(TString setFName="",
	       TString setHistoNameBase="",
	       TString setHistoNameBaseV2="",
	       int nMaxSamples=-1) :
      fFName(setFName),
      fHistoNameBase(setHistoNameBase),
      fHistoNameBaseV2(setHistoNameBaseV2),
      fMaxSamples(nMaxSamples)

  {}
  TString fname() const { return fFName; }
  TString histoNameBase() const { return fHistoNameBase; }
  TString histoNameBaseV2() const { return fHistoNameBaseV2; }
  int maxSamples() const { return fMaxSamples; }

  void assign(TString setFName,
	      TString setHistoNameBase, TString setHistoNameBaseV2)
  {
    fFName=setFName;
    fHistoNameBase=setHistoNameBase;
    fHistoNameBaseV2=setHistoNameBaseV2;
  }

  TString histoName(int iVar) const
  { return fHistoNameBase + Form("%d",iVar); }
  TString histoNameV2(int iVar) const
  { return fHistoNameBaseV2 + Form("%d",iVar); }

  int loadHistos(int sampleSize,
		 std::vector<TH1D*> &rndVec,
		 std::vector<TH1D*> *rndVec2=NULL) const;
};

// -----------------------------------------------------------

class CrossSection_t {
  friend class MuonCrossSection_t;
 protected:
  TString fName, fTag;
  TVersion_t fVersion;
  TH1D *fh1Yield, *fh1Bkg;
  RooUnfoldResponse *fDetRes, *fFSRRes;
  TH1D *fh1Eff, *fh1Rho, *fh1Acc;
  TH1D* fh1EffAcc;
  TH1D* fh1Theory;
  TVaried_t fVar;
  TCSType_t fCSType;
  int fNItersDetRes, fNItersFSR;
  double fLumi;
 protected:
  TH1D *fh1Signal, *fh1Unf, *fh1UnfRhoCorr, *fh1UnfRhoEffCorr;
  TH1D *fh1UnfRhoEffAccCorr, *fh1PreFsr, *fh1PreFsrCS;
  TH1D *fh1Varied;
  RooUnfoldResponse *fResVaried;
  RooUnfoldBayes *fDetResBayes, *fFSRBayes;
 public:
  CrossSection_t(TString setName="xs", TString setTag="",
		 TCSType_t setCSType=_csPreFsrFullSp,
		 TVersion_t setVersion=_verUndef);
  CrossSection_t(const CrossSection_t &cs,
		 TString setName="xsNew", TString setTag="newTag",
		 TVersion_t newVersion=_verUndef);
  void clearAllPtrs();

  const TH1D *h1Yield() const { return fh1Yield; }
  void h1Yield(const TH1D* h1) { fh1Yield=copy(h1,"h1Yield",fTag); }
  const TH1D *h1Bkg() const { return fh1Bkg; }
  TH1D *editH1Bkg() { return fh1Bkg; }
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
  void removeEffAcc() { if (fh1EffAcc) { delete fh1EffAcc; fh1EffAcc=NULL; } }
  const RooUnfoldResponse* detRes() const { return fDetRes; }
  RooUnfoldResponse* detResPtr() { return fDetRes; }
  void detRes(const RooUnfoldResponse &rs) { fDetRes= new RooUnfoldResponse(rs); }
  const RooUnfoldResponse* fsrRes() const { return fFSRRes; }
  RooUnfoldResponse* fsrResPtr() { return fFSRRes; }
  void fsrRes(const RooUnfoldResponse &rs) { fFSRRes= new RooUnfoldResponse(rs); }
  TString name() const { return fName; }
  void name(TString new_name) { fName=new_name; }
  TString tag() const { return fTag; }
  void tag(TString new_tag) { fTag=new_tag; }
  TVersion_t version() const { return fVersion; }
  TString versionAsString() const { return versionName(fVersion); }
  void version(TVersion_t setVersion) { fVersion=setVersion; }
  TVaried_t var() const { return fVar; }
  void var(TVaried_t new_var) { fVar=new_var; }
  TCSType_t csType() const { return fCSType; }
  void csType(TCSType_t setCS) { fCSType=setCS; }
  void setNIters_internal(TVersion_t);
  int nItersDetRes() const { return fNItersDetRes; }
  void nItersDetRes(int setIters) { fNItersDetRes=setIters; }
  int nItersFSR() const { return fNItersFSR; }
  void nItersFSR(int setIters) { fNItersFSR=setIters; }
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
  void removeBkg(); // set bkg to 0
  void removeRho(); // set rho factor to 1

  TH1D* calcCrossSection(int removeNegativeSignal=0);
  TH1D* calcCrossSection(TVaried_t new_var, int idx, int removeNegativeSignal=0);
  TH1D* calcCrossSection(TVaried_t new_var, int idx, const TH1D* set_h1Var,
			 int removeNegativeSignal=0);
  TCanvas* plotCrossSection(TString canvName="cs", int removeNegativeSignal=0);
  TCanvas* plotCrossSection_StepByStep(TString canvName="csStep");

  int sampleRndVec(TVaried_t new_var, int sampleSize, int removeNegativeSignal,
		   std::vector<TH1D*> &rndCS);
  int sampleRndVec(TVaried_t new_var,
		   const std::vector<TH1D*> &rndHistos,
		   int removeNegativeSignal,
		   std::vector<TH1D*> &rndCS);
  int sampleRndVec(TVaried_t new_var, int sampleSize,
		   const RndVecInfo_t &info,
		   int removeNegativeSignal,
		   std::vector<TH1D*> &rndCS);
  int sampleRndResponse(TVaried_t new_var, int sampleSize,
			int removeNegativeSignal,
			std::vector<TH1D*> &rndCS,
			std::vector<RooUnfoldResponse*> *rndRespVPtr);
  int sampleRndRespVec(TVaried_t new_var,
		       const std::vector<RooUnfoldResponse*> &rndRespV,
		       int removeNegativeSignal,
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
  TVersion_t fVersion;
  std::vector<TH1D*> fh1BkgV;
  TVectorD fBkgWeightUnc; // bkg weight uncertainty
  TH1D *fh1Bkg;
  TH1D *fh1PostFSR, *fh1PreFSR;
  TH1D *fh1CS;
  TH1D *fh1Theory;
 protected:
  RooUnfoldBayes *fFSRBayes;
 public:
  MuonCrossSection_t(TString setName, TString setTag,
		     double lumiA=0, double lumiB=0,
		     TVersion_t setVersion=_verMu1);

  void clear() { fh1BkgV.clear(); }

  const CrossSection_t& csA() const { return fCSa; }
  CrossSection_t& editCSa() { return fCSa; }
  const CrossSection_t& csB() const { return fCSb; }
  CrossSection_t& editCSb() { return fCSb; }

  TString tag() const { return fTag; }
  TVersion_t version() const { return fVersion; }
  TString versionAsString() const { return versionName(fVersion); }
  void version(TVersion_t setVersion) { fVersion=setVersion; }
  void lumi(double lumiA, double lumiB) { fCSa.lumi(lumiA); fCSb.lumi(lumiB); }
  double lumiTot() const { return (fCSa.lumi() + fCSb.lumi()); }

  const TH1D* h1Bkg() const { return fh1Bkg; }
  // The latest version of muon results used backgrounds that were partially
  // scaled by a SF.
  // Instead of old h1Bkg(), call recalcBkg()
  void h1Bkg_approximate(const TH1D *h1, int clearVec=1);
  const std::vector<TH1D*>& bkgV() const { return fh1BkgV; }
  const TVectorD& bkgWUnc() const { return fBkgWeightUnc; }
  int setBkgV(const std::vector<TH1D*> &setBkg, const std::vector<double> &set_bkgWUncertainty);
  int recalcBkg(const std::vector<double> *weights,
		TH1D **h1bkg_fromMC=NULL, TH1D **h1bkg_fromData=NULL);

  const TH1D *h1PostFsr() const { return fh1PostFSR; }
  const TH1D *h1PreFsr() const { return fh1PreFSR; }
  const TH1D* h1CS() const { return fh1CS; }
  const TH1D* h1Theory() const { return fh1Theory; }
  void h1Theory(const TH1D* h1) { fh1Theory=fCSa.copy(h1,"h1Theory",fTag); }

  TH1D* calcCrossSection(int removeNegativeSignal=0);
  TH1D* calcCrossSection(TVaried_t new_var, int idx,
			 const TH1D *h1_setVarA, const TH1D *h1_setVarB,
			 int removeNegativeSignal=0);
  TCanvas* plotCrossSection(TString canvName="cs", int recalculate=0,
			    int removeNegativeSignal=0);

  int sampleRndVec(TVaried_t new_var, int sampleSize,
		   int removeNegativeSignal,
		   std::vector<TH1D*> &rndCS,
		   std::vector<TH1D*> *rndCSa_out=NULL,
		   std::vector<TH1D*> *rndCSb_out=NULL,
		   std::vector<TH1D*> *rndVarVec_out=NULL);

  int sampleRndVec(TVaried_t new_var, int sampleSize,
		   const RndVecInfo_t &info,
		   int removeNegativeSignal,
		   std::vector<TH1D*> &rndCS,
		   std::vector<TH1D*> *rndCSa_out=NULL,
		   std::vector<TH1D*> *rndCSb_out=NULL,
		   std::vector<TH1D*> *rndVarVec1_out=NULL,
		   std::vector<TH1D*> *rndVarVec2_out=NULL);

  int deriveCov(const std::vector<TH1D*> &rndCS, TVaried_t var,
		TH1D **h1avgCS_out, TH2D **h2cov_out);

  int save(TString fnameBase);
  int load(TString fnameBase, TString setTag);

 protected:
  TH1D* calcPreFsrCS_sumAB(const TH1D *h1a, const TH1D *h1b,TString useTag);
};

// -----------------------------------------------------------

namespace DYtools {

#ifdef def_41massBin
  const TVersion_t activeVer= _verEl3mb41;
#else
#  ifdef def_42massBin
  const TVersion_t activeVer= _verEl3mb42;
#  else
  //const TVersion_t activeVer= _verEl3;
  const TVersion_t activeVer= _verElMay2017;
#  endif
#endif


  inline
  TH1D *rebin_43toLess_perMBW(const TH1D *h1_inp_perMBW, int printMBW=0) {
    TH1D* h1_noMBW= ::timesMassBinWidth(h1_inp_perMBW,printMBW);
    TH1D* h1_4X= rebin_43toLess(h1_noMBW);
    TH1D* h1_4XperMBW= ::perMassBinWidth(h1_4X,printMBW);
    delete h1_noMBW;
    delete h1_4X;
    return h1_4XperMBW;
  }
};

// -----------------------------------------------------------


#endif
