#include "crossSection.h"
#include <TKey.h>
//#include <TList.h>

// --------------------------------------------------------------
// --------------------------------------------------------------

TString variedVarName(TVaried_t v)
{
  TString name="variedUNKNOWN";
  switch(v) {
  case _varNone: name="none"; break;
  case _varYield: name="varYield"; break;
  case _varBkg: name="varBkg"; break;
  case _varSig: name="varSig"; break;
  case _varDetRes: name="varDetRes"; break;
  case _varFSRRes: name="varFSRRes"; break;
  case _varEff: name="varEff"; break;
  case _varRho: name="varRho"; break;
  case _varAcc: name="varAcc"; break;
  case _varEffAcc: name="varEffAcc"; break;
  case _varLast: name="varLast"; break;
  default:
    std::cout << "variedVarName is not ready for this var\n";
  }
  return name;
}

// --------------------------------------------------------------

TString bkgName(TBkg_t b)
{
  TString name="bkgUNKNOWN";
  switch(b) {
  case _bkgZZ: name="ZZ"; break;
  case _bkgWZ: name="WZ"; break;
  case _bkgWW: name="WW"; break;
  case _bkgTTbar: name="ttbar"; break;
  case _bkgDYtautau: name="DYtautau"; break;
  case _bkgTW: name="tW"; break;
  case _bkgWJets: name="WJets"; break;
  case _bkgQCD: name="QCD"; break;
  default:
    std::cout << "bkgName is not ready for this bkg\n";
  }
  return name;
}

// --------------------------------------------------------------

TString csTypeName(TCSType_t cs)
{
  TString name="csTypeUNKNOWN";
  switch(cs) {
  case _csPostFsrInAcc: name="postFSRInAcc"; break;
  case _csPostFsrFullSp: name="postFSRFullSp"; break;
  case _csPreFsrInAcc: name="preFSRInAcc"; break;
  case _csPreFsrFullSp: name="preFSRFullSp"; break;
  default:
    std::cout << "csTypeName is not ready for this cs type\n";
  }
  return name;
}

// --------------------------------------------------------------
// --------------------------------------------------------------

RooUnfoldResponse *loadRooUnfoldResponse(TString fname, TString fieldName,
					 TString setName)
{
  TFile fin(fname);
  RooUnfoldResponse *rs= loadRooUnfoldResponse(fin,fieldName,setName);
  fin.Close();
  return rs;
}

// --------------------------------------------------------------

RooUnfoldResponse *loadRooUnfoldResponse(TFile &fin, TString fieldName,
					 TString setName)
{
  if (!fin.IsOpen()) return NULL;
  TObject *chk= fin.Get(fieldName);
  if (!chk) return NULL;
  delete chk;
  RooUnfoldResponse *rs= new RooUnfoldResponse();
  rs->Read(fieldName);
  rs->SetName(setName);
  return rs;
}

// --------------------------------------------------------------

RooUnfoldBayes *loadRooUnfoldBayes(TFile &fin, TString fieldName,
				   TString setName)
{
  if (!fin.IsOpen()) return NULL;
  TObject *chk= fin.Get(fieldName);
  if (!chk) return NULL;
  delete chk;
  RooUnfoldBayes *rs= new RooUnfoldBayes();
  rs->Read(fieldName);
  rs->SetName(setName);
  return rs;
}

// --------------------------------------------------------------

void plotHisto(RooUnfoldResponse &rs, TString cNameBase, int logx, int logy)
{
  plotHisto(rs.Hresponse(),cNameBase + "_Response",logx,logx); // both logx!
  plotHisto(rs.Hmeasured(),cNameBase + "_Meas",logx,logy);
  plotHisto(rs.Hfakes(), cNameBase + "_Fakes",logx,logy);
  plotHisto(rs.Htruth(), cNameBase + "_Truth",logx,logy);
}

// --------------------------------------------------------------
// --------------------------------------------------------------

CrossSection_t::CrossSection_t(TString setName, TString setTag,
			       TCSType_t set_csType) :
  fName(setName), fTag(setTag),
  fh1Yield(NULL), fh1Bkg(NULL),
  fDetRes(NULL), fFSRRes(NULL),
  fh1Eff(NULL), fh1Rho(NULL), fh1Acc(NULL),
  fh1EffAcc(NULL),
  fh1Theory(NULL),
  fVar(_varNone), fCSType(set_csType),
  fNIters(4), fLumi(1.),
  fh1Signal(NULL), fh1Unf(NULL), fh1UnfRhoCorr(NULL),
  fh1UnfRhoEffCorr(NULL), fh1UnfRhoEffAccCorr(NULL),
  fh1PreFsr(NULL), fh1PreFsrCS(NULL),
  fh1Varied(NULL), fDetResBayes(NULL), fFSRBayes(NULL)
{
}

// --------------------------------------------------------------

CrossSection_t::CrossSection_t(const CrossSection_t &cs,
			       TString setName, TString setTag) :
  fName(setName), fTag(setTag),
  fh1Yield(NULL), fh1Bkg(NULL),
  fDetRes(NULL), fFSRRes(NULL),
  fh1Eff(NULL), fh1Rho(NULL), fh1Acc(NULL),
  fh1EffAcc(NULL),
  fh1Theory(NULL),
  fVar(_varNone), fCSType(cs.fCSType),
  fNIters(4), fLumi(1.),
  fh1Signal(NULL), fh1Unf(NULL), fh1UnfRhoCorr(NULL),
  fh1UnfRhoEffCorr(NULL), fh1UnfRhoEffAccCorr(NULL),
  fh1PreFsr(NULL), fh1PreFsrCS(NULL),
  fh1Varied(NULL), fDetResBayes(NULL), fFSRBayes(NULL)
{
  if (!this->assign(cs)) {
    std::cout << "error in CrossSection_t::CrossSection_t(CrossSection_t)\n";
  }
}

// --------------------------------------------------------------

void CrossSection_t::clearAllPtrs()
{
  if (fh1Yield) { delete fh1Yield; fh1Yield=NULL; }
  if (fh1Bkg) { delete fh1Bkg; fh1Bkg=NULL; }
  if (fDetRes) { delete fDetRes; fDetRes=NULL; }
  if (fFSRRes) { delete fFSRRes; fFSRRes=NULL; }
  if (fh1Eff) { delete fh1Eff; fh1Eff=NULL; }
  if (fh1Rho) { delete fh1Rho; fh1Rho=NULL; }
  if (fh1Acc) { delete fh1Acc; fh1Acc=NULL; }
  if (fh1EffAcc) { delete fh1EffAcc; fh1EffAcc=NULL; }
  if (fh1Theory) { delete fh1Theory; fh1Theory=NULL; }
  if (fh1Signal) { delete fh1Signal; fh1Signal=NULL; }
  if (fh1Unf) { delete fh1Unf; fh1Unf=NULL; }
  if (fh1UnfRhoCorr) { delete fh1UnfRhoCorr; fh1UnfRhoCorr=NULL; }
  if (fh1UnfRhoEffCorr) { delete fh1UnfRhoEffCorr; fh1UnfRhoEffCorr=NULL; }
  if (fh1UnfRhoEffAccCorr) { delete fh1UnfRhoEffAccCorr; fh1UnfRhoEffAccCorr=NULL; }
  if (fh1PreFsr) { delete fh1PreFsr; fh1PreFsr=NULL; }
  if (fh1PreFsrCS) { delete fh1PreFsrCS; fh1PreFsrCS=NULL; }
  if (fh1Varied) { delete fh1Varied; fh1Varied=NULL; }
  if (fDetResBayes) { delete fDetResBayes; fDetResBayes=NULL; }
  if (fFSRBayes) { delete fFSRBayes; fFSRBayes=NULL; }
}

// --------------------------------------------------------------

TH1D* CrossSection_t::getVariedHisto(TVaried_t new_var)
{
  fVar= new_var;
  TH1D* h1=NULL;
  switch(fVar) {
  case _varYield: h1=fh1Yield; break;
  case _varBkg: h1=fh1Bkg; break;
  case _varSig: h1=fh1Signal; break;
  case _varDetRes:
  case _varFSRRes:
    std::cout << "getVariedHisto should not be called for "
	      << variedVarName(new_var) << "\n";
    break;
  case _varEff: h1=fh1Eff; break;
  case _varRho: h1=fh1Rho; break;
  case _varAcc: h1=fh1Acc; break;
  case _varEffAcc: h1=fh1EffAcc; break;
  default:
    std::cout << "getVariedHisto should not be called for "
	      << variedVarName(new_var) << " (default warning)\n";
  }
  return h1;
}

// --------------------------------------------------------------

int CrossSection_t::checkPtrs(int *failCode) const {
  int fail=0;
  if (!fh1Yield || !fh1Bkg) fail=1;
  else if (!fDetRes || !fFSRRes) fail=2;
  else if (!fh1Rho) fail=3;
  else if ((!fh1Eff || !fh1Acc) && !fh1EffAcc) fail=4;
  else if (!fFSRRes &&
	   ((fCSType==_csPreFsrInAcc) || (fCSType==_csPreFsrFullSp))) fail=5;
  else if (!fh1Signal || !fh1Unf || !fh1UnfRhoCorr) fail=6;
  else if (!fh1UnfRhoEffCorr && !fh1UnfRhoEffAccCorr) fail=7;
  //else if (!fh1PreFsr) fail=8; // allow failure of this histo
  else if (!fh1Varied) fail=9;
  if (fail) std::cout << "CrossSection_t::checkPtrs fail=" << fail << "\n";
  if (failCode) *failCode= fail;
  return (!fail) ? 1:0;
}

// --------------------------------------------------------------

int CrossSection_t::assign(const CrossSection_t &cs)
{
  if (&cs == this) return 1;

  clearAllPtrs();

  fh1Yield= copy(cs.fh1Yield);
  fh1Bkg= copy(cs.fh1Bkg);
  fDetRes= new RooUnfoldResponse(*cs.fDetRes);
  fFSRRes= new RooUnfoldResponse(*cs.fFSRRes);
  fh1Eff= copy(cs.fh1Eff);
  fh1Rho= copy(cs.fh1Rho);
  fh1Acc= copy(cs.fh1Acc);
  fh1EffAcc= copy(cs.fh1EffAcc);
  fh1Theory= copy(cs.fh1Theory);
  fVar= cs.fVar;
  fNIters= cs.fNIters;
  fLumi= cs.fLumi;
  fh1Signal= copy(cs.fh1Signal);
  fh1Unf= copy(cs.fh1Unf);
  fh1UnfRhoCorr= copy(cs.fh1UnfRhoCorr);
  fh1UnfRhoEffCorr= copy(cs.fh1UnfRhoEffCorr);
  fh1UnfRhoEffAccCorr= copy(cs.fh1UnfRhoEffAccCorr);
  fh1Varied= copy(cs.fh1Varied);

  int res= checkPtrs();
  std::cout << "assign res=" << res << "\n";
  return res;
}

// --------------------------------------------------------------

TH1D* CrossSection_t::calcCrossSection()
{
  if (fh1PreFsr) { delete fh1PreFsr; fh1PreFsr=NULL; }
  if (fh1PreFsrCS) { delete fh1PreFsrCS; fh1PreFsrCS=NULL; }
  int code=0;
  int res= checkPtrs(&code);
  if (!res && (code<6)) return NULL;

  fh1Signal= copy(fh1Yield,"h1Signal",fTag);
  fh1Signal->Add(fh1Bkg,-1);
  if (fDetResBayes) delete fDetResBayes;
  fDetResBayes= new RooUnfoldBayes( fDetRes, fh1Signal, fNIters, false );
  fh1Unf= copy(fDetResBayes->Hreco(), "h1Unf_" + fTag,
	       "h1Unf_" + fTag
	       + TString(";") + fh1Yield->GetXaxis()->GetTitle()
	       + TString(";unfolded yield")
	       ,
	       fh1Yield); // set x binning
  fh1UnfRhoCorr= copy(fh1Unf,"h1UnfRho",fTag);
  fh1UnfRhoCorr->Divide(fh1Rho);
  if (!fh1EffAcc) {
    fh1UnfRhoEffCorr= copy(fh1UnfRhoCorr,"h1UnfRhoEff",fTag);
    fh1UnfRhoEffCorr->Divide(fh1Eff);
    fh1UnfRhoEffAccCorr= copy(fh1UnfRhoEffCorr,"h1UnfRhoEffAcc",fTag);
    fh1UnfRhoEffAccCorr->Divide(fh1Acc);
  }
  else {
    if (fh1UnfRhoEffCorr) { delete fh1UnfRhoEffCorr; fh1UnfRhoEffCorr=NULL; }
    fh1UnfRhoEffAccCorr= copy(fh1UnfRhoCorr,"h1UnfRhoEffAcc",fTag);
    fh1UnfRhoEffAccCorr->Divide(fh1EffAcc);
  }
  fh1UnfRhoEffAccCorr->Scale(1/fLumi);
  TH1D* h1Out=NULL;
  if (fCSType==_csPostFsrFullSp) {
    h1Out= fh1UnfRhoEffAccCorr;
  }
  else if (fCSType==_csPreFsrFullSp) {
    if (fFSRBayes) delete fFSRBayes;
    fFSRBayes= new RooUnfoldBayes( fFSRRes, fh1UnfRhoEffAccCorr,fNIters,false);
    fh1PreFsr= copy(fFSRBayes->Hreco(), "h1preFsr_" + fTag,
		    "h1preFsr_" + fTag
		    + TString(";") + fh1Yield->GetXaxis()->GetTitle()
		    + TString(";pre-FSR full space yield")
		    ,
		    fh1Yield); // set x binning
    TH1D* h1PreFsrCS_raw= copy(fh1PreFsr,"h1PreFsrCS",fTag);
    fh1PreFsrCS= perMassBinWidth(h1PreFsrCS_raw);
    h1Out= fh1PreFsrCS;
  }
  return h1Out;
}

// --------------------------------------------------------------

TH1D* CrossSection_t::calcCrossSection(TVaried_t new_var, int idx)
{

  if (fh1UnfRhoEffAccCorr)
    { delete fh1UnfRhoEffAccCorr; fh1UnfRhoEffAccCorr=NULL; }
  if (fh1PreFsr) { delete fh1PreFsr; fh1PreFsr=NULL; }
  if (fh1PreFsrCS) { delete fh1PreFsrCS; fh1PreFsrCS=NULL; }

  int code=0;
  int res= checkPtrs(&code);
  if (!res && (code<6)) return NULL;

  TString useTag= fTag + variedVarName(new_var) + Form("_%d",idx);

  if (new_var <= _varBkg) {
    if (new_var==_varYield) {
      fh1Signal= copy(fh1Varied,"h1Signal_var",useTag);
      fh1Signal->Add(fh1Bkg,-1);
    }
    else if (new_var==_varBkg) {
      fh1Signal= copy(fh1Yield,"h1Signal_var",useTag);
      fh1Signal->Add(fh1Varied,-1);
    }
    if (fDetResBayes) delete fDetResBayes;
    fDetResBayes= new RooUnfoldBayes( fDetRes, fh1Signal, fNIters, false );
    fh1Unf= copy(fDetResBayes->Hreco(), "h1Unf_" + useTag,
		 "h1Unf_" + useTag
		 + TString(";") + fh1Yield->GetXaxis()->GetTitle()
		 + TString(";unfolded yield")
		 ,
		 fh1Yield); // set x binning
    fh1UnfRhoCorr= copy(fh1Unf,"h1UnfRho",useTag);
    fh1UnfRhoCorr->Divide(fh1Rho);
    if (!fh1EffAcc) {
      fh1UnfRhoEffCorr= copy(fh1UnfRhoCorr,"h1UnfRhoEff",useTag);
      fh1UnfRhoEffCorr->Divide(fh1Eff);
      fh1UnfRhoEffAccCorr= copy(fh1UnfRhoEffCorr,"h1UnfRhoEffAcc",useTag);
      fh1UnfRhoEffAccCorr->Divide(fh1Acc);
    }
    else {
      if (fh1UnfRhoEffCorr) { delete fh1UnfRhoEffCorr; fh1UnfRhoEffCorr=NULL; }
      fh1UnfRhoEffAccCorr= copy(fh1UnfRhoCorr,"h1UnfRhoEffAcc",useTag);
      fh1UnfRhoEffAccCorr->Divide(fh1EffAcc);
    }
    fh1UnfRhoEffAccCorr->Scale(1/fLumi);
    if (fCSType==_csPreFsrFullSp) {
      if (fFSRBayes) delete fFSRBayes;
      fFSRBayes= new RooUnfoldBayes( fFSRRes, fh1UnfRhoEffAccCorr, fNIters, false);
      fh1PreFsr= copy(fFSRBayes->Hreco(), "h1preFsr_" + useTag,
		      "h1preFsr_" + useTag
		      + TString(";") + fh1Yield->GetXaxis()->GetTitle()
		      + TString(";pre-FSR full space yield")
		      ,
		      fh1Yield); // set x binning

      TH1D* h1PreFsrCS_raw= copy(fh1PreFsr,"h1PreFsrCS",useTag);
      fh1PreFsrCS= perMassBinWidth(h1PreFsrCS_raw);
    }
  }
  if (new_var == _varSig) {
    if (fDetResBayes) delete fDetResBayes;
    fDetResBayes= new RooUnfoldBayes( fDetRes, fh1Varied, fNIters, false );
    fh1Unf= copy(fDetResBayes->Hreco(), "h1Unf_" + useTag,
		 "h1Unf_" + useTag
		 + TString(";") + fh1Yield->GetXaxis()->GetTitle()
		 + TString(";unfolded yield")
		 ,
		 fh1Yield); // set x binning
    fh1UnfRhoCorr= copy(fh1Unf,"h1UnfRho",useTag);
    fh1UnfRhoCorr->Divide(fh1Rho);
    if (!fh1EffAcc) {
      fh1UnfRhoEffCorr= copy(fh1UnfRhoCorr,"h1UnfRhoEff",useTag);
      fh1UnfRhoEffCorr->Divide(fh1Eff);
      fh1UnfRhoEffAccCorr= copy(fh1UnfRhoEffCorr,"h1UnfRhoEffAcc",useTag);
      fh1UnfRhoEffAccCorr->Divide(fh1Acc);
    }
    else {
      if (fh1UnfRhoEffCorr) { delete fh1UnfRhoEffCorr; fh1UnfRhoEffCorr=NULL; }
      fh1UnfRhoEffAccCorr= copy(fh1UnfRhoCorr,"h1UnfRhoEffAcc",useTag);
      fh1UnfRhoEffAccCorr->Divide(fh1EffAcc);
    }
    fh1UnfRhoEffAccCorr->Scale(1/fLumi);
    if (fCSType==_csPreFsrFullSp) {
      if (fFSRBayes) delete fFSRBayes;
      fFSRBayes= new RooUnfoldBayes( fFSRRes, fh1UnfRhoEffAccCorr, fNIters, false);
      fh1PreFsr= copy(fFSRBayes->Hreco(), "h1preFsr_" + useTag,
		      "h1preFsr_" + useTag
		      + TString(";") + fh1Yield->GetXaxis()->GetTitle()
		      + TString(";pre-FSR full space yield")
		      ,
		      fh1Yield); // set x binning

      TH1D* h1PreFsrCS_raw= copy(fh1PreFsr,"h1PreFsrCS",useTag);
      fh1PreFsrCS= perMassBinWidth(h1PreFsrCS_raw);
    }
  }
  else if (new_var == _varRho) {
    fh1UnfRhoCorr= copy(fh1Unf,"h1UnfRho",useTag);
    fh1UnfRhoCorr->Divide(fh1Varied);
    if (!fh1EffAcc) {
      fh1UnfRhoEffCorr= copy(fh1UnfRhoCorr,"h1UnfRhoEff",useTag);
      fh1UnfRhoEffCorr->Divide(fh1Eff);
      fh1UnfRhoEffAccCorr= copy(fh1UnfRhoEffCorr,"h1UnfRhoEffAcc",useTag);
      fh1UnfRhoEffAccCorr->Divide(fh1Acc);
    }
    else {
      if (fh1UnfRhoEffCorr) { delete fh1UnfRhoEffCorr; fh1UnfRhoEffCorr=NULL; }
      fh1UnfRhoEffAccCorr= copy(fh1UnfRhoCorr,"h1UnfRhoEffAcc",useTag);
      fh1UnfRhoEffAccCorr->Divide(fh1EffAcc);
    }
    fh1UnfRhoEffAccCorr->Scale(1/fLumi);
    if (fCSType==_csPreFsrFullSp) {
      if (fFSRBayes) delete fFSRBayes;
      fFSRBayes= new RooUnfoldBayes( fFSRRes, fh1UnfRhoEffAccCorr, fNIters, false);
      fh1PreFsr= copy(fFSRBayes->Hreco(), "h1preFsr_" + useTag,
		      "h1preFsr_" + useTag
		      + TString(";") + fh1Yield->GetXaxis()->GetTitle()
		      + TString(";pre-FSR full space yield")
		      ,
		      fh1Yield); // set x binning

      TH1D* h1PreFsrCS_raw= copy(fh1PreFsr,"h1PreFsrCS",useTag);
      fh1PreFsrCS= perMassBinWidth(h1PreFsrCS_raw);
    }
  }
  else if (new_var == _varEff) {
    fh1UnfRhoEffCorr= copy(fh1UnfRhoCorr,"h1UnfRhoEff",useTag);
    fh1UnfRhoEffCorr->Divide(fh1Varied);
    fh1UnfRhoEffAccCorr= copy(fh1UnfRhoEffCorr,"h1UnfRhoEffAcc",useTag);
    fh1UnfRhoEffAccCorr->Divide(fh1Acc);
    fh1UnfRhoEffAccCorr->Scale(1/fLumi);
    if (fCSType==_csPreFsrFullSp) {
      if (fFSRBayes) delete fFSRBayes;
      fFSRBayes= new RooUnfoldBayes( fFSRRes, fh1UnfRhoEffAccCorr, fNIters, false);
      fh1PreFsr= copy(fFSRBayes->Hreco(), "h1preFsr_" + useTag,
		      "h1preFsr_" + useTag
		      + TString(";") + fh1Yield->GetXaxis()->GetTitle()
		      + TString(";pre-FSR full space yield")
		      ,
		      fh1Yield); // set x binning

      TH1D* h1PreFsrCS_raw= copy(fh1PreFsr,"h1PreFsrCS",useTag);
      fh1PreFsrCS= perMassBinWidth(h1PreFsrCS_raw);
    }
  }
  else if (new_var == _varAcc) {
    if (!fh1UnfRhoEffCorr || !fh1Varied) {
      std::cout << "null ptr\n";
    }
    fh1UnfRhoEffAccCorr= copy(fh1UnfRhoEffCorr,"h1UnfRhoEffAcc",useTag);
    fh1UnfRhoEffAccCorr->Divide(fh1Varied);
    fh1UnfRhoEffAccCorr->Scale(1/fLumi);
    if (fCSType==_csPreFsrFullSp) {
      if (fFSRBayes) delete fFSRBayes;
      fFSRBayes= new RooUnfoldBayes( fFSRRes, fh1UnfRhoEffAccCorr, fNIters, false);
      fh1PreFsr= copy(fFSRBayes->Hreco(), "h1preFsr_" + useTag,
		      "h1preFsr_" + useTag
		      + TString(";") + fh1Yield->GetXaxis()->GetTitle()
		      + TString(";pre-FSR full space yield")
		      ,
		      fh1Yield); // set x binning

      TH1D* h1PreFsrCS_raw= copy(fh1PreFsr,"h1PreFsrCS",useTag);
      fh1PreFsrCS= perMassBinWidth(h1PreFsrCS_raw);
    }
  }
  else if (new_var == _varEffAcc) {
    fh1UnfRhoEffAccCorr= copy(fh1UnfRhoCorr,"h1UnfRhoEffAcc",useTag);
    fh1UnfRhoEffAccCorr->Divide(fh1Varied);
    fh1UnfRhoEffAccCorr->Scale(1/fLumi);
    if (fCSType==_csPreFsrFullSp) {
      if (fFSRBayes) delete fFSRBayes;
      fFSRBayes= new RooUnfoldBayes( fFSRRes, fh1UnfRhoEffAccCorr, fNIters, false);
      fh1PreFsr= copy(fFSRBayes->Hreco(), "h1preFsr_" + useTag,
		      "h1preFsr_" + useTag
		      + TString(";") + fh1Yield->GetXaxis()->GetTitle()
		      + TString(";pre-FSR full space yield")
		      ,
		      fh1Yield); // set x binning

      TH1D* h1PreFsrCS_raw= copy(fh1PreFsr,"h1PreFsrCS",useTag);
      fh1PreFsrCS= perMassBinWidth(h1PreFsrCS_raw);
      //fh1PreFsrCS->GetXaxis()->SetTitle("M_{ee} [GeV]");
      fh1PreFsrCS->GetYaxis()->SetTitle("pre-FSR full space #sigma");
    }
  }

  TH1D *h1Out=NULL;
  if (fCSType==_csPostFsrFullSp) {
    h1Out= fh1UnfRhoEffAccCorr;
  }
  else if (fCSType==_csPreFsrFullSp) {
    h1Out= fh1PreFsrCS;
  }
  return h1Out;
}


// --------------------------------------------------------------

TCanvas* CrossSection_t::plotCrossSection(TString canvName)
{
  TH1D *h1=NULL;
  if (fCSType==_csPreFsrFullSp) h1=fh1PreFsrCS;
  else if (fCSType==_csPostFsrFullSp) h1=fh1UnfRhoEffAccCorr;
  if (!h1) h1=calcCrossSection();
  if (!h1) {
    std::cout << "plotCrossSection: pre-calculated cross-section "
	      << "of type " << csTypeName(fCSType)
	      << "is not available and could not be computed\n";
    return NULL;
  }
  TCanvas *c= new TCanvas(canvName,canvName,600,600);
  c->SetLogx();
  c->SetLogy();
  TString drawOpt="LPE";
  if (fh1Theory) {
    fh1Theory->SetLineColor(kBlue);
    fh1Theory->Draw("hist");
    drawOpt.Append("same");
  }
  h1->SetMarkerStyle(24);
  h1->Draw(drawOpt);
  c->Update();
  return c;
}

// --------------------------------------------------------------

int CrossSection_t::sampleRndVec(TVaried_t new_var, int sampleSize,
				 std::vector<TH1D*> &rndCS)
{
  rndCS.clear();
  rndCS.reserve(sampleSize);
  if (!calcCrossSection()) return 0;
  const TH1D *h1Src= getVariedHisto(new_var);
  if (!h1Src) return 0;
  int trackRnd=1;
  TH1D *h1Track=NULL;
  if (trackRnd) {
    h1Track= copy(h1Src,"h1Track",fTag);
    h1Track->SetMarkerStyle(24);
    plotHisto(h1Track, "cVaried",1);
    plotHistoSame(h1Track,"cVaried","LPE");
  }
  if (!fh1Varied) {
    fh1Varied= copy(h1Src,"hVar" + variedVarName(new_var), fTag);
    fh1Varied->Reset();
  }
  for (int i=0; i<sampleSize; i++) {
    std::cout << "\n\nrandomize: i=" << i << ", " << variedVarName(new_var) << "\n";
    const int nonNegative=0;
    randomizeWithinErr(h1Src, fh1Varied, nonNegative);
    if (trackRnd) {
      //fh1Varied->SetMarkerStyle(24);
      int color=(i%9)+1;
      TH1D *h1var= copy(fh1Varied,"h1var" + TString(Form("_%d",i)), fTag);
      h1var->SetLineColor(color);
      h1var->SetMarkerColor(color);
      plotHistoSame(h1var, "cVaried", "hist");
    }
    TH1D* h1=copy(calcCrossSection(new_var,i));
    rndCS.push_back(h1);
  }
  return 1;
}

// --------------------------------------------------------------

int CrossSection_t::sampleRndVec(TVaried_t new_var,
				 const std::vector<TH1D*> &rndHistos,
				 std::vector<TH1D*> &rndCS)
{
  rndCS.clear();
  rndCS.reserve(rndHistos.size());
  if (!calcCrossSection()) return 0;
  const TH1D *h1Src= getVariedHisto(new_var);
  if (!h1Src) return 0;
  int trackRnd=1;
  TH1D *h1Track=NULL;
  if (trackRnd) {
    h1Track= copy(h1Src,"h1Track",fTag);
    h1Track->SetMarkerStyle(24);
    plotHisto(h1Track, "cVaried",1);
    plotHistoSame(h1Track,"cVaried","LPE");
  }
  if (!fh1Varied) {
    fh1Varied= copy(h1Src,"hVar" + variedVarName(new_var), fTag);
    fh1Varied->Reset();
  }
  for (unsigned int i=0; i<rndHistos.size(); i++) {
    std::cout << "\n\nuse: i=" << i << ", " << variedVarName(new_var) << "\n";
    fh1Varied= copy(rndHistos[i],"h1var" + TString(Form("_copy%d",i)), fTag);
    if (trackRnd) {
      //fh1Varied->SetMarkerStyle(24);
      int color=(i%9)+1;
      TH1D *h1var= copy(fh1Varied,"h1var" + TString(Form("_%d",i)), fTag);
      h1var->SetLineColor(color);
      h1var->SetMarkerColor(color);
      plotHistoSame(h1var, "cVaried", "hist");
    }
    TH1D* h1=copy(calcCrossSection(new_var,i));
    rndCS.push_back(h1);
    delete fh1Varied; fh1Varied=NULL;
  }
  return 1;
}
				 

// --------------------------------------------------------------

int CrossSection_t::deriveCov(const std::vector<TH1D*> &rndCS,
			      TH1D **h1avgCS_out, TH2D **h2cov_out)
{
  int res= deriveCovariance(rndCS,
			    variedVarName(fVar) + fTag,
			    "avg " + csTypeName(fCSType),
			    h1avgCS_out, h2cov_out);
  if (!res) std::cout << "error in CrossSection_t::deriveCov\n";
  return res;

  /*
  int nSize= int(rndCS.size());
  if (nSize<2) return 0;

  TH1D *h1a= copy(rndCS[0], "h1avgCS" + variedVarName(fVar), fTag);
  h1a->Reset();

  // accumulate the average value
  for (unsigned int i=0; i<rndCS.size(); i++) {
    h1a->Add(rndCS[i]);
  }
  h1a->Scale(1/double(nSize));

  // accumulate sqr sums
  int dim=h1a->GetNbinsX();
  TMatrixD sum2(dim+1,dim+1);
  sum2.Zero();
  for (unsigned int i=0; i<rndCS.size(); i++) {
    const TH1D *h1= rndCS[i];
    for (int ibin=1; ibin<=h1->GetNbinsX(); ibin++) {
      for (int jbin=1; jbin<=h1->GetNbinsX(); jbin++) { // ! 1D histo
	sum2(ibin,jbin) += 
	  ( h1->GetBinContent(ibin) - h1a->GetBinContent(ibin) ) *
	  ( h1->GetBinContent(jbin) - h1a->GetBinContent(jbin) );
      }
    }
  }

  sum2 *= 1/double(nSize);

  // create the covariance histogram
  const TArrayD *xb= h1a->GetXaxis()->GetXbins();
  const Double_t *x= xb->GetArray();
  TString h2name= "h2cov_" + variedVarName(fVar) + fTag;
  TH2D *h2= new TH2D(h2name,h2name, xb->GetSize()-1,x, xb->GetSize()-1,x);
  h2->SetDirectory(0);
  h2->Sumw2();
  for (int ibin=1; ibin<=h2->GetNbinsX(); ibin++) {
    for (int jbin=1; jbin<=h2->GetNbinsY(); jbin++) { // ! 2D histo
      h2->SetBinContent( ibin,jbin, sum2(ibin,jbin) );
      h2->SetBinError  ( ibin,jbin, 0.);
    }
  }

  // assign uncertainties, taking into account that the uncertainty
  // cannot be negative
  for (int ibin=1; ibin<=h2->GetNbinsX(); ibin++) {
    double unc2= h2->GetBinContent(ibin,ibin);
    if (unc2<0) {
      std::cout << "negative central covariance in ibin=" << ibin << "\n";
      unc2=0;
    }
    h1a->SetBinError(ibin, sqrt(unc2));
  }

  // assign results
  if (h1avgCS_out) *h1avgCS_out= h1a;
  if (h2cov_out) *h2cov_out= h2;
  return 1;
  */
}

// --------------------------------------------------------------


int CrossSection_t::save(TString fname) const
{
  TFile fout(fname,"RECREATE");
  if (!fout.IsOpen()) {
    std::cout << "failed to create the file <" << fname << ">\n";
    return 0;
  }
  if (fh1Yield) fh1Yield->Write("h1Yield");
  if (fh1Bkg) fh1Bkg->Write("h1Bkg");
  if (fDetRes) fDetRes->Write("detRes");
  if (fFSRRes) fFSRRes->Write("FSRRes");
  if (fh1Eff) fh1Eff->Write("h1Eff");
  if (fh1Rho) fh1Rho->Write("h1Rho");
  if (fh1Acc) fh1Acc->Write("h1Acc");
  if (fh1EffAcc) fh1EffAcc->Write("h1EffAcc");
  if (fh1Signal) fh1Signal->Write("h1Signal");
  if (fh1Unf) fh1Unf->Write("h1Unf");
  if (fh1UnfRhoCorr) fh1UnfRhoCorr->Write("h1UnfRho");
  if (fh1UnfRhoEffCorr) fh1UnfRhoEffCorr->Write("h1UnfRhoEff");
  if (fh1UnfRhoEffAccCorr) fh1UnfRhoEffAccCorr->Write("h1UnfRhoEffAcc");
  if (fh1PreFsr) fh1PreFsr->Write("h1PreFSR");
  if (fh1PreFsrCS) fh1PreFsrCS->Write("h1PreFSRCS");
  if (fh1Theory) fh1Theory->Write("h1Theory");
  if (fh1Varied) fh1Varied->Write("h1Varied");
  if (fDetResBayes) fDetResBayes->Write("detResBayes");
  if (fFSRBayes) fFSRBayes->Write("fsrBayes");

  TObjString infoName(fName);
  infoName.Write("name");
  TObjString infoTag(fTag);
  infoTag.Write("tag");
  TObjString infoVar(Form("%d",fVar));
  infoVar.Write("variedVar");
  TObjString infoCSType(Form("%d",fCSType));
  infoCSType.Write("CSType");
  infoCSType.Write(csTypeName(fCSType));
  TObjString infoNIter(Form("%d",fNIters));
  infoNIter.Write("infoNIter");
  TObjString infoLumi(Form("%lf",fLumi));
  infoLumi.Write("lumi");
  TObjString timeTag(DayAndTimeTag(0));
  timeTag.Write("timeTag");
  TObjString explain("productionTime");
  explain.Write(timeTag.String());

  fout.Close();

  return 1;
}

// --------------------------------------------------------------

int CrossSection_t::load(TString fname, TString setTag)
{
  TFile fin(fname);
  if (!fin.IsOpen()) {
    std::cout << "failed to open the file <" << fname << ">\n";
    return 0;
  }
  if (setTag.Length()) fTag=setTag;
  clearAllPtrs();

  TH1D *h1dummy=NULL;
  fh1Yield= loadHisto(fin,"h1Yield","h1Yield" + fTag, 1, h1dummy);
  fh1Bkg  = loadHisto(fin,"h1Bkg", "h1Bkg" + fTag, 1, h1dummy);
  fDetRes = loadRooUnfoldResponse(fin, "detRes", "detRes" + fTag );
  fFSRRes = loadRooUnfoldResponse(fin, "FSRRes", "FSRRes" + fTag );
  fh1Eff  = loadHisto(fin, "h1Eff", "h1Eff" + fTag, 0, h1dummy);
  fh1Rho  = loadHisto(fin, "h1Rho", "h1Rho" + fTag, 1, h1dummy);
  fh1Acc  = loadHisto(fin, "h1Acc", "h1Acc" + fTag, 0, h1dummy);
  fh1EffAcc=loadHisto(fin, "h1EffAcc", "h1EffAcc" + fTag, 0, h1dummy);
  fh1Signal=loadHisto(fin, "h1Signal", "h1Signal" + fTag, 0, h1dummy);
  fh1Unf  = loadHisto(fin, "h1Unf", "h1Unf" + fTag, 0, h1dummy);
  fh1UnfRhoCorr= loadHisto(fin, "h1UnfRho", "h1UnfRho" + fTag, 0, h1dummy);
  fh1UnfRhoEffCorr= loadHisto(fin, "h1UnfRhoEffAcc", "h1UnfRhoEffAcc" + fTag, 0, h1dummy);
  fh1UnfRhoEffAccCorr= loadHisto(fin, "h1UnfRhoEffAcc", "h1UnfRhoEffAcc" + fTag, 0, h1dummy);
  fh1PreFsr = loadHisto(fin, "h1PreFSR", "h1PreFSR" + fTag, 0, h1dummy);
  fh1PreFsrCS=loadHisto(fin, "h1PreFSRCS","h1PreFSRCS" + fTag, 0, h1dummy);
  fh1Theory = loadHisto(fin, "h1Theory", "h1Theory" + fTag, 0, h1dummy);
  fh1Varied = loadHisto(fin, "h1Varied", "h1Varied" + fTag, 0, h1dummy);
  fDetResBayes= loadRooUnfoldBayes(fin, "detResBayes", "detResBayes" + fTag);
  fFSRBayes = loadRooUnfoldBayes(fin, "fsrBayes", "fsrBayes" + fTag);

  HERE("load objString");

  int res=1;
  TObjString *infoName=(TObjString*)fin.Get("name");
  if (!infoName) { res=0; HERE("infoName is null"); }
  else fName= infoName->String();
  //TObjString *infoTag=(TObjString*)fin.Get("tag");
  //if (!infoTag) { res=0; HERE("infoTag is null"); }
  //else fTag= infoTag->String();
  TObjString *infoVar=(TObjString*)fin.Get("variedVar");
  if (!infoVar) { res=0; HERE("infoVar is null"); }
  else fVar=TVaried_t(atoi(infoVar->String().Data()));
  TObjString *infoCSType=(TObjString*)fin.Get("CSType");
  if (!infoCSType) { res=0; HERE("infoCSType is null"); }
  else fCSType=TCSType_t(atoi(infoCSType->String().Data()));
  TObjString *infoNIter=(TObjString*)fin.Get("infoNIter");
  if (!infoNIter) { res=0; HERE("infoNIter is null"); }
  else fNIters= atoi(infoNIter->String().Data());
  TObjString *infoLumi= (TObjString*)fin.Get("lumi");
  if (!infoLumi) { res=0; HERE("infoLumi is null"); }
  else fLumi= atof(infoLumi->String().Data());

  if (res) {
    delete infoName;
    //delete infoTag;
    delete infoVar;
    delete infoCSType;
    delete infoNIter;
    delete infoLumi;
  }

  if (res) HERE("all loaded"); else HERE("some fields were missing");

  fin.Close();
  if (!res) return 0;

  int failCode=0;
  res=checkPtrs(&failCode);
  //if (res) std::cout << "load res=" << res << "\n";
  return (failCode<6) ? 0 : 1;
}

// --------------------------------------------------------------

TH1D* CrossSection_t::copy(const TH1D* h1_orig,
			   TString newName, TString useTag) const
{
  if (!h1_orig) {
    std::cout << "error CrossSection_t::copy(NULL,newName,tag)\n";
    return NULL;
  }
  TH1D* h1=(TH1D*) h1_orig->Clone(newName + useTag);
  h1->SetDirectory(0);
  h1->SetStats(0);
  if (!h1->GetSumw2()) h1->Sumw2();
  return h1;
}

// --------------------------------------------------------------

TH1D* CrossSection_t::copy(const TH1D* h1_orig, TString useTag) const
{
  if (!h1_orig) {
    std::cout << "error CrossSection_t::copy(NULL,tag)\n";
    return NULL;
  }
  TH1D* h1=(TH1D*) h1_orig->Clone(h1_orig->GetName() + useTag);
  h1->SetDirectory(0);
  h1->SetStats(0);
  if (!h1->GetSumw2()) h1->Sumw2();
  return h1;
}

// --------------------------------------------------------------

TH1D* CrossSection_t::copy(const TH1 *h1unf, TString newName, TString newTitle,
			   const TH1D *h1protoType) const
{
  if (h1unf->GetNbinsX() != h1protoType->GetNbinsX()) {
    std::cout << "copy(TH1*): the number of bins is different: "
	      << h1unf->GetNbinsX()<<" vs "<<h1protoType->GetNbinsX() << "\n";
    std::cout << "h1unf->GetName=" << h1unf->GetName() << "\n";
    std::cout << "h1protoType->GetName=" << h1protoType->GetName() << "\n";
    return NULL;
  }
  TH1D *h1=(TH1D*)h1protoType->Clone(newName);
  h1->SetTitle(newTitle);
  h1->SetDirectory(0);
  h1->SetStats(0);
  if (!h1->GetSumw2()) h1->Sumw2();

  for (int ibin=1; ibin<=h1unf->GetNbinsX(); ibin++) {
    h1->SetBinContent(ibin, h1unf->GetBinContent(ibin));
    h1->SetBinError  (ibin, h1unf->GetBinError  (ibin));
  }
  return h1;
}

// --------------------------------------------------------------
// --------------------------------------------------------------

MuonCrossSection_t::MuonCrossSection_t(TString setName, TString setTag,
				       double lumiA, double lumiB) :
  fCSa(setName + TString("a"), setTag + TString("a"), _csPostFsrFullSp),
  fCSb(setName + TString("b"), setTag + TString("b"), _csPostFsrFullSp),
  fTag(setTag),
  fh1BkgV(), fBkgWeights(),
  fh1Bkg(NULL),
  fh1CS(NULL),
  fh1Theory(NULL),
  fFSRBayes(NULL)
{
  fCSa.lumi(lumiA);
  fCSb.lumi(lumiB);
}

// --------------------------------------------------------------

void MuonCrossSection_t::h1Bkg(const TH1D *h1, int clearVec)
{
  if (clearVec) {
    fh1BkgV.clear();
    fBkgWeights.Clear();
  }

  double lumiTot= fCSa.lumi() + fCSb.lumi();
  TH1D *h1bkgA= fCSa.copy(h1);
  h1bkgA->Scale( fCSa.lumi()/lumiTot );
  fCSa.h1Bkg(h1bkgA);
  TH1D *h1bkgB= fCSb.copy(h1);
  h1bkgB->Scale( fCSb.lumi()/lumiTot );
  fCSb.h1Bkg(h1bkgB);

  fh1Bkg=fCSa.copy(h1, "h1BkgTot", fTag);
}

// --------------------------------------------------------------

void MuonCrossSection_t::setBkgV(const std::vector<TH1D*> &setBkg,
				 const std::vector<double> &weights)
{
  if (fh1Bkg) { delete fh1Bkg; fh1Bkg=NULL; }
  fh1BkgV.clear();
  fBkgWeights.Clear();
  if (setBkg.size()==0) {
    std::cout << "setBkgV: size is 0\n";
    return;
  }
  TH1D *h1BkgSum= fCSa.copy(setBkg[0],fTag);
  for (unsigned int i=0; i<setBkg.size(); i++) {
    TH1D *h1=fCSa.copy(setBkg[i],fTag);
    h1->Scale(weights[i]);
    fh1BkgV.push_back(h1);
    if (i>0) h1BkgSum->Add(h1);
  }
  fBkgWeights.ResizeTo(weights.size());
  for (unsigned int i=0; i<weights.size(); i++) {
    fBkgWeights[i]= weights[i];
  }
  h1Bkg(h1BkgSum,0);
}

// --------------------------------------------------------------

int MuonCrossSection_t::recalcBkg(const std::vector<double> &weights)
{
  if (fh1Bkg) { delete fh1Bkg; fh1Bkg=NULL; }
  if (fh1BkgV.size()==0) {
    std::cout << "recalcBkg: fh1BkgV.size=0\n";
    return 0;
  }
  if (fh1BkgV.size()!=weights.size()) {
    std::cout << "recalcBkg: fh1BkgV.size=" << fh1BkgV.size()
	      << ", weights.size=" << weights.size() << "\n";
    return 0;
  }
  TH1D *h1BkgSum= fCSa.copy(fh1BkgV[0],fTag);
  h1BkgSum->Scale(weights[0]);
  for (unsigned int i=1; i<fh1BkgV.size(); i++) {
    TH1D *h1=fCSa.copy(fh1BkgV[i],Form("_tmp%d_%s",i,fTag.Data()));
    h1->Scale(weights[i]);
    h1BkgSum->Add(h1);
    delete h1;
  }
  h1Bkg(h1BkgSum,0);
  return 1;
}

// --------------------------------------------------------------

TH1D* MuonCrossSection_t::calcCrossSection()
{
  if (!fh1Bkg) {
    std::cout << "MuonCrossSection_t::calcCrossSection: object is not ready\n";
    return NULL;
  }
  TH1D *h1postFSR= fCSa.copy(fCSa.calcCrossSection(),"h1postFSR",fTag);
  TH1D *h1CSb= fCSb.copy(fCSb.calcCrossSection());
  if (!h1postFSR || !h1CSb) {
    std::cout << "MuonCrossSection_t::calcCrossSection - failure\n";
    return NULL;
  }
  HERE("scale");
  h1postFSR->Scale( fCSa.lumi() );
  h1CSb->Scale( fCSb.lumi() );
  h1postFSR->Add(h1CSb);

  if (fFSRBayes) delete fFSRBayes;
  fFSRBayes= new RooUnfoldBayes( fCSa.fFSRRes, h1postFSR, fCSa.fNIters, false);
  TH1D* h1CS_raw= (TH1D*) fFSRBayes->Hreco()->Clone("h1preFSR_" + fTag);
  h1CS_raw->SetDirectory(0);
  h1CS_raw->SetTitle("h1preFSR_" + fTag +
		     TString(";") + fh1Bkg->GetXaxis()->GetTitle()
		     + TString("; pre-FSR full space cross section"));
  double lumiTot= fCSa.lumi() + fCSb.lumi();
  h1CS_raw->Scale(1/lumiTot);
  fh1CS= perMassBinWidth(h1CS_raw,0);
  return fh1CS;
}

// --------------------------------------------------------------

TCanvas* MuonCrossSection_t::plotCrossSection(TString canvName,int recalculate)
{
  TH1D *h1= NULL;
  if (recalculate && fh1CS) { delete fh1CS; fh1CS=NULL; }
  if (!fh1CS) calcCrossSection();
  if (fh1CS) h1= cloneHisto(fh1CS, fh1CS->GetName() + canvName, fh1CS->GetTitle());
  if (!h1) {
    std::cout << "MuonCrossSection_t::plotCrossSection: "
	      << "pre-calculated cross-section "
	      << "is not available and could not be computed\n";
    return NULL;
  }
  TCanvas *c= new TCanvas(canvName,canvName,600,600);
  c->SetLogx();
  c->SetLogy();
  TString drawOpt="LPE";
  if (fh1Theory) {
    fh1Theory->SetLineColor(kBlue);
    fh1Theory->Draw("hist");
    drawOpt.Append("same");
  }
  h1->SetMarkerStyle(24);
  h1->Draw(drawOpt);
  c->Update();
  return c;
}

// --------------------------------------------------------------

int MuonCrossSection_t::save(TString fnameBase)
{
  TString fnameA= fnameBase + "_csA.root";
  TString fnameB= fnameBase + "_csB.root";
  int res= (fCSa.save(fnameA) && fCSb.save(fnameB)) ? 1:0;
  if (!res) {
    std::cout << "MuonCrossSection_t::save: failed to save individual "
	      << "CS's\n";
    return 0;
  }

  TString fname= fnameBase + "_cs.root";
  TFile fout(fname,"RECREATE");
  if (!fout.IsOpen()) {
    std::cout << "failed to create the file <" << fname << ">\n";
    return 0;
  }
  fBkgWeights.Write("bkgWeights");
  fout.mkdir("bkg");
  fout.cd("bkg");
  for (unsigned int i=0; i<fh1BkgV.size(); i++) fh1BkgV[i]->Write();
  fout.cd();
  if (fh1Bkg) fh1Bkg->Write("h1Bkg");
  if (fh1CS) fh1CS->Write("h1CS");
  if (fh1Theory) fh1Theory->Write("h1Theory");
  if (fFSRBayes) fFSRBayes->Write("fsrBayes");

  TObjString infoTag(fTag);
  infoTag.Write("tag");
  TObjString timeTag(DayAndTimeTag(0));
  timeTag.Write("timeTag");
  TObjString explain("productionTime");
  explain.Write(timeTag.String());

  fout.Close();
  return 1;
}

// --------------------------------------------------------------

int MuonCrossSection_t::load(TString fnameBase, TString setTag)
{
  this->clear();
  fTag= setTag;

  TString fnameA= fnameBase + "_csA.root";
  TString fnameB= fnameBase + "_csB.root";
  int res= (fCSa.load(fnameA, setTag + "a")
	    && fCSb.load(fnameB, setTag + "b")) ? 1:0;
  if (!res) {
    std::cout << "MuonCrossSection_t::load: failed to load individual "
	      << "CS's\n";
    return 0;
  }

  TString fname= fnameBase + "_cs.root";
  TFile fin(fname);
  if (!fin.IsOpen()) {
    std::cout << "failed to open the file <" << fname << ">\n";
    return 0;
  }
  fBkgWeights.Read("bkgWeights");
  if (fBkgWeights.GetNoElements()>0) {
    fh1BkgV.reserve(fBkgWeights.GetNoElements());
    fin.cd("bkg");
    TList* objList= gDirectory->GetListOfKeys();
    if (!objList) {
      std::cout << "objList is null\n";
      return 0;
    }
    TIter objs(objList);
    TKey* key=NULL;
    while ( (key= (TKey*)objs()) ) {
      TObject *obj= key->ReadObj();
      if (obj->IsA()->InheritsFrom(TH1D::Class())) {
	TH1D* h1=(TH1D*)obj;
	fh1BkgV.push_back(cloneHisto(h1, h1->GetName()+fTag, h1->GetTitle()));
	//HERE("delete h1");
	delete h1;
      }
      else {
	//HERE("delete obj");
	//delete obj;
      }
    }
    //HERE("delete objList");
    //delete objList;
  }
  fin.cd();
  TH1D* h1dummy=NULL;
  fh1Bkg= loadHisto(fin, "h1Bkg", "h1Bkg" + fTag, 1, h1dummy);
  fh1CS = loadHisto(fin, "h1CS", "h1CS" + fTag, 1, h1dummy);
  fh1Theory= loadHisto(fin, "h1Theory", "h1Theory"+fTag, 1, h1dummy);
  fFSRBayes= loadRooUnfoldBayes(fin, "fsrBayes","fsrBayes" + fTag);

  //TObjString *infoTag=(TObjString*)fin.Get("tag");
  //if (!infoTag) { res=0; HERE("infoTag is null"); }
  //else fTag= infoTag->String();

  fin.Close();
  return 1;
}

// --------------------------------------------------------------
// --------------------------------------------------------------
// --------------------------------------------------------------