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
  case _varBkgXS: name="varBkgXS"; break;
  case _varSig: name="varSig"; break;
  case _varDetRes: name="varDetRes"; break;
  case _varFSRRes: name="varFSRRes"; break;
  case _varFSRRes_Poisson: name="varFSRResPoisson"; break;
  case _varEff: name="varEff"; break;
  case _varRho: name="varRho"; break;
  case _varRhoFile: name="varRhoFile"; break;
  case _varAcc: name="varAcc"; break;
  case _varEffAcc: name="varEffAcc"; break;
  case _varLast: name="varLast"; break;
  case _varRhoSyst: name="varRhoSyst"; break;
  case _varTheory: name="varTheory"; break;
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

TString tnpSystName(TTnPSystType_t s, int shortName)
{
  TString name="tnpSystUNKNOWN";
  if (shortName) {
    switch(s) {
    case _tnpSyst_none: name="tnpSystNone"; break;
    case _tnpSyst_bkgPdf: name="bkgPdf"; break;
    case _tnpSyst_sigPdf: name="sigPdf"; break;
    case _tnpSyst_NLOvsLO: name="NLOvsLO"; break;
    case _tnpSyst_tagDef: name="tag"; break;
    case _tnpSyst_last: name="tnpSystLast"; break;
    default:
      std::cout << "tnpSystName (shortName=1) is not ready for this "
		<< "TTnpSystType value\n";
    }
  }
  else {
    switch(s) {
    case _tnpSyst_none: name="tnpSystNone"; break;
    case _tnpSyst_bkgPdf: name="tnpSystBkgPdf"; break;
    case _tnpSyst_sigPdf: name="tnpSystSigPdf"; break;
    case _tnpSyst_NLOvsLO: name="tnpSystNLOvsLO"; break;
    case _tnpSyst_tagDef: name="tnpSystTagDef"; break;
    case _tnpSyst_last: name="tnpSystLast"; break;
    default:
      std::cout << "tnpSystName (shortName=0) is not ready for this "
		<< "TTnpSystType value\n";
    }
  }
  return name;
}

// --------------------------------------------------------------
// --------------------------------------------------------------

RooUnfoldResponse *loadRooUnfoldResponse(TString fname, TString fieldName,
					 TString setName, int warnIfAbsent)
{
  TFile fin(fname);
  RooUnfoldResponse *rs= loadRooUnfoldResponse(fin,fieldName,setName,warnIfAbsent);
  fin.Close();
  return rs;
}

// --------------------------------------------------------------

RooUnfoldResponse *loadRooUnfoldResponse(TFile &fin, TString fieldName,
					 TString setName, int warnIfAbsent)
{
  if (!fin.IsOpen()) {
    std::cout << "loadRooUnfoldResponse(TFile&): file is not open\n";
    return NULL;
  }
  TObject *chk= fin.Get(fieldName);
  if (!chk) {
    if (warnIfAbsent) {
      std::cout << "loadRooUnfoldResponse(TFile&): field with name=<"
		<< fieldName << "> does not exist in file <"
		<< fin.GetName() << ">\n";
    }
    return NULL;
  }
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

RooUnfoldResponse* randomizeWithinErr(const RooUnfoldResponse *R, TString tag,
				      int poissonRnd)
{
  TH2D *h2MigRnd= cloneHisto(R->Hresponse(), "h2MigRnd_"+tag, R->Hresponse()->GetTitle(),h2dummy);
  if (!h2MigRnd) {
    std::cout << "error in randomizeWithinErr(RooUnfoldResponse)\n";
    return NULL;
  }
  int nonNegative=0;
  randomizeWithinErr((TH2D*)R->Hresponse(), h2MigRnd, nonNegative, poissonRnd);

  TH1D* h1measRnd= cloneHisto(R->Hmeasured(),"h1measRnd"+tag,"h1measRnd",h1dummy);
  TH1D* h1trueRnd= cloneHisto(R->Htruth(),"h1trueRnd"+tag,"h1trueRnd",h1dummy);


  h1measRnd->Reset();
  h1trueRnd->Reset();

  for (int ibin=1; ibin<= h2MigRnd->GetNbinsX(); ibin++) {
    double sum=0;
    double sumErr2=0;
    for (int jbin=1; jbin<= h2MigRnd->GetNbinsY(); jbin++) {
      sum+= h2MigRnd->GetBinContent(ibin,jbin);
      sumErr2+= pow(h2MigRnd->GetBinError(ibin,jbin),2);
    }
    h1measRnd->SetBinContent(ibin, sum);
    h1measRnd->SetBinError  (ibin, sqrt(sumErr2));
  }

  for (int jbin=1; jbin<= h2MigRnd->GetNbinsY(); jbin++) {
    double sum=0;
    double sumErr2=0;
    for (int ibin=1; ibin<= h2MigRnd->GetNbinsX(); ibin++) {
      sum+= h2MigRnd->GetBinContent(ibin,jbin);
      sumErr2+= pow(h2MigRnd->GetBinError(ibin,jbin),2);
    }
    h1trueRnd->SetBinContent(jbin, sum);
    h1trueRnd->SetBinError  (jbin, sqrt(sumErr2));
  }

  RooUnfoldResponse *Rrnd= new RooUnfoldResponse(h1measRnd,h1trueRnd,h2MigRnd,
						 "RespRnd"+tag,
						 R->Hresponse()->GetTitle());
  delete h1measRnd;
  delete h1trueRnd;
  delete h2MigRnd;
  return Rrnd;
}

// --------------------------------------------------------------
// --------------------------------------------------------------

int RndVecInfo_t::loadHistos(int sampleSize,
			     std::vector<TH1D*> &rndVec,
			     std::vector<TH1D*> *rndVec2) const
{
  rndVec.clear();
  if (rndVec2) rndVec2->clear();

  TFile fin(fFName);
  if (!fin.IsOpen()) {
    std::cout << "failed to open the file <" << fin.GetName() << ">\n";
    return 0;
  }

  for (int iRnd=0; iRnd<sampleSize; iRnd++) {
    if ((fMaxSamples>0) && (iRnd>=fMaxSamples)) {
      std::cout << "loaded the requested " << fMaxSamples << " samples\n";
      break;
    }
    TH1D *h1r= loadHisto(fin,histoName(iRnd),histoName(iRnd),
			 1,h1dummy);
    rndVec.push_back(h1r);
    if (!h1r) return 0;
    if (rndVec2) {
      TH1D *h1rV2= loadHisto(fin,histoNameV2(iRnd),histoNameV2(iRnd),
			     1,h1dummy);
      rndVec2->push_back(h1rV2);
      if (!h1rV2) return 0;
    }
  }
  fin.Close();

  return 1;
}

// --------------------------------------------------------------
// --------------------------------------------------------------

CrossSection_t::CrossSection_t(TString setName, TString setTag,
			       TCSType_t set_csType, TVersion_t setVersion) :
  fName(setName), fTag(setTag),
  fVersion(setVersion),
  fh1Yield(NULL), fh1Bkg(NULL),
  fDetRes(NULL), fFSRRes(NULL),
  fh1Eff(NULL), fh1Rho(NULL), fh1Acc(NULL),
  fh1EffAcc(NULL),
  fh1Theory(NULL),
  fVar(_varNone), fCSType(set_csType),
  fNItersDetRes(4), fNItersFSR(4), fLumi(1.),
  fh1Signal(NULL), fh1Unf(NULL), fh1UnfRhoCorr(NULL),
  fh1UnfRhoEffCorr(NULL), fh1UnfRhoEffAccCorr(NULL),
  fh1PreFsr(NULL), fh1PreFsrCS(NULL),
  fh1Varied(NULL), fResVaried(NULL),
  fDetResBayes(NULL), fFSRBayes(NULL)
{
  this->setNIters_internal(fVersion);
}

// --------------------------------------------------------------

CrossSection_t::CrossSection_t(const CrossSection_t &cs,
			       TString setName, TString setTag,
			       TVersion_t setVersion) :
  fName(setName), fTag(setTag),
  fVersion(cs.fVersion),
  fh1Yield(NULL), fh1Bkg(NULL),
  fDetRes(NULL), fFSRRes(NULL),
  fh1Eff(NULL), fh1Rho(NULL), fh1Acc(NULL),
  fh1EffAcc(NULL),
  fh1Theory(NULL),
  fVar(_varNone), fCSType(cs.fCSType),
  fNItersDetRes(4), fNItersFSR(4), fLumi(1.),
  fh1Signal(NULL), fh1Unf(NULL), fh1UnfRhoCorr(NULL),
  fh1UnfRhoEffCorr(NULL), fh1UnfRhoEffAccCorr(NULL),
  fh1PreFsr(NULL), fh1PreFsrCS(NULL),
  fh1Varied(NULL), fResVaried(NULL),
  fDetResBayes(NULL), fFSRBayes(NULL)
{
  if (!this->assign(cs)) {
    std::cout << "error in CrossSection_t::CrossSection_t(CrossSection_t)\n";
  }
  if (setVersion!=_verUndef) fVersion=setVersion;
  this->setNIters_internal(fVersion);
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
  if (fResVaried) { delete fResVaried; fResVaried=NULL; }
  if (fDetResBayes) { delete fDetResBayes; fDetResBayes=NULL; }
  if (fFSRBayes) { delete fFSRBayes; fFSRBayes=NULL; }
}

// --------------------------------------------------------------

void CrossSection_t::setNIters_internal(TVersion_t)
{
  switch(fVersion) {
  case _verMu76X: fNItersDetRes=17; fNItersFSR=100; break;
  case _verMuApproved: fNItersDetRes=17; fNItersFSR=100; break;
  case _verEl2skim: fNItersDetRes=15; fNItersFSR=15; break;
  default: ; // nothing
  }
}


// --------------------------------------------------------------

TH1D* CrossSection_t::getVariedHisto(TVaried_t new_var)
{
  fVar= new_var;
  TH1D* h1=NULL;
  switch(fVar) {
  case _varYield: h1=fh1Yield; break;
  case _varBkg:
  case _varBkgXS: h1=fh1Bkg; break;
  case _varSig: h1=fh1Signal; break;
  case _varDetRes:
  case _varFSRRes:
  case _varFSRRes_Poisson:
    std::cout << "getVariedHisto should not be called for "
	      << variedVarName(new_var) << "\n";
    break;
  case _varEff: h1=fh1Eff; break;
  case _varRho: h1=fh1Rho; break;
  case _varRhoFile: h1=fh1Rho; break;
  case _varAcc: h1=fh1Acc; break;
  case _varEffAcc: h1=fh1EffAcc; break;
  case _varRhoSyst: // no break
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
  fNItersDetRes= cs.fNItersDetRes;
  fNItersFSR= cs.fNItersFSR;
  fLumi= cs.fLumi;
  fh1Signal= copy(cs.fh1Signal);
  fh1Unf= copy(cs.fh1Unf);
  fh1UnfRhoCorr= copy(cs.fh1UnfRhoCorr);
  fh1UnfRhoEffCorr= copy(cs.fh1UnfRhoEffCorr);
  fh1UnfRhoEffAccCorr= copy(cs.fh1UnfRhoEffAccCorr);
  fh1Varied= copy(cs.fh1Varied);
  if (cs.fResVaried) fResVaried= new RooUnfoldResponse(*cs.fResVaried);

  int res= checkPtrs();
  std::cout << "assign res=" << res << "\n";
  return res;
}

// --------------------------------------------------------------

void CrossSection_t::removeBkg()
{
  if (fh1Bkg) {
    for (int ibin=1; ibin<=fh1Bkg->GetNbinsX(); ibin++) {
      fh1Bkg->SetBinContent(ibin, 0.);
      fh1Bkg->SetBinError  (ibin, 0.);
    }
  }
  else std::cout << "removeBkg: fh1Bkg is null\n";
}

// --------------------------------------------------------------

void CrossSection_t::removeRho()
{
  if (fh1Rho) {
    for (int ibin=1; ibin<=fh1Rho->GetNbinsX(); ibin++) {
      fh1Rho->SetBinContent(ibin, 1.);
      fh1Rho->SetBinError  (ibin, 0.);
    }
  }
  else std::cout << "removeRho: fh1Rho is null\n";
}

// --------------------------------------------------------------

TH1D* CrossSection_t::calcCrossSection(int removeNegativeSignal)
{
  if (fh1PreFsr) { delete fh1PreFsr; fh1PreFsr=NULL; }
  if (fh1PreFsrCS) { delete fh1PreFsrCS; fh1PreFsrCS=NULL; }
  int code=0;
  int res= checkPtrs(&code);
  if (!res && (code<6)) return NULL;

  fh1Signal= copy(fh1Yield,"h1Signal",fTag);
  fh1Signal->Add(fh1Bkg,-1);
  if (removeNegativeSignal) removeNegatives(fh1Signal);
  if (fDetResBayes) delete fDetResBayes;
  fDetResBayes= new RooUnfoldBayes( fDetRes, fh1Signal, fNItersDetRes, false );
  fh1Unf= copy(fDetResBayes->Hreco(), "h1Unf_" + fTag,
	       "h1Unf_" + fTag
	       + TString(";") + fh1Yield->GetXaxis()->GetTitle()
	       + TString(";unfolded yield")
       ,
	       fh1Yield); // set x binning
  //HERE("divide by fh1Rho");
  fh1UnfRhoCorr= copy(fh1Unf,"h1UnfRho",fTag);
  //printHisto(fh1Unf);
  //printHisto(fh1Rho);
  //printHisto(fh1EffAcc);
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
  //HERE("divide by lumi");
  TH1D* h1Out=NULL;
  if (fCSType==_csPostFsrFullSp) {
    fh1UnfRhoEffAccCorr->Scale(1/fLumi);
    h1Out= fh1UnfRhoEffAccCorr;
  }
  else if (fCSType==_csPreFsrFullSp) {
    if (fFSRBayes) delete fFSRBayes;

    TH1D *h1UnfFSR_loc= cloneHisto(fh1Yield, "h1UnfFSR_loc","h1UnfFSR_loc");
    h1UnfFSR_loc->Reset();
    if ((fVersion!=_verMu76X) && (fVersion!=_verMuApproved)) {
      std::cout << "fsr unfold Bayes\n";
      fFSRBayes= new RooUnfoldBayes( fFSRRes, fh1UnfRhoEffAccCorr,fNItersFSR,false);
      copyContents(h1UnfFSR_loc,fFSRBayes->Hreco());
    }
    else {
      std::cout << "fsr unfold Invert\n";
      RooUnfoldInvert *rooInv= new RooUnfoldInvert( fFSRRes, fh1UnfRhoEffAccCorr, "rooInv" );
      copyContents(h1UnfFSR_loc, rooInv->Hreco());
    }

    //printRatio(fh1UnfRhoEffAccCorr,h1UnfFSR_loc);

    fh1PreFsr= copy(h1UnfFSR_loc, "h1preFsr_" + fTag,
		    "h1preFsr_" + fTag
		    + TString(";") + fh1Yield->GetXaxis()->GetTitle()
		    + TString(";pre-FSR full space yield")
		    ,
		    fh1Yield); // set x binning
    TH1D* h1PreFsrCS_raw= copy(fh1PreFsr,"h1PreFsrCS",fTag);
    fh1PreFsrCS= perMassBinWidth(h1PreFsrCS_raw);
    fh1PreFsrCS->Scale(1/fLumi);
    h1Out= fh1PreFsrCS;
  }
  return h1Out;
}

// --------------------------------------------------------------

TH1D* CrossSection_t::calcCrossSection(TVaried_t new_var, int idx,
				       int removeNegativeSignal)
{

  int fsrResponse=((new_var==_varFSRRes) || (new_var==_varFSRRes_Poisson))? 1:0;
  int rndResponse=((new_var==_varDetRes) || fsrResponse) ? 1:0;
  if (!fh1Varied && !rndResponse) {
    HERE("entered CrossSection_t::calcCrossSection(TVaried_t). fh1Varied is null");
    return NULL;
  }

  if (fh1UnfRhoEffAccCorr && !fsrResponse)
    { delete fh1UnfRhoEffAccCorr; fh1UnfRhoEffAccCorr=NULL; }
  if (fh1PreFsr) { delete fh1PreFsr; fh1PreFsr=NULL; }
  if (fh1PreFsrCS) { delete fh1PreFsrCS; fh1PreFsrCS=NULL; }

  int code=0;
  int res= checkPtrs(&code);
  if (!res && (code<6)) return NULL;

  TString useTag= fTag + variedVarName(new_var) + Form("_%d",idx);

  if (new_var <= _varBkg) {
    if (new_var==_varYield) {
      if (0) {
	TString cname="cTest" + useTag;
	plotHisto(fh1Yield,cname,1,1,"hist");
	plotHistoSame(fh1Varied,cname,"LPE");
	printHisto(fh1Varied);
      }
      if (1) { // prevent compiler from malfunctioning optimization
	TH1D *h1tmp=fh1Varied;
	HERE(h1tmp->GetName());
	fh1Varied=copy(fh1Yield,"h1Var_xx",useTag);
	for (int ibin=1; ibin<=fh1Varied->GetNbinsX(); ibin++) {
	  fh1Varied->SetBinContent(ibin, h1tmp->GetBinContent(ibin));
	  fh1Varied->SetBinError(ibin,h1tmp->GetBinError(ibin));
	}
	//delete h1tmp; // do not delete - causes trouble
      }
      fh1Signal= copy(fh1Varied,"h1Signal_var",useTag);
      fh1Signal->Add(fh1Bkg,-1);
    }
    else if (new_var==_varBkg) {
      fh1Signal= copy(fh1Yield,"h1Signal_var",useTag);
      fh1Signal->Add(fh1Varied,-1);
    }
    if (removeNegativeSignal) removeNegatives(fh1Signal);
    if (fDetResBayes) delete fDetResBayes;
    fDetResBayes= new RooUnfoldBayes( fDetRes, fh1Signal, fNItersDetRes, false );
    fh1Unf= copy(fDetResBayes->Hreco(), "h1Unf_" + useTag,
		 "h1Unf_" + useTag
		 + TString(";") + fh1Yield->GetXaxis()->GetTitle()
		 + TString(";unfolded yield")
		 ,
		 fh1Yield); // set x binning
    //printHisto(fh1Unf);
    fh1UnfRhoCorr= copy(fh1Unf,"h1UnfRho",useTag);
    //printHisto(fh1UnfRhoCorr);
    fh1UnfRhoCorr->Divide(fh1Rho);
    //printHisto(fh1UnfRhoCorr);
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
      fFSRBayes= new RooUnfoldBayes( fFSRRes, fh1UnfRhoEffAccCorr, fNItersFSR, false);
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
  if ((new_var == _varSig) || (new_var == _varDetRes)) {
    if (fDetResBayes) delete fDetResBayes;
    if (new_var == _varSig) {
      fDetResBayes= new RooUnfoldBayes( fDetRes, fh1Varied, fNItersDetRes, false );
    }
    else {
      fDetResBayes= new RooUnfoldBayes( fResVaried, fh1Signal, fNItersDetRes, false );
    }
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
      fFSRBayes= new RooUnfoldBayes( fFSRRes, fh1UnfRhoEffAccCorr, fNItersFSR, false);
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
  else if ((new_var == _varRho) || (new_var == _varRhoFile)) {
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
      fFSRBayes= new RooUnfoldBayes( fFSRRes, fh1UnfRhoEffAccCorr, fNItersFSR, false);
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
      fFSRBayes= new RooUnfoldBayes( fFSRRes, fh1UnfRhoEffAccCorr, fNItersFSR, false);
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
      std::cout << "varAcc: either fh1UnfRhoEffCorr or fh1Varied are null ptr\n";
      return NULL;
    }
    fh1UnfRhoEffAccCorr= copy(fh1UnfRhoEffCorr,"h1UnfRhoEffAcc",useTag);
    fh1UnfRhoEffAccCorr->Divide(fh1Varied);
    fh1UnfRhoEffAccCorr->Scale(1/fLumi);
    if (fCSType==_csPreFsrFullSp) {
      if (fFSRBayes) delete fFSRBayes;
      fFSRBayes= new RooUnfoldBayes( fFSRRes, fh1UnfRhoEffAccCorr, fNItersFSR, false);
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
      fFSRBayes= new RooUnfoldBayes( fFSRRes, fh1UnfRhoEffAccCorr, fNItersFSR, false);
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
  else if ((new_var == _varFSRRes) || (new_var == _varFSRRes_Poisson)) {
    if (fCSType!=_csPreFsrFullSp) {
      std::cout << "varied varFSRRes, but the cross section is not csPreFsrFullSp\n";
      return NULL;
    }
    if (fFSRBayes) delete fFSRBayes;
    HERE("calling FSR unfold. fResVaried is null %d",int(!fResVaried));
    HERE("fh1UnfRhoEffAccCorr is null %d",int(!fh1UnfRhoEffAccCorr));
    HERE("fNItersFSR=%d",fNItersFSR);
    fFSRBayes= new RooUnfoldBayes( fResVaried, fh1UnfRhoEffAccCorr, fNItersFSR, false);
    HERE("copying unfolded null ptr %d",int(!fFSRBayes->Hreco()));
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

TCanvas* CrossSection_t::plotCrossSection(TString canvName,
					  int removeNegativeSignal)
{
  TH1D *h1=NULL;
  if (fCSType==_csPreFsrFullSp) h1=fh1PreFsrCS;
  else if (fCSType==_csPostFsrFullSp) h1=fh1UnfRhoEffAccCorr;
  if (!h1) h1=calcCrossSection(removeNegativeSignal);
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
  std::cout << "plotCrossSection: theory is blue, calculation is circles\n";
  return c;
}
// --------------------------------------------------------------

TH1D* CrossSection_t::calcCrossSection(TVaried_t new_var, int idx,
				       const TH1D *set_h1varied,
				       int removeNegativeSignal)
{
  int fsrResponse=((new_var==_varFSRRes) || (new_var==_varFSRRes_Poisson))? 1:0;
  int detResponse=((new_var==_varDetRes)) ? 1:0;
  if (fsrResponse || detResponse) {
    std::cout << "calcCrossSection(new_var,idx,TH1D*) should not be called for "
	      << "unfolding task\n";
    return NULL;
  }
  if (!fh1Varied || 1) {
    TH1D *h1=getVariedHisto(new_var);
    fh1Varied=cloneHisto(h1,h1->GetName() + TString("_var%d",idx),
			 h1->GetTitle());
  }
  if (!copyContents(fh1Varied,set_h1varied)) {
    std::cout << "copying the histo failed\n";
    return NULL;
  }
  return calcCrossSection(new_var,idx,removeNegativeSignal);
}

// --------------------------------------------------------------

TCanvas* CrossSection_t::plotCrossSection_StepByStep(TString canvNameBase)
{
  TH1D *h1=NULL;
  if (fCSType==_csPreFsrFullSp) h1=fh1PreFsrCS;
  else if (fCSType==_csPostFsrFullSp) h1=fh1UnfRhoEffAccCorr;
  if (!h1) h1=calcCrossSection();
  if (!h1) {
    std::cout << "plotCrossSection_StepByStep: pre-calculated cross-section "
	      << "of type " << csTypeName(fCSType)
	      << "is not available and could not be computed\n";
    return NULL;
  }

  histoStyle(fh1Yield,46,5,2);
  histoStyle(fh1Bkg,kRed,26,2);
  histoStyle(fh1Signal,kBlack,27,2);
  histoStyle(fh1Unf,kBlue,24,2);
  TString cn= canvNameBase + "_signalUnf";
  plotHisto(fh1Yield,cn,1,1,"LPE","yield");
  plotHistoSame(fh1Bkg,cn,"LPE","bkg");
  plotHistoSame(fh1Signal,cn,"LPE","signal");
  plotHistoSame(fh1Unf,cn,"LPE","unf.signal");

  histoStyle(fh1UnfRhoCorr,kGreen+1,32,3);
  histoStyle(fh1UnfRhoEffAccCorr,kRed,26,2);
  cn= canvNameBase + "_preFSR";
  plotHisto(fh1Unf,cn,1,1,"LPE","unf.signal");
  plotHistoSame(fh1UnfRhoCorr,cn,"LPE","unf. x #rho");
  plotHistoSame(fh1UnfRhoEffAccCorr,cn,"LPE","unf. x #rho x #epsilon x A");
  plotHistoSame(fh1PreFsr,cn,"LPE","pre-FSR");

  histoStyle(fh1PreFsrCS,kBlue,24,1);
  cn= canvNameBase + "_final";
  plotHisto(fh1PreFsr,cn,1,1,"LPE","pre-FSR");
  TCanvas *c=plotHistoSame(fh1PreFsrCS,cn,"LPE","pre-FSR cs");
  return c;
}

// --------------------------------------------------------------

int CrossSection_t::sampleRndVec(TVaried_t new_var, int sampleSize,
				 std::vector<TH1D*> &rndCS)
{
  if ((new_var==_varDetRes) || (new_var==_varFSRRes)
      || (new_var==_varFSRRes_Poisson)) {
    if (!sampleRndResponse(new_var,sampleSize,rndCS,NULL)) {
      std::cout << " called from sampleRndVec(var,sampleSize,rndCS)\n";
      return 0;
    }
    return 1;
  }
  if (new_var==_varAcc) this->removeEffAcc();

  rndCS.clear();
  rndCS.reserve(sampleSize);
  if (!calcCrossSection()) return 0;
  const TH1D *h1Src= getVariedHisto(new_var);
  if (!h1Src) return 0;
  //printHisto(h1Src);
  int trackRnd=1;
  TH1D *h1Track=NULL;
  TString canvVaried="cVaried_" + variedVarName(new_var) + "_" + fTag;
  if (trackRnd) {
    h1Track= copy(h1Src,"h1Track",fTag);
    h1Track->SetMarkerStyle(24);
    plotHisto(h1Track, canvVaried,1);
    plotHistoSame(h1Track,canvVaried,"LPE");
  }
  if (fh1Varied) { delete fh1Varied; fh1Varied=NULL; }
  fh1Varied= copy(h1Src,"hVar" + variedVarName(new_var), fTag);
  fh1Varied->Reset();

  for (int i=0; i<sampleSize; i++) {
    std::cout << "\n\nrandomize: i=" << i << ", " << variedVarName(new_var) << "\n";
    const int nonNegative=0;
    if (0) {
      if (fh1Varied) { delete fh1Varied; fh1Varied=NULL; }
      fh1Varied= copy(h1Src,"hVar" + variedVarName(new_var) + Form("_%d",i), fTag);
      fh1Varied->Reset();
    }
    randomizeWithinErr(h1Src, fh1Varied, nonNegative);
    if (trackRnd) {
      //fh1Varied->SetMarkerStyle(24);
      int color=(i%9)+1;
      TH1D *h1var= copy(fh1Varied,"h1varTmp" + TString(Form("_%d",i)), fTag);
      h1var->SetLineColor(color);
      h1var->SetMarkerColor(color);
      h1var->SetMarkerStyle(5);
      plotHistoSame(h1var, canvVaried, "LPE");
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
  //printHisto(h1Src);
  int trackRnd=1;
  TH1D *h1Track=NULL;
  if (trackRnd) {
    h1Track= copy(h1Src,"h1Track",fTag);
    h1Track->SetMarkerStyle(24);
    plotHisto(h1Track, "cVaried" + fTag,1);
    plotHistoSame(h1Track,"cVaried" + fTag,"LPE");
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
      //printHisto(h1var);
      h1var->SetLineColor(color);
      h1var->SetMarkerColor(color);
      plotHistoSame(h1var, "cVaried" + fTag, "hist");
    }
    TH1D* h1=copy(calcCrossSection(new_var,i));
    rndCS.push_back(h1);
    delete fh1Varied; fh1Varied=NULL;
  }
  return 1;
}

// --------------------------------------------------------------

int CrossSection_t::sampleRndVec(TVaried_t new_var, int sampleSize,
				 const RndVecInfo_t &info,
				 std::vector<TH1D*> &rndCS)
{
  if ((new_var==_varDetRes) || (new_var==_varFSRRes)
      || (new_var==_varFSRRes_Poisson)) {
    std::cout << "CrossSection_t::sampleRndVec(RndVecInfo) should not be "
	      << "called for " << variedVarName(new_var) << "\n";
    return 0;
  }

  rndCS.clear();

  std::vector<TH1D*> rndVec;
  if (!info.loadHistos(sampleSize,rndVec,NULL) ||
      !sampleRndVec(new_var,rndVec,rndCS)) {
    std::cout << "error from CrossSection_t::sampleRndVec(RndVecInfo)\n";
    return 0;
  }
  clearVec(rndVec);
  return 1;
}
				 
// --------------------------------------------------------------

int CrossSection_t::sampleRndResponse(TVaried_t new_var, int sampleSize,
				      std::vector<TH1D*> &rndCS,
				      std::vector<RooUnfoldResponse*> *rndV_out)
{
  rndCS.clear();
  rndCS.reserve(sampleSize);
  TH1D *hCS=calcCrossSection();
  if (!hCS) return 0;
  int code=0;
  std::cout << "sampleRndResponse(" << variedVarName(new_var) << "): "
	    << "checkPtrs=" << checkPtrs(&code) << ", code=" << code << "\n";

  int notRandom= (((new_var==_varFSRRes) || (new_var==_varFSRRes_Poisson)) &&
		  (fCSType!=_csPreFsrInAcc) &&
		  (fCSType!=_csPreFsrFullSp)) ? 1:0;
  if (notRandom) {
      std::cout << "CrossSection_t::sampleRndResponse: the cross section "
		<< "is not pre-FSR, the result will not be random.\n";
      return 0;
  }

  if (new_var!=_varDetRes) {
    if (!fh1UnfRhoEffAccCorr) {
      std::cout << "the needed histogram, h1UnfRhoEffAccCorr, is null\n";
      return 0;
    }
    fh1UnfRhoEffAccCorr->Scale(1/fLumi);
  }

  const RooUnfoldResponse *r=(new_var==_varDetRes) ? fDetRes : fFSRRes;
  if (!r) {
    std::cout << "sampleRndResponse(" << variedVarName(new_var)
	      << ") the object is not ready\n";
    return 0;
  }
  if (fResVaried) { delete fResVaried; fResVaried=NULL; }

  int poissonRnd=0;
  if (new_var==_varFSRRes_Poisson) poissonRnd=1;
  for (int i=0; i<sampleSize; i++) {
    std::cout << "\n\nrandomize: i=" << i << ", " << variedVarName(new_var) << "\n";
    TString respName= variedVarName(new_var) + Form("_%d",i) + fTag;
    fResVaried= randomizeWithinErr(r, respName, poissonRnd);
    fResVaried->UseOverflow(false); // ensure correct behavior
    TH1D* h1=copy(calcCrossSection(new_var,i));
    rndCS.push_back(h1);
    if (rndV_out) {
      rndV_out->push_back(fResVaried);
    }
    else delete fResVaried;
    fResVaried=NULL;
  }

  if (new_var!=_varDetRes) {
    fh1UnfRhoEffAccCorr->Scale(fLumi);
  }


  return 1;
}

// --------------------------------------------------------------

int CrossSection_t::sampleRndRespVec(TVaried_t new_var,
		     const std::vector<RooUnfoldResponse*> &rndRespV,
				     std::vector<TH1D*> &rndCS)
{
  rndCS.clear();
  rndCS.reserve(rndRespV.size());
  if (!calcCrossSection()) return 0;

  int poissonRnd=0;
  if (new_var==_varFSRRes_Poisson) poissonRnd=1;

  for (unsigned int i=0; i<rndRespV.size(); i++) {
    const RooUnfoldResponse *r= rndRespV[i];
    std::cout << "\n\nsampeRndRespV: i=" << i << ", " << variedVarName(new_var) << "\n";

    fResVaried= new RooUnfoldResponse(*r);
    //TString respName= r->GetName() + variedVarName(new_var) + Form("_%d",i) + fTag;
    //fResVaried= randomizeWithinErr(r, respName, poissonRnd);
    TH1D* h1=copy(calcCrossSection(new_var,i));
    rndCS.push_back(h1);
    delete fResVaried;
    fResVaried=NULL;
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
  if (fResVaried) fResVaried->Write("respVaried");
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
  TObjString infoNIterDetRes(Form("%d",fNItersDetRes));
  infoNIterDetRes.Write("infoNIterDetRes");
  TObjString infoNIterFSR(Form("%d",fNItersFSR));
  infoNIterFSR.Write("infoNIterFSR");
  TObjString infoLumi(Form("%lf",fLumi));
  infoLumi.Write("lumi");
  TObjString timeTag(DayAndTimeTag(0));
  timeTag.Write("timeTag");
  TObjString explain("productionTime");
  explain.Write(timeTag.String());

  fout.Close();

  std::cout << "file <" << fout.GetName() << "> created\n";
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
  if (fDetRes) fDetRes->UseOverflow(false);
  if (fFSRRes) fFSRRes->UseOverflow(false);
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
  fResVaried= loadRooUnfoldResponse(fin, "respVaried","respVaried"+fTag,0);
  if (fResVaried) fResVaried->UseOverflow(false);
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

  int noNIter=0;
  TObjString *infoNIterDetRes=(TObjString*)fin.Get("infoNIterDetRes");
  if (!infoNIterDetRes) { noNIter=1; HERE("infoNIterDetRes is null"); }
  else fNItersDetRes= atoi(infoNIterDetRes->String().Data());
  TObjString *infoNIterFSR=(TObjString*)fin.Get("infoNIterFSR");
  if (!infoNIterFSR) { noNIter+=2; HERE("infoNIterFSR is null"); }
  else fNItersFSR= atoi(infoNIterFSR->String().Data());
  if (noNIter) {
    TObjString *infoNIter=(TObjString*)fin.Get("infoNIter");
    if (!infoNIter) { res=0; HERE("infoNIter is null"); }
    else {
      fNItersDetRes= atoi(infoNIter->String().Data());
      fNItersFSR= fNItersDetRes;
      delete infoNIter;
    }
  }

  TObjString *infoLumi= (TObjString*)fin.Get("lumi");
  if (!infoLumi) { res=0; HERE("infoLumi is null"); }
  else fLumi= atof(infoLumi->String().Data());

  if (res) {
    delete infoName;
    //delete infoTag;
    delete infoVar;
    delete infoCSType;
    if (infoNIterDetRes) delete infoNIterDetRes;
    if (infoNIterFSR) delete infoNIterFSR;
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
				       double lumiA, double lumiB,
				       TVersion_t setVersion) :
  fCSa(setName + TString("a"), setTag + TString("a"), _csPostFsrFullSp, setVersion),
  fCSb(setName + TString("b"), setTag + TString("b"), _csPostFsrFullSp, setVersion),
  fTag(setTag),
  fVersion(setVersion),
  fh1BkgV(), fBkgWeightUnc(),
  fh1Bkg(NULL),
  fh1PostFSR(NULL), fh1PreFSR(NULL),
  fh1CS(NULL),
  fh1Theory(NULL),
  fFSRBayes(NULL)
{
  fCSa.lumi(lumiA);
  fCSb.lumi(lumiB);
}

// --------------------------------------------------------------

void MuonCrossSection_t::h1Bkg_approximate(const TH1D *h1, int clearVec)
{

  if (clearVec) {
    fh1BkgV.clear();
    fBkgWeightUnc.Clear();
  }

  double loc_lumiTot= fCSa.lumi() + fCSb.lumi();
  TH1D *h1bkgA= fCSa.copy(h1);
  h1bkgA->Scale( fCSa.lumi()/loc_lumiTot );
  fCSa.h1Bkg(h1bkgA);
  TH1D *h1bkgB= fCSb.copy(h1);
  h1bkgB->Scale( fCSb.lumi()/loc_lumiTot );
  fCSb.h1Bkg(h1bkgB);

  fh1Bkg=fCSa.copy(h1, "h1BkgTot", fTag);
}

// --------------------------------------------------------------

int MuonCrossSection_t::setBkgV(const std::vector<TH1D*> &setBkg,
				 const std::vector<double> &weightUnc)
{
  if (fh1Bkg) { delete fh1Bkg; fh1Bkg=NULL; }
  fh1BkgV.clear();
  fBkgWeightUnc.Clear();
  if (setBkg.size()==0) {
    std::cout << "setBkgV: size is 0\n";
    return 0;
  }
  //TH1D *h1BkgSum= fCSa.copy(setBkg[0],fTag);
  for (unsigned int i=0; i<setBkg.size(); i++) {
    TH1D *h1=fCSa.copy(setBkg[i],fTag);
    fh1BkgV.push_back(h1);
    //if (i>0) h1BkgSum->Add(h1);
  }
  fBkgWeightUnc.ResizeTo(weightUnc.size());
  for (unsigned int i=0; i<weightUnc.size(); i++) {
    fBkgWeightUnc[i]= weightUnc[i];
  }
  //h1Bkg(h1BkgSum,0);
  return recalcBkg(NULL);
}

// --------------------------------------------------------------

int MuonCrossSection_t::recalcBkg(const std::vector<double> *weights,
				  TH1D **h1bkg_fromMC_out,
				  TH1D **h1bkg_fromData_out)
{
  if (fh1Bkg) { delete fh1Bkg; fh1Bkg=NULL; }
  if (fh1BkgV.size()==0) {
    std::cout << "recalcBkg: fh1BkgV.size=0\n";
    return 0;
  }
  if (weights && (fh1BkgV.size()!=weights->size())) {
    std::cout << "recalcBkg: fh1BkgV.size=" << fh1BkgV.size()
	      << ", weights->size=" << weights->size() << "\n";
    return 0;
  }
  TH1D *h1Bkg_fromMC= fCSa.copy(fh1BkgV[_bkgZZ],fTag);
  h1Bkg_fromMC->Reset();
  TH1D *h1Bkg_fromData= fCSa.copy(h1Bkg_fromMC,fTag);
  double wZZ=(weights) ? weights->at(_bkgZZ) : 1.;
  double wWZ=(weights) ? weights->at(_bkgWZ) : 1.;
  //printBin(h1Bkg_fromMC,0,0); std::cout << " + ";
  //printBin(fh1BkgV[_bkgZZ],0,0); std::cout << " + ";
  //printBin(fh1BkgV[_bkgWZ],0,1);
  h1Bkg_fromMC->Add(fh1BkgV[_bkgZZ],wZZ);
  h1Bkg_fromMC->Add(fh1BkgV[_bkgWZ],wWZ);
  //std::cout << "   = "; printBin(h1Bkg_fromMC,0,1);
  for (unsigned int i=0; i<fh1BkgV.size(); i++) {
    if ((i==_bkgZZ) || (i==_bkgWZ)) continue;
    double w=(weights) ? weights->at(i) : 1.;
    h1Bkg_fromData->Add(fh1BkgV[i],w);
  }

  double loc_lumiTot= fCSa.lumi() + fCSb.lumi();
  if (fVersion==_verMuApproved) {
    if (!fCSa.h1Rho() || !fCSb.h1Rho()) {
      std::cout << "MuonCrossSection_t::recalcBkg -- fCSv.h1Rho() is lacking\n";
      return 0;
    }
    TH1D *h1bkgA= fCSa.copy(h1Bkg_fromMC);
    h1bkgA->Multiply(fCSa.h1Rho());
    h1bkgA->Add(h1Bkg_fromData);
    h1bkgA->Scale( fCSa.lumi()/loc_lumiTot );
    fCSa.h1Bkg(h1bkgA);
    TH1D *h1bkgB= fCSb.copy(h1Bkg_fromMC);
    h1bkgB->Multiply(fCSb.h1Rho());
    h1bkgB->Add(h1Bkg_fromData);
    h1bkgB->Scale( fCSb.lumi()/loc_lumiTot );
    fCSb.h1Bkg(h1bkgB);
    // set the main bkg histogram
    if (fh1Bkg) delete fh1Bkg;
    fh1Bkg = fCSa.copy(h1Bkg_fromMC,"h1BkgSum_",fTag);
    fh1Bkg->Add(h1Bkg_fromData);
    // clean-up
    delete h1bkgA;
    delete h1bkgB;
  }
  else {
    // set the main bkg histogram
    if (fh1Bkg) delete fh1Bkg;
    fh1Bkg= fCSa.copy(h1Bkg_fromMC,"h1BkgSum_",fTag);
    fh1Bkg->Add(h1Bkg_fromData);
    TH1D *h1bkgA= fCSa.copy(fh1Bkg);
    h1bkgA->Scale( fCSa.lumi()/loc_lumiTot );
    fCSa.h1Bkg(h1bkgA);
    TH1D *h1bkgB= fCSb.copy(fh1Bkg);
    h1bkgB->Scale( fCSb.lumi()/loc_lumiTot );
    fCSb.h1Bkg(h1bkgB);
    // cleanup
    delete h1bkgA;
    delete h1bkgB;
  }

  if (h1bkg_fromMC_out) *h1bkg_fromMC_out= h1Bkg_fromMC;
  else delete h1Bkg_fromMC;
  if (h1bkg_fromData_out) *h1bkg_fromData_out= h1Bkg_fromData;
  else delete h1Bkg_fromData;

  //h1Bkg(h1BkgSum,0); // the latest muon results uses partially scaled Bkg
  return 1;
}

// --------------------------------------------------------------

TH1D* MuonCrossSection_t::calcCrossSection(int removeNegativeSignal)
{
  if (!fh1Bkg) {
    std::cout << "MuonCrossSection_t::calcCrossSection: object is not ready\n";
    return NULL;
  }
  TH1D *h1a= fCSa.calcCrossSection(removeNegativeSignal);
  TH1D *h1b= fCSb.calcCrossSection(removeNegativeSignal);
  /*
  TH1D *h1postFSR= fCSa.copy(fCSa.calcCrossSection(),"h1postFSR",fTag);
  TH1D *h1CSb= fCSb.copy(fCSb.calcCrossSection());
  if (!h1postFSR || !h1CSb) {
    std::cout << "MuonCrossSection_t::calcCrossSection - failure\n";
    return NULL;
  }
  h1postFSR->Scale( fCSa.lumi() );
  h1CSb->Scale( fCSb.lumi() );
  h1postFSR->Add(h1CSb);

  if (fFSRBayes) delete fFSRBayes;
  fFSRBayes= new RooUnfoldBayes( fCSa.fFSRRes, h1postFSR, fCSa.fNItersFSR, false);
  TH1D* h1CS_raw= (TH1D*) fFSRBayes->Hreco()->Clone("h1preFSR_" + fTag);
  h1CS_raw->SetDirectory(0);
  h1CS_raw->SetTitle("h1preFSR_" + fTag +
		     TString(";") + fh1Bkg->GetXaxis()->GetTitle()
		     + TString("; pre-FSR full space cross section"));
  double lumiTot= fCSa.lumi() + fCSb.lumi();
  h1CS_raw->Scale(1/lumiTot);
  */
  fh1CS= this->calcPreFsrCS_sumAB(h1a,h1b,fTag);
  return fh1CS;
}

// --------------------------------------------------------------

TH1D* MuonCrossSection_t::calcCrossSection(TVaried_t new_var, int idx,
					   const TH1D *h1_setVarA,
					   const TH1D *h1_setVarB,
					   int removeNegativeSignal)
{
  if (!fh1Bkg) {
    std::cout << "MuonCrossSection_t::calcCrossSection: object is not ready\n";
    return NULL;
  }
  if (!h1_setVarA || !h1_setVarB) {
    std::cout << "MuonCrossSection_t::calcCrossSection(TVaried): the provide "
	      << "histograms are null\n";
    return NULL;
  }
  TString idxStr= Form("_%d",idx);
  if (!fCSa.fh1Varied) {
    TH1D *h1a= fCSa.getVariedHisto(new_var);
    fCSa.fh1Varied= cloneHisto(h1a,h1a->GetName()+idxStr,h1a->GetTitle()+idxStr);
  }
  if (!fCSb.fh1Varied) {
    TH1D *h1b= fCSb.getVariedHisto(new_var);
    fCSb.fh1Varied= cloneHisto(h1b,h1b->GetName()+idxStr,h1b->GetTitle()+idxStr);
  }
  if (!fCSa.fh1Varied || !fCSb.fh1Varied) {
    std::cout << "MuonCrossSection_t::calcCrossSection(TVaried): fh1Varied is null\n";

    return NULL;
  }
  copyContents(fCSa.fh1Varied, h1_setVarA);
  copyContents(fCSb.fh1Varied, h1_setVarB);
  TH1D *h1a= fCSa.calcCrossSection(new_var,idx,removeNegativeSignal);
  TH1D *h1b= fCSb.calcCrossSection(new_var,idx,removeNegativeSignal);
  fh1CS= this->calcPreFsrCS_sumAB(h1a,h1b,fTag+idxStr);
  return fh1CS;
}

// --------------------------------------------------------------

TCanvas* MuonCrossSection_t::plotCrossSection(TString canvName,int recalculate,
					      int removeNegativeSignal)
{
  TH1D *h1= NULL;
  if (recalculate && fh1CS) { delete fh1CS; fh1CS=NULL; }
  if (!fh1CS) this->calcCrossSection(removeNegativeSignal);
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

int MuonCrossSection_t::sampleRndVec(TVaried_t new_var, int sampleSize,
				     std::vector<TH1D*> &rndCS,
				     std::vector<TH1D*> *rndCSa_out,
				     std::vector<TH1D*> *rndCSb_out,
				     std::vector<TH1D*> *rndVarVec_out)
{

  if (new_var==_varAcc) {
    int recalc=0;
    if (fCSa.h1EffAcc()!=NULL) { fCSa.removeEffAcc(); recalc=1; }
    if (fCSb.h1EffAcc()!=NULL) { fCSb.removeEffAcc(); recalc=1; }
    if (recalc) {
      if (!this->calcCrossSection()) {
	std::cout << "MuonCrossSection_t::sampleRndVec: recalculation failed\n";
	return 0;
      }
    }
  }

  int trackCSa=1;
  if (trackCSa) {
    TH1D *h1a=fCSa.copy(fCSa.h1UnfRhoEffAccCorr(),fCSa.h1UnfRhoEffAccCorr()->GetName()+TString("tmp"),fTag + TString("tmp"));
    h1a->SetTitle("UnfRhoEffAccCorr " + fCSa.tag() + " (postFSR cs)");
    h1a->GetXaxis()->SetNoExponent();
    h1a->GetXaxis()->SetMoreLogLabels();
    plotHisto(h1a, "muCS_a",1,1,"LPE");
  }
  int trackCSb=1;
  if (trackCSb) {
    TH1D *h1b=fCSb.copy(fCSb.h1UnfRhoEffAccCorr(),fCSb.h1UnfRhoEffAccCorr()->GetName()+TString("tmp"),fTag + TString("tmp"));
    h1b->SetTitle("UnfRhoEffAccCorr " + fCSb.tag() + " (postFSR cs)");
    h1b->GetXaxis()->SetNoExponent();
    h1b->GetXaxis()->SetMoreLogLabels();
    plotHisto(h1b, "muCS_b",1,1,"LPE");
  }

  int trackCS=1;
  if (trackCS) {
    plotHisto(fh1CS,"muCSVar",1,1,"LPE");
  }

  int res=1;
  rndCS.clear();
  std::vector<TH1D*> rndCSa, rndCSb;
  rndCS.reserve(sampleSize);
  rndCSa.reserve(sampleSize);
  rndCSb.reserve(sampleSize);
  if ((new_var==_varYield) ||
      (new_var==_varBkg) ||
      (new_var==_varBkgXS) ||
      (new_var==_varRho)) {
    if (new_var==_varRho) std::cout << "\n\tIt is not good to use _varRho int the muon channel\n";
    // Prepare post-FSR cross sections
    if ((new_var==_varYield) || (new_var==_varRho)) {
      res= (fCSa.sampleRndVec(new_var,sampleSize,rndCSa)
	    && fCSb.sampleRndVec(new_var,sampleSize,rndCSb)
	    ) ? 1:0;
    }
    else if (new_var==_varBkg) {
      std::vector<TH1D*> rndBkgVa, rndBkgVb;
      TH1D *h1rnd= NULL;
      double loc_lumiTot= fCSa.lumi() + fCSb.lumi();
      if (fVersion!=_verMuApproved) {
	for (int i=0; i<sampleSize; i++) {
	  TString useTag= Form("_rnd%d",i) + fTag;
	  h1rnd=cloneHisto(fh1Bkg,
			   "h1rndBkg" + useTag,
			   fh1Bkg->GetTitle());
	  int nonNegative=0;
	  randomizeWithinErr(fh1Bkg,h1rnd,nonNegative);
	  TH1D* h1rndA= fCSa.copy(h1rnd);
	  h1rndA->Scale( fCSa.lumi()/loc_lumiTot );
	  TH1D* h1rndB= fCSb.copy(h1rnd);
	  h1rndB->Scale( fCSb.lumi()/loc_lumiTot );

	  if (rndVarVec_out) rndVarVec_out->push_back(h1rnd);
	  else delete h1rnd;
	  rndBkgVa.push_back(h1rndA);
	  rndBkgVb.push_back(h1rndB);
	}
      }
      else { // special case -- multiplication by rho is needed
	TH1D *h1bkg_fromMC=NULL, *h1bkg_fromData=NULL;
	if (!recalcBkg(NULL,&h1bkg_fromMC,&h1bkg_fromData)) {
	  std::cout << "failed to prepare bkg vectors\n";
	}
	TH1D *h1rndData=NULL;
	for (int i=0; i<sampleSize; i++) {
	  TString useTag= Form("_rnd%d",i) + fTag;
	  h1rnd= cloneHisto(h1bkg_fromMC,
			    "h1rndBkgMC"+useTag, fh1Bkg->GetTitle());
	  h1rndData=cloneHisto(h1bkg_fromData,
			       "h1rndBkgData" + useTag,fh1Bkg->GetTitle());
	  int nonNegative=0;
	  randomizeWithinErr(h1bkg_fromMC,h1rnd,nonNegative);
	  randomizeWithinErr(h1bkg_fromData,h1rndData,nonNegative);
	  TH1D *h1rndA= fCSa.copy(h1rnd);
	  h1rndA->Multiply(fCSa.h1Rho());
	  h1rndA->Add(h1bkg_fromData);
	  h1rndA->Scale( fCSa.lumi()/loc_lumiTot );
	  TH1D *h1rndB= fCSb.copy(h1rnd);
	  h1rndB->Multiply(fCSb.h1Rho());
	  h1rndB->Add(h1bkg_fromData);
	  h1rndB->Scale( fCSb.lumi()/loc_lumiTot );

	  delete h1rnd;
	  if (rndVarVec_out) {
	    std::cout << "sampleRndVec(varBkg) does not return full rndBkg\n";
	    rndVarVec_out->push_back(h1rndData);
	  }
	  else delete h1rndData;
	  rndBkgVa.push_back(h1rndA);
	  rndBkgVb.push_back(h1rndB);
	}
      }

      res= (fCSa.sampleRndVec(_varBkg, rndBkgVa, rndCSa) &&
	    fCSb.sampleRndVec(_varBkg, rndBkgVb, rndCSb)) ? 1:0;
    }
    else if (new_var==_varBkgXS) {
      TH1D *h1BkgCopy= cloneHisto(fh1Bkg, fh1Bkg->GetName() + TString("_copy"),
				  fh1Bkg->GetTitle());

      std::vector<TH1D*> rndBkgVa, rndBkgVb;
      TH1D *h1rnd= NULL;
      double loc_lumiTot= fCSa.lumi() + fCSb.lumi();
      std::vector<double> rndV;
      for (int i=0; i<fBkgWeightUnc.GetNoElements(); i++) rndV.push_back(1.);
      for (int i=0; i<sampleSize; i++) {
	TString useTag= Form("_rnd%d",i) + fTag;
	for (int i=0; i<fBkgWeightUnc.GetNoElements(); i++) {
	  rndV[i] = gRandom->Gaus(1.,fBkgWeightUnc[i]);
	}
	if (!this->recalcBkg(&rndV)) {
	  std::cout << "error recalculating the background\n";
	  return 0;
	}
	h1rnd=cloneHisto(fh1Bkg,
			 "h1rndBkg" + useTag,
			 fh1Bkg->GetTitle());
	TH1D* h1rndA= fCSa.copy(h1rnd);
	h1rndA->Scale( fCSa.lumi()/loc_lumiTot );
	TH1D* h1rndB= fCSb.copy(h1rnd);
	h1rndB->Scale( fCSb.lumi()/loc_lumiTot );

	if (rndVarVec_out) rndVarVec_out->push_back(h1rnd);
	else delete h1rnd;
	rndBkgVa.push_back(h1rndA);
	rndBkgVb.push_back(h1rndB);
      }

      res= (fCSa.sampleRndVec(_varBkg, rndBkgVa, rndCSa) &&
	    fCSb.sampleRndVec(_varBkg, rndBkgVb, rndCSb)) ? 1:0;
      res=res && copyContents(fh1Bkg, h1BkgCopy);
      delete h1BkgCopy;
    }
    else {
      std::cout << "case is not ready !!\n";
      return 0;
    }
  }

  else if ((new_var==_varEff) || (new_var==_varAcc) || (new_var==_varEffAcc)) {
    std::vector<TH1D*> rndVarVec;
    TH1D *h1rnd= NULL;
    TH1D *h1varied= fCSa.getVariedHisto(new_var);
    if (!h1varied) {
      std::cout << "failed to get the varied histo\n";
      return 0;
    }
    for (int i=0; i<sampleSize; i++) {
      TString useTag= Form("_rnd%d",i) + fTag;
      h1rnd=cloneHisto(h1varied,
		       "h1varied" + useTag,
		       h1varied->GetTitle() + useTag);
      int nonNegative=0;
      randomizeWithinErr(h1varied,h1rnd,nonNegative);
      rndVarVec.push_back(h1rnd);
      if (rndVarVec_out) rndVarVec_out->push_back(h1rnd);
    }

    res= (fCSa.sampleRndVec(new_var, rndVarVec, rndCSa) &&
	  fCSb.sampleRndVec(new_var, rndVarVec, rndCSb)) ? 1:0;
  }

  else if (new_var==_varDetRes) {
    std::vector<RooUnfoldResponse*> rndRespV;
    if (!fCSa.sampleRndResponse(new_var,sampleSize,rndCSa,&rndRespV) ||
	!fCSb.sampleRndRespVec(new_var,rndRespV,rndCSb)) {
      std::cout << " error in MuonCrossSection_t::sampleRndVec\n";
      return 0;
    }
  }

  /*
  else if (new_var==_varBkg) {
    std::vector<TH1D*> rndBkgV;
    TH1D *h1rnd= NULL;
    double lumiTot= fCSa.lumi() + fCSb.lumi();
    for (int i=0; i<sampleSize; i++) {
      TString useTag= Form("_rnd%d",i) + fTag;
      h1rnd=cloneHisto(fh1Bkg,
		       "h1rndBkg" + useTag,
		       fh1Bkg->GetTitle());
      int nonNegative=0;
      randomizeWithinErr(fh1Bkg,h1rnd,nonNegative);
      rndBkgV.push_back(h1rnd);
      fCSa.h1Bkg(h1rnd, fCSa.lumi()/lumiTot);
      fCSb.h1Bkg(h1rnd, fCSb.lumi()/lumiTot);
      TH1D *h1csA_raw= fCSa.calcCrossSection();
      TH1D *h1csA = fCSa.copy(h1csA_raw,
			      h1csA_raw->GetName() + useTag + "a",
			      fCSa.tag());
      TH1D *h1csB_raw= fCSb.calcCrossSection();
      TH1D* h1csB= fCSb.copy(h1csB_raw,
			     h1csB_raw->GetName() + useTag + "b",
			     fCSb.tag());
      TH1D *h1cs = this->calcPreFsrCS_sumAB(h1csA,h1csB,useTag);
      rndCSa.push_back(h1csA);
      rndCSb.push_back(h1csB);
      rndCS .push_back(h1cs);

      if (trackCSa) {
	h1csA->SetLineColor((i%9)+1);
	plotHistoSame(fCSa.copy(h1csA,h1csA->GetName()+TString("tmp"),fTag + TString("tmp")),"muCS_a","hist");
      }
      if (trackCSb) {
	h1csB->SetLineColor((i%9)+1);
	plotHistoSame(fCSb.copy(h1csB,h1csB->GetName()+TString("tmp"),fTag + TString("tmp")),"muCS_b","hist");
      }
      if (trackCS) {
	h1cs->SetLineColor((i%9)+1);
	plotHistoSame(h1cs,"muCSVar","hist");
      }
    }
    if (rndVarVec_out) (*rndVarVec_out)= rndBkgV;
  }
  */
  else if ((new_var==_varFSRRes) || (new_var==_varFSRRes_Poisson)) {
    double loc_lumiTot= fCSa.lumi() + fCSb.lumi();
    const TH1D *h1aCS= fCSa.h1UnfRhoEffAccCorr();
    TH1D* h1postFSRYield= cloneHisto( h1aCS, "h1postFSRYield", "postFSR yield" );
    h1postFSRYield->Scale( fCSa.lumi() );
    h1postFSRYield->Add( fCSb.h1UnfRhoEffAccCorr(), fCSb.lumi() );
    const RooUnfoldResponse *r= fCSa.fFSRRes;
    int poissonRnd=0;
    if (new_var==_varFSRRes_Poisson) poissonRnd=1;
    for (int i=0; i<sampleSize; i++) {
      TString respName= variedVarName(new_var) + Form("_%d",i) + fTag;
      RooUnfoldResponse *rRnd= randomizeWithinErr(r, respName, poissonRnd);
      if (!rRnd) {
	std::cout << "failed to get randomized response\n";
	return 0;
      }
      rRnd->UseOverflow(false); // ensure correct behavior
      if (fFSRBayes) { delete fFSRBayes; fFSRBayes=NULL; }
      TH1D* h1tmp= cloneHisto(h1postFSRYield,"h1postFSRYield" + respName, "h1postFSRYield"+respName);

      //fFSRBayes= new RooUnfoldBayes( rRnd, h1tmp, fCSa.fNItersFSR, false);
      //TH1D* h1preFSRYield=
      //cloneHisto(fFSRBayes->Hreco(), "h1preFSRYield_" + respName,
      //	   "h1preFSR_" + respName +
      //	   TString(";") + fh1Bkg->GetXaxis()->GetTitle()
      //	   +TString("; preFSR full space cs"), h1dummy);

      TH1D *h1preFSRYield= cloneHisto(h1postFSRYield,
				      "h1preFSRYield_" + respName,
				      "h1preFSR_" + respName +
			      TString(";") + fh1Bkg->GetXaxis()->GetTitle()
		   +TString("; preFSR full space cs"), h1dummy);
      if ((fVersion!=_verMu76X) && (fVersion!=_verMuApproved)) {
	//std::cout << "\n\taaaaa\n"; return 0;
	fFSRBayes= new RooUnfoldBayes( rRnd, h1tmp, fCSa.fNItersFSR, false);
	copyContents(h1preFSRYield, fFSRBayes->Hreco());
      }
      else {
	if (i%1000==0) std::cout << "\tdoing matrix inversion (2). i=" << i << "\n";
	RooUnfoldInvert *rooInv = new RooUnfoldInvert( rRnd, h1tmp );
	copyContents(h1preFSRYield, rooInv->Hreco());
	delete rooInv;
      }

      TH1D* h1cs= perMassBinWidth(h1preFSRYield,0);
      h1cs->Scale(1/loc_lumiTot);
      delete h1tmp;
      delete h1preFSRYield;
      rndCS.push_back(h1cs);

      if (trackCS) {
	h1cs->SetLineColor((i%9)+1);
	plotHistoSame(h1cs,"muCSVar","hist");
      }
    }
  }
  else {
    std::cout << "sampleRndVec is not ready for new_var="
	      << variedVarName(new_var) << "\n";
    res=0;
  }


  // calculate final cross section, if it was not the case of varFSRRes*
  if (res && (new_var!=_varFSRRes) && (new_var!=_varFSRRes_Poisson)) {
    for (int i=0; i<sampleSize; i++) {
      TString useTag=Form("_rnd%d",i) + fTag;
      std::cout << "iRnd=" << i << "\n";
      TH1D *h1cs= this->calcPreFsrCS_sumAB(rndCSa[i],rndCSb[i],useTag);
      rndCS.push_back(h1cs);

      if (trackCSa) {
	rndCSa[i]->SetLineColor((i%9)+1);
	plotHistoSame(fCSa.copy(rndCSa[i],rndCSa[i]->GetName()+TString("tmp"),fTag + TString("tmp")),"muCS_a","hist");
      }
      if (trackCSb) {
	rndCSb[i]->SetLineColor((i%9)+1);
	plotHistoSame(fCSb.copy(rndCSb[i],rndCSb[i]->GetName()+TString("tmp"),fTag + TString("tmp")),"muCS_b","hist");
      }
      if (trackCS) {
	h1cs->SetLineColor((i%9)+1);
	plotHistoSame(h1cs,"muCSVar","hist");
      }
    }
  }


  if (rndCSa_out) (*rndCSa_out) = rndCSa;
  if (rndCSb_out) (*rndCSb_out) = rndCSb;
  if (!res) std::cout << "Error in MuonCrossSection_t::sampleRndVec\n";
  return res;
}

// --------------------------------------------------------------

int MuonCrossSection_t::sampleRndVec(TVaried_t new_var, int sampleSize,
				     const RndVecInfo_t &info,
				     std::vector<TH1D*> &rndCS,
				     std::vector<TH1D*> *rndCSa_out,
				     std::vector<TH1D*> *rndCSb_out,
				     std::vector<TH1D*> *rndVarVec1_out,
				     std::vector<TH1D*> *rndVarVec2_out)
{

  if (info.histoNameBaseV2().Length()==0) {
    std::cout << "MuonCrossSection_t::sampleRndVec requires two histoNameBases\n";
    return 0;
  }

  if (new_var==_varAcc) {
    int recalc=0;
    if (fCSa.h1EffAcc()!=NULL) { fCSa.removeEffAcc(); recalc=1; }
    if (fCSb.h1EffAcc()!=NULL) { fCSb.removeEffAcc(); recalc=1; }
    if (recalc) {
      if (!this->calcCrossSection()) {
	std::cout << "MuonCrossSection_t::sampleRndVec: recalculation failed\n";
	return 0;
      }
    }
  }

  int trackCSa=1;
  if (trackCSa) {
    TH1D *h1a=fCSa.copy(fCSa.h1UnfRhoEffAccCorr(),fCSa.h1UnfRhoEffAccCorr()->GetName()+TString("tmp"),fTag + TString("tmp"));
    h1a->SetTitle("UnfRhoEffAccCorr " + fCSa.tag() + " (postFSR cs)");
    h1a->GetXaxis()->SetNoExponent();
    h1a->GetXaxis()->SetMoreLogLabels();
    plotHisto(h1a, "muCS_a",1,1,"LPE");
  }
  int trackCSb=1;
  if (trackCSb) {
    TH1D *h1b=fCSb.copy(fCSb.h1UnfRhoEffAccCorr(),fCSb.h1UnfRhoEffAccCorr()->GetName()+TString("tmp"),fTag + TString("tmp"));
    h1b->SetTitle("UnfRhoEffAccCorr " + fCSb.tag() + " (postFSR cs)");
    h1b->GetXaxis()->SetNoExponent();
    h1b->GetXaxis()->SetMoreLogLabels();
    plotHisto(h1b, "muCS_b",1,1,"LPE");
  }

  int trackCS=1;
  if (trackCS) {
    plotHisto(fh1CS,"muCSVar",1,1,"LPE");
  }

  int res=1;
  rndCS.clear();
  std::vector<TH1D*> rndVarA,rndVarB;
  std::vector<TH1D*> rndCSa, rndCSb;
  rndVarA.reserve(sampleSize);
  rndVarB.reserve(sampleSize);
  rndCS.reserve(sampleSize);
  rndCSa.reserve(sampleSize);
  rndCSb.reserve(sampleSize);

  res= (info.loadHistos(sampleSize,rndVarA,&rndVarB) &&
	fCSa.sampleRndVec(new_var, rndVarA, rndCSa) &&
	fCSb.sampleRndVec(new_var, rndVarB, rndCSb)) ? 1 : 0;

  if (rndVarVec1_out) *rndVarVec1_out= rndVarA;
  else clearVec(rndVarA);
  if (rndVarVec2_out) *rndVarVec2_out= rndVarB;
  else clearVec(rndVarB);

  // calculate final cross section, if it was not the case of varFSRRes*
  if (res) {
    for (int i=0; i<sampleSize; i++) {
      TString useTag=Form("_rnd%d",i) + fTag;
      std::cout << "iRnd=" << i << "\n";
      TH1D *h1cs= this->calcPreFsrCS_sumAB(rndCSa[i],rndCSb[i],useTag);
      rndCS.push_back(h1cs);

      if (trackCSa) {
	rndCSa[i]->SetLineColor((i%9)+1);
	plotHistoSame(fCSa.copy(rndCSa[i],rndCSa[i]->GetName()+TString("tmp"),fTag + TString("tmp")),"muCS_a","hist");
      }
      if (trackCSb) {
	rndCSb[i]->SetLineColor((i%9)+1);
	plotHistoSame(fCSb.copy(rndCSb[i],rndCSb[i]->GetName()+TString("tmp"),fTag + TString("tmp")),"muCS_b","hist");
      }
      if (trackCS) {
	h1cs->SetLineColor((i%9)+1);
	plotHistoSame(h1cs,"muCSVar","hist");
      }
    }
  }


  if (rndCSa_out) (*rndCSa_out) = rndCSa;
  if (rndCSb_out) (*rndCSb_out) = rndCSb;
  if (!res) std::cout << "Error in MuonCrossSection_t::sampleRndVec(RndVecInfo)\n";
  return res;
}

// --------------------------------------------------------------

int MuonCrossSection_t::deriveCov(const std::vector<TH1D*> &rndCS,
				  TVaried_t var,
				  TH1D **h1avgCS_out, TH2D **h2cov_out)
{
  int res= deriveCovariance(rndCS,
			    variedVarName(var) + fTag,
			    "avg " + csTypeName(_csPreFsrFullSp),
			    h1avgCS_out, h2cov_out);
  if (!res) std::cout << "error in MuonCrossSection_t::deriveCov\n";
  return res;
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
  fBkgWeightUnc.Write("bkgWeightUnc");
  fout.mkdir("bkg");
  fout.cd("bkg");
  for (unsigned int i=0; i<fh1BkgV.size(); i++) fh1BkgV[i]->Write();
  fout.cd();
  if (fh1Bkg) fh1Bkg->Write("h1Bkg");
  if (fh1PostFSR) fh1PostFSR->Write("h1PostFSR");
  if (fh1PreFSR) fh1PreFSR->Write("h1PreFSR");
  if (fh1CS) fh1CS->Write("h1CS");
  if (fh1Theory) fh1Theory->Write("h1Theory");
  //HERE("save fsrBayes");
  //if (fFSRBayes) fFSRBayes->Write("fsrBayes"); // crashes

  TObjString infoTag(fTag);
  infoTag.Write("tag");
  TString strVersion(Form("%d ",int(fVersion)) + versionName(fVersion));
  TObjString infoVersion(strVersion);
  infoVersion.Write("version");
  TString strLumiTot(Form("%8.6lf",this->lumiTot()));
  TObjString lumi(strLumiTot);
  lumi.Write("lumiTot");
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
  fBkgWeightUnc.Read("bkgWeightUnc");
  if (fBkgWeightUnc.GetNoElements()>0) {
    fh1BkgV.reserve(fBkgWeightUnc.GetNoElements());
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
  fh1PostFSR= loadHisto(fin, "h1PostFSR", "h1PostFSR" + fTag, 1, h1dummy);
  fh1PreFSR= loadHisto(fin, "h1PreFSR", "h1PreFSR" + fTag, 1, h1dummy);
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

TH1D* MuonCrossSection_t::calcPreFsrCS_sumAB(const TH1D *h1a,
					     const TH1D *h1b,
					     TString useTag)
{
  //HERE("entered calcPreFsrCS_sumAB");
  std::cout << "entered calcPreFSRCS_sumAB\n";
  if (!h1a || !h1b) {
    std::cout << "calcPreFsrCS_sumAB : null ptrs\n";
    return NULL;
  }
  if (fh1PostFSR) { delete fh1PostFSR; }
  fh1PostFSR= fCSa.copy(h1a,"h1PostFSR",fTag);
  fh1PostFSR->Scale( fCSa.lumi() );
  fh1PostFSR->Add(h1b, fCSb.lumi());
  //printHisto(h1postFSR);

  if (fFSRBayes) { delete fFSRBayes; fFSRBayes=NULL; }
  //fFSRBayes= new RooUnfoldBayes( fCSa.fFSRRes, h1postFSR, fCSa.fNItersFSR, false);
  //TH1D* h1CS_raw= (TH1D*) fFSRBayes->Hreco()->Clone("h1preFSR_" + useTag);

  if (fh1PreFSR) { delete fh1PreFSR; }
  fh1PreFSR = cloneHisto(fh1PostFSR, "h1Unf_loc","h1Unf_loc");
  if ((fVersion!=_verMu76X) && (fVersion!=_verMuApproved)) {
    fFSRBayes= new RooUnfoldBayes( fCSa.fFSRRes, fh1PostFSR, fCSa.fNItersFSR, false);
    copyContents(fh1PreFSR,fFSRBayes->Hreco());
  }
  else {
    std::cout << "\tdoing matrix inversion\n";
    RooUnfoldInvert *rooInv= new RooUnfoldInvert( fCSa.fFSRRes, fh1PostFSR, "rooInv" );
    if (0) {
      std::cout << "check whether fh1PostFSR has changed\n";
      printRatio(fh1PreFSR,fh1PostFSR);
    }
    copyContents(fh1PreFSR, rooInv->Hreco());
    delete rooInv;
  }

  fh1PreFSR->SetDirectory(0);
  fh1PreFSR->SetTitle("h1preFSR_" + useTag +
		     TString(";") + fh1Bkg->GetXaxis()->GetTitle()
		     + TString("; pre-FSR full space cross section"));
  TH1D *h1CS_tmp = perMassBinWidth(fh1PreFSR,0);
  h1CS_tmp->Scale(1/this->lumiTot());
  //printHisto(h1CS_tmp);
  return h1CS_tmp;
}

// --------------------------------------------------------------
// --------------------------------------------------------------
// --------------------------------------------------------------
// --------------------------------------------------------------
