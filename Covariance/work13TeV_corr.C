#include "inputs.h"
#include "crossSection.h"
#include "Blue.h"
#include "analyseBLUEResult.h"
//#include "work13TeV.C" // making it obsolete
#include "DYbinning.h"

// -------------------------------------------------------
// -------------------------------------------------------

struct CovStruct_t {
  typedef enum { _varCov_none=0, _varCov_statYield, _varCov_statNonYield,
		 _varCov_systAcc, _varCov_systNonAcc } TCovPart_t;

  TMatrixD covStat_Yield, covStat_nonYield;
  TMatrixD covSyst_Acc, covSyst_nonAcc;
public:
  CovStruct_t(TMatrixD &M) :
    covStat_Yield(TMatrixD::kZero,M), covStat_nonYield(TMatrixD::kZero,M),
    covSyst_Acc(TMatrixD::kZero,M), covSyst_nonAcc(TMatrixD::kZero,M) {}

  CovStruct_t(const TMatrixD &set_covStatYield,
	      const TMatrixD &set_covStatNonYield,
	      const TMatrixD &set_covSystAcc,
	      const TMatrixD &set_covSystNonAcc) :
    covStat_Yield(set_covStatYield),
    covStat_nonYield(set_covStatNonYield),
    covSyst_Acc(set_covSystAcc),
    covSyst_nonAcc(set_covSystNonAcc)
  {}

  CovStruct_t(const CovStruct_t &s) :
    covStat_Yield(s.covStat_Yield), covStat_nonYield(s.covStat_nonYield),
    covSyst_Acc(s.covSyst_Acc), covSyst_nonAcc(s.covSyst_nonAcc)
  {}

  TMatrixD& editStatYield() { return covStat_Yield; }
  TMatrixD& editStatNonYield() { return covStat_nonYield; }
  TMatrixD& editSystAcc() { return covSyst_Acc; }
  TMatrixD& editSystNonAcc() { return covSyst_nonAcc; }

  TMatrixD getStatYield() const { return covStat_Yield; }
  TMatrixD getNonStatYield() const { return covStat_nonYield; }
  TMatrixD getSystAcc() const { return covSyst_Acc; }
  TMatrixD getSystNonAcc() const { return covSyst_nonAcc; }

  void Zero() {
    covStat_Yield.Zero(); covStat_nonYield.Zero();
    covSyst_Acc.Zero(); covSyst_nonAcc.Zero();
  }

  void ReduceCorrelations(double factor) {
    covStat_Yield= reduceCorrelations(covStat_Yield,factor);
    covStat_nonYield= reduceCorrelations(covStat_nonYield, factor);
    covSyst_Acc= reduceCorrelations(covSyst_Acc,factor);
    covSyst_nonAcc= reduceCorrelations(covSyst_nonAcc,factor);
  }

  TMatrixD Sum() const
  { return (covStat_Yield + covStat_nonYield + covSyst_Acc + covSyst_nonAcc); }

  TMatrixD GetPart(int i) const { return GetPart(TCovPart_t(i)); }

  int SetPart(int i, const TMatrixD &m)
  { return SetPart(TCovPart_t(i),m); }

  TString GetPartName(int i) const {
    TString name="UNKNOWN";
    switch(i) {
    case 0: name="none"; break;
    case 1: name="statYield"; break;
    case 2: name="statNonYield"; break;
    case 3: name="systAcc"; break;
    case 4: name="systNonAcc"; break;
    default: ;
    }
    return name;
  }

  TMatrixD GetPart(TCovPart_t p) const
  {
    TMatrixD a(covStat_Yield);
    switch (p) {
    case _varCov_none: a.Zero(); break;
    case _varCov_statYield: a=covStat_Yield; break;
    case _varCov_statNonYield: a=covStat_nonYield; break;
    case _varCov_systAcc: a=covSyst_Acc; break;
    case _varCov_systNonAcc: a=covSyst_nonAcc; break;
    default: a.Zero();
    }
    return a;
  }

  int SetPart(TCovPart_t p, const TMatrixD &m)
  {
    int res=0;
    switch (p) {
    case _varCov_none: break;
    case _varCov_statYield: covStat_Yield=m; res=1; break;
    case _varCov_statNonYield: covStat_nonYield=m; res=1; break;
    case _varCov_systAcc: covSyst_Acc=m; res=1; break;
    case _varCov_systNonAcc: covSyst_nonAcc=m; res=1; break;
    default: res=0;
    }
    return res;
  }

  int ZrangeCov(int idxMin, int idxMax_p1, const TH1D *h1BinDef,
		std::vector<double> &cov, int printNumbers=0) const
  {
    if (!h1BinDef) {
      std::cout << "CovStruct_t::ZrangeCov: h1BinDef is null\n";
      return 0;
    }
    cov.clear();
    cov.push_back(sumCov(covStat_Yield,idxMin,idxMax_p1,h1BinDef,printNumbers));
    cov.push_back(sumCov(covStat_nonYield,idxMin,idxMax_p1,h1BinDef));
    cov.push_back(sumCov(covSyst_Acc,idxMin,idxMax_p1,h1BinDef));
    cov.push_back(sumCov(covSyst_nonAcc,idxMin,idxMax_p1,h1BinDef));
    return 1;
  }

  int PrintZRangeUnc(TString label, const TH1D *h1cs_perMBW,
		     int printNumbers=0) const
  {
    int idxMin=-1, idxMax_p1=-1;
    if (!DYtools::ZmassRange(60,120,idxMin,idxMax_p1)) {
      std::cout << "in CovStruct_t::ZrangeCov\n";
      return 0;
    }
    std::cout << "ZmassRange: " << idxMin << " .. " << idxMax_p1 << "\n";
    std::vector<double> zCov;
    TH1D *h1cs_abs= timesMassBinWidth(h1cs_perMBW);
    if (!this->ZrangeCov(idxMin,idxMax_p1,h1cs_abs, zCov,printNumbers)) {
      std::cout << "CovStruct_t::PrintZRangeUnc: failed to get the ZrangeCov\n";
      return 0;
    }
    std::cout << "integrated uncertainties:\n";
    //printHisto(h1cs_abs);
    std::cout << "  - " << label << ": "
	      << h1cs_abs->Integral(idxMin+1,idxMax_p1) << " ";
    for (int ii=0; ii<4; ii++) {
      std::cout << "+- " << zCov[ii] << " (" << this->GetPartName(ii+1) << ")";
    }
    std::cout << "\n";
    return 1;
  }

  void Plot(TString tag, TH1D *h1binning) {
    plotCovCorr(covStat_Yield,h1binning,"h2CovStatYield_"+tag,"c2CovStatYield_"+tag);
    plotCovCorr(covStat_nonYield,h1binning,"h2CovStatNonYield_"+tag,"c2CovStatNonYield_"+tag);
    plotCovCorr(covSyst_Acc,h1binning,"h2CovSystAcc_"+tag,"c2CovSystAcc_"+tag);
    plotCovCorr(covSyst_nonAcc,h1binning,"h2CovSystNonAcc_"+tag,"c2CovSystNonAcc_"+tag);
  }

  void PlotUnc(TString tag, TH1D *h1binning, TH1D *h1centralVal=NULL,
	       double yrangeMin=0, double yrangeMax=0) {
    TH1D *h1statYield= uncFromCov(covStat_Yield,"h1statYield_"+tag,
				  h1binning,h1centralVal,1);
    TH1D *h1statNonYield= uncFromCov(covStat_nonYield,"h1statNonYield_"+tag,
				     h1binning,h1centralVal,1);
    TH1D *h1systAcc= uncFromCov(covSyst_Acc,"h1systAcc_"+tag,
				h1binning,h1centralVal,1);
    TH1D *h1systNonAcc= uncFromCov(covSyst_nonAcc,"h1systNonAcc_"+tag,
				   h1binning,h1centralVal,1);
    TH1D *h1tot= cloneHisto(h1statYield,"h1tot_"+tag,"h1tot_"+tag);
    addInQuadrature(h1tot, h1statNonYield);
    addInQuadrature(h1tot, h1systAcc);
    addInQuadrature(h1tot, h1systNonAcc);
    hsBlack.SetStyle(h1statYield);
    hsColor46.SetStyle(h1statNonYield);
    hsGreen.SetStyle(h1systAcc);
    hsBlue.SetStyle(h1systNonAcc);
    HistoStyle_t(kRed,7,2,0.8,2.0).SetStyle(h1tot);
    if ((yrangeMin!=0) || (yrangeMax!=0)) {
      h1statYield->GetYaxis()->SetRangeUser(yrangeMin,yrangeMax);
    }
    TString cName="cCovStructUnc_"+tag;
    logAxis(h1statYield,1+8);
    plotHisto(h1statYield,cName,1,1,"LP","stat Yield");
    plotHistoSame(h1statNonYield,cName,"LP","stat nonYield");
    plotHistoSame(h1systAcc,cName,"LP","syst Acc");
    plotHistoSame(h1systNonAcc,cName,"LP","syst nonAcc");
    plotHistoSame(h1tot,cName,"LP","total");
  }

};

// -------------------------------------------------------
// -------------------------------------------------------

TMatrixD constructEMAccCov(const CovStruct_t &eeCovS, const CovStruct_t &mmCovS,
			   double corrFactor=1.);
TMatrixD constructMeasCov(const CovStruct_t &eeCovS, const CovStruct_t &mmCovS,
			  CovStruct_t::TCovPart_t covFlag, double accCorrFactor,
			  BLUEResult_t *blue);

int loadUncData(int isEE,
		const finfomap_t &fNames,
		const finfomap_t &h1names,
		const fvarweightmap_t &relUncW,
		const TH1D *h1binning,
		const TH1D *h1centralVal, // if the uncertainty is relative
		std::vector<TH1D*> &h1uncV,
		TString tag);

int adjustMMUnc(const finfomap_t &mmCovFNames, // needed for first value
		std::vector<TMatrixD> &mmCovV,
		CovStruct_t &covM, int showChangedCov,
		int plotCmp);

int adjustEEUnc(const finfomap_t &eeCovFNames, // needed for first value
		std::vector<TMatrixD> &eeCovV,
		CovStruct_t &eeCovM, int showChangedCov,
		int plotCmp);

int aggregateUnc(const finfomap_t &fnames, // needed for first value
		 const std::vector<TMatrixD> &covStatV,
		 const std::vector<TH1D*> &h1SystV,
		 CovStruct_t &covM);

int adjustUnc(int isEE,
	      const finfomap_t &covFNames,
	      std::vector<TMatrixD> &covV,
	      const finfomap_t &uncFiles,
	      const finfomap_t &uncStat,
	      const finfomap_t &uncSyst,
	      const fvarweightmap_t &relUncW,
	      const TH1D *h1binning,
	      const TH1D *h1central,
	      int plotCmp, // 1 - plot uncertainties, 2 - and covariances
	      int change_covV, // 1 - change values, 2 - verify change by plot
	      TString tag,
	      CovStruct_t &covM);

int compareUnc(const std::vector<TString> &labelV,
	       const std::vector<TH1D*> &h1V_target,
	       std::vector<TMatrixD> &covV,
	       int plotCmp,
	       int change_covV,
	       TString tag);

// -------------------------------------------------------
/*
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
*/
// -------------------------------------------------------

const int useUncorrYieldUnc=0;
TString eeCSFName="cs_DYee_13TeV_El3.root";
TString eeCSH1Name="h1PreFSRCS";
TString mmCSFName="cs_DYmm_13TeVMuApproved_cs.root";
TString mmCSH1Name="h1CS";

// -------------------------------------------------------

void work13TeV_corr(int printCanvases=0, int corrCase=1, int includeLumiUnc=0,
		    int plotChangedCov=0,
		    int excludeSyst=0, // 1 - all, 2 - Acc only
		    TString saveTag="")
{
  closeCanvases(5);
  const int changeNames=1;

  std::string showCanvs;
  //showCanvs+=" plotChannelInputCov";
 // whether relative uncertainty in Acc in the channels is similar
  //showCanvs+=" plotChannelRelSystAccUnc";
  // compare randomized vs provided uncertainties
  //showCanvs+=" plotAdjUnc";
  //showCanvs+=" plotAdjCov";
  //showCanvs+=" plotInputContributedCovs"; // 4 canvases in each channel
  showCanvs+= " plotFinalCovByType"; // 4 canvases in the combined channel

  std::string showCombinationCanvases;
  //showCombinationCanvases="ALL";

  std::vector<int> eeCovIdx, mmCovIdx;
  std::vector<TString> eeOldFName, eeNewFName; // for file name changing
  std::vector<TString> mmOldFName, mmNewFName; // for file name changing
  if (1) {
    addToVector(eeCovIdx,7, _varYield, _varBkg, _varDetRes, _varEff,
		_varRhoFile, _varAcc, _varFSRRes);
    addToVector(mmCovIdx,7, _varYield, _varBkg, _varDetRes, _varEff,
		_varRhoFile, _varAcc, _varFSRRes);
  }
  else {
    addToVector(eeCovIdx,6,  _varBkg, _varDetRes, _varEff,
		_varRhoFile, _varAcc, _varFSRRes);
    addToVector(mmCovIdx,6,  _varBkg, _varDetRes, _varEff,
		_varRhoFile, _varAcc, _varFSRRes);
  }

  if (changeNames) {
    eeOldFName.push_back("cov_ee_varYield_2000.root");
    eeNewFName.push_back("cov_ee_varYield_5000.root");
    eeOldFName.push_back("cov_ee_varRhoFile_2000.root");
    eeNewFName.push_back("cov_ee_varRhoFile_5000.root");
    mmOldFName.push_back("cov_mumu_varYield_2000.root");
    mmNewFName.push_back("dymm_cov_varYield5000_20170202-covOnly.root");
    mmOldFName.push_back("cov_mumu_varRhoFile_2000.root");
    mmNewFName.push_back("cov_mumu_varRhoFile_5000.root");
  }

  TString fileTag="_corr";
  TString plotTag=" (corr)";

  finfomap_t eeCovFNames,mmCovFNames; // covariance derived by randomization
  finfomap_t eeUncTarget_stat, eeUncTarget_syst; // uncertainties from other
  finfomap_t mmUncTarget_stat, mmUncTarget_syst; // ... groups

  for (unsigned int i=0; i<eeCovIdx.size(); i++) {
    TVaried_t var= TVaried_t(eeCovIdx[i]);
    TString fname = "cov_ee_" + variedVarName(var) + TString("_2000.root");
    if (eeOldFName.size()>0) {
      for (unsigned int i=0; i<eeOldFName.size(); i++) {
	if (fname==eeOldFName[i]) fname=eeNewFName[i];
      }
    }
    eeCovFNames[var] = fname;
  }
  for (unsigned int i=0; i<mmCovIdx.size(); i++) {
    TVaried_t var= TVaried_t(mmCovIdx[i]);
    TString fname = "cov_mumu_" + variedVarName(var) + TString("_2000.root");
    if (mmOldFName.size()>0) {
      for (unsigned int i=0; i<mmOldFName.size(); i++) {
	if (fname==mmOldFName[i]) fname=mmNewFName[i];
      }
    }
    mmCovFNames[var] = fname;
  }

  std::cout << "ee covariance names:\n";
  for (finfomap_t::const_iterator it= eeCovFNames.begin();
       it!=eeCovFNames.end(); it++) {
    TString fname= it->second;
    std::cout << "  -- " << variedVarName(it->first) << "\t\t"
	      << it->second << "\n";
  }
  std::cout << "mm covariance names:\n";
  for (finfomap_t::const_iterator it= mmCovFNames.begin();
       it!=mmCovFNames.end(); it++) {
    TString fname= it->second;
    std::cout << "  -- " << variedVarName(it->first) << "\t\t"
	      << it->second << "\n";
  }

  std::vector<TH2D*> eeCovH2D, mmCovH2D;
  std::vector<TMatrixD> eeCovV, mmCovV;
  TMatrixD eeCovStat(43,43);
  TMatrixD eeCovSyst(eeCovStat);
  TMatrixD mmCovStat(eeCovStat);
  TMatrixD mmCovSyst(eeCovStat);

  // xxCovStat is the covariance of measured yield
  // xxCovSyst equals to other covariances
  // xxCovV contains all covariance matrices
  if (!loadCovData("ee",eeCovFNames,eeCovH2D,eeCovV,eeCovStat,eeCovSyst) ||
      !loadCovData("mm",mmCovFNames,mmCovH2D,mmCovV,mmCovStat,mmCovSyst)) {
    std::cout << "failed to load covariances\n";
    return;
  }

  // Adjust components of the uncertainties
  CovStruct_t eeCovS(eeCovStat); // creates and nullifies
  CovStruct_t mmCovS(mmCovStat);

  int plotCmpUnc= (hasValue("plotAdjUnc",showCanvs)) ? 1 : 0;
  if (hasValue("plotAdjCov",showCanvs)) plotCmpUnc=2;

  adjustMMUnc(mmCovFNames,mmCovV,mmCovS,plotChangedCov,plotCmpUnc);
  adjustEEUnc(eeCovFNames,eeCovV,eeCovS,plotChangedCov,plotCmpUnc);

  if (excludeSyst) {
    if (excludeSyst==1) {
      mmCovS.editSystNonAcc().Zero();
      eeCovS.editSystNonAcc().Zero();
    }
    mmCovS.editSystAcc().Zero();
    eeCovS.editSystAcc().Zero();
  }

  TH1D* h1csEE_tmp=loadHisto(eeCSFName, eeCSH1Name, "h1csEE_tmp", 1,h1dummy);
  TH1D* h1csMM_tmp=loadHisto(mmCSFName, mmCSH1Name, "h1csMM_tmp", 1,h1dummy);
  if (!h1csEE_tmp || !h1csMM_tmp) return;

  // plot input covariance in ee and mm channels
  if (hasValue("plotChannelInputCov",showCanvs)) {
    double ymin=0, ymax=0;
    ymin=1e-9; ymax=20.;
    mmCovS.Plot("mmCov",h1csMM_tmp);
    mmCovS.PlotUnc("mmCovUnc_",h1csMM_tmp,NULL,ymin,ymax);
    eeCovS.Plot("eeCov",h1csEE_tmp);
    eeCovS.PlotUnc("eeCovUnc_",h1csEE_tmp,NULL,ymin,ymax);
    //return;
  }

  // check whether mmCovS.Syst_Acc = eeCovS.Syst_Acc
  if (hasValue("plotChannelRelSystAccUnc",showCanvs)) {
    TH1D *h1mmAcc_rel= uncFromCov(mmCovS.covSyst_Acc,"h1mmAcc_rel",
				  h1csMM_tmp,h1csMM_tmp,0);
    TH1D *h1eeAcc_rel= uncFromCov(eeCovS.covSyst_Acc,"h1eeAcc_rel",
				  h1csEE_tmp,h1csEE_tmp,0);
    hsRed.SetStyle(h1eeAcc_rel);
    plotHisto(h1mmAcc_rel,"cAccUnc",1,1,"LP","#mu#mu");
    plotHistoSame(h1eeAcc_rel,"cAccUnc","LP","ee");
    //std::cout << "stopping\n"; return;
  }


  // create Acc theory-correlated uncertainty
  TMatrixD emCovTot(TMatrixD::kZero, eeCovS.covSyst_Acc);
  if (!excludeSyst) {
    double accCorrFactor=1.;
    emCovTot= constructEMAccCov(eeCovS,mmCovS,accCorrFactor);
  }

  if (corrCase==0) {
    double reductionFactor=1.;
    eeCovS.ReduceCorrelations(reductionFactor);
    mmCovS.ReduceCorrelations(reductionFactor);
    if (includeLumiUnc!=2) emCovTot.Zero();
  }

  if (0) {
    if (!eeCovS.PrintZRangeUnc("ee",h1csEE_tmp,0) ||
	!mmCovS.PrintZRangeUnc("mm",h1csMM_tmp,0)) {
      std::cout << "error\n";
    }
    return;
  }

  TMatrixD eeCovTot(eeCovS.Sum());
  TMatrixD mmCovTot(mmCovS.Sum());

  double lumiUnc=0.026;
  TString fileTagExtra,plotTagExtra;
  if (includeLumiUnc) {
    fileTagExtra="_wLumi"; plotTagExtra=" (wLumi)";
    if (includeLumiUnc==2) {
      eeCovTot.Zero(); mmCovTot.Zero(); emCovTot.Zero();
      lumiUnc=0.26;
      fileTagExtra="_wLumi10x"; plotTagExtra=" (wLumi10x)";
    }
    if (!addLumiCov(eeCovTot,mmCovTot,emCovTot, lumiUnc, h1csEE_tmp,h1csMM_tmp))
      return;
  }

  if (corrCase==0) {
    fileTag="_uncorr";
    plotTag=" (uncorr)";
  }
  else {
    fileTag="_corr";
    plotTag=" (corr)";
  }
  if (fileTagExtra.Length()) fileTag+=fileTagExtra;
  if (plotTagExtra.Length()) plotTag+=plotTagExtra;
  if (plotChangedCov) fileTag.Append("_plotChCov");

  if (!changeNames) fileTag.Append("_noChangeNames");
  if (useUncorrYieldUnc) fileTag.Append("_uncorrYieldUnc");
  if (excludeSyst) fileTag.Append(Form("_excludeSyst%d",excludeSyst));
  if (saveTag.Length()) fileTag.Append(saveTag);


  BLUEResult_t *blue=
    combineData(&eeCovTot,&mmCovTot,&emCovTot,fileTag,plotTag,printCanvases,
		showCombinationCanvases);
  if (!blue) { std::cout << "failed to get BLUE result\n"; return; }

  if (includeLumiUnc) {
    TMatrixD totFinalCov(*blue->getCov());
    TH1D *h1csLL_test= convert2histo1D(*blue->getEst(),totFinalCov,
				       "h1csLL_test","h1csLL test",NULL);
    if (!addLumiCov(totFinalCov,-lumiUnc,h1csLL_test)) return;
    plotCovCorr(totFinalCov,NULL,"h2finCovNoLumi","cFinCovNoLumi");
  }

  if (hasValue("plotInputContributedCovs",showCanvs)) {
    PlotCovCorrOpt_t optCC(1,1,1,1.8,0.15,0.15);
    TH1D* h1csEE_tmp=loadHisto(eeCSFName, eeCSH1Name, "h1csEE_tmp", 1,h1dummy);
    TH1D* h1csMM_tmp=loadHisto(mmCSFName, mmCSH1Name, "h1csMM_tmp", 1,h1dummy);
    if (1)
    for (int i=1; i<=4; i++) {
      TString tagName="ee_" + eeCovS.GetPartName(i);
      TMatrixD partM= eeCovS.GetPart(i);
      //plotCovCorr(partM,h1csEE_tmp,"h2_"+tagName,"c2_"+tagName,optCC);
      TMatrixD corr = blue->contributedCorrPartial(partM,1);
      TString cName= "cContrCov_" + tagName;
      TH2D *h2= convert2histo(corr,h1csEE_tmp,"h2"+tagName,"h2"+tagName);
      logAxis(h2);
      plotHisto(h2,cName,optCC);
    }
    for (int i=1; i<=4; i++) {
      TString tagName="mm_" + mmCovS.GetPartName(i);
      TMatrixD partM= mmCovS.GetPart(i);
      //plotCovCorr(partM,h1csMM_tmp,"h2_"+tagName,"c2_"+tagName,optCC);
      TMatrixD corr = blue->contributedCorrPartial(partM,0);
      TString cName= "cContrCov_" + tagName;
      TH2D *h2= convert2histo(corr,h1csMM_tmp,"h2"+tagName,"h2"+tagName);
      logAxis(h2);
      plotHisto(h2,cName,optCC);
    }
  }

  TH1D* h1csLL_tmp= cloneHisto(h1csEE_tmp, "h1csLL_tmp", "h1csLL_tmp");
  h1csLL_tmp->GetXaxis()->SetTitle( niceMassAxisLabel(2,"",0) );
  h1csLL_tmp->GetYaxis()->SetTitle( niceMassAxisLabel(2,"",1) );

  TH1D *h1csLL_badBinning= convert2histo1D(*blue->getEst(),*blue->getCov(),
			   "h1csLL_badBinning","h1csLL_badBinning",NULL);
  if (!h1csLL_badBinning) return;
  if (!copyContents(h1csLL_tmp, h1csLL_badBinning)) return;

  TH1D* h1deltacsLL_tmp= cloneHisto(h1csLL_tmp,
				    "h1deltacsLL_tmp", "h1deltacsLL_tmp");
  h1deltacsLL_tmp->GetXaxis()->SetTitle( niceMassAxisLabel(2,"",0) );
  h1deltacsLL_tmp->GetYaxis()->SetTitle( "\\delta" + niceMassAxisLabel(2,"",1) );

  PlotCovCorrOpt_t optCCNoLog;
  optCCNoLog.setLogScale(0,0);

  if (0) { // check function constructMeasCov
    TMatrixD covM_ee_yield= blue->measCovMatrix(eeCovS.getStatYield(),1);
    TMatrixD covM_mm_yield= blue->measCovMatrix(mmCovS.getStatYield(),0);
    //covM_ee_yield.Print();
    //covM_mm_yield.Print();
    TMatrixD covM_yield(covM_ee_yield,TMatrixD::kPlus,covM_mm_yield);
    plotCovCorr(covM_yield,NULL,"h2covMeasYield_chk","cCovMeasYield_chk");
    TMatrixD covFin_yield= blue->contributedCov(covM_yield);
    TH2D* h2cov_fromYield=NULL;
    plotCovCorr(covFin_yield,h1csLL_tmp,
		"h2covFin_fromYield_chk","cCovFin_fromYield_chk",
		PlotCovCorrOpt_t(),&h2cov_fromYield);
    TH1D *h1_dCS_fromYield= uncFromCov(covFin_yield,
				       "h1_dCS_fromYield_chk",
				       h1deltacsLL_tmp,NULL,0);
    printHisto(h2cov_fromYield);
    printHisto(h1_dCS_fromYield);
  }

  TMatrixD sumContributedCov(TMatrixD::kZero,eeCovTot);
  CovStruct_t llCovS(eeCovStat); // creates and nullifies
  for (int iFlag=1; iFlag<=4; iFlag++) {
    CovStruct_t::TCovPart_t part=CovStruct_t::TCovPart_t(iFlag);
    double accCorrFlag= (corrCase==1) ? 1. : 0.;
    TMatrixD measCov= constructMeasCov(eeCovS,mmCovS,part,accCorrFlag,blue);
    TString label= eeCovS.GetPartName(iFlag);
    if (0) plotCovCorr(measCov,NULL,"h2covMeas"+label,"cCovMeas_"+label,
		       optCCNoLog,NULL);
    TMatrixD covFin_contrib= blue->contributedCov(measCov);
    llCovS.SetPart(part,covFin_contrib);
    sumContributedCov+=covFin_contrib;
    if (hasValue("plotFinalCovByType",showCanvs)) {
      TH2D* h2cov_contrib=NULL;
      plotCovCorr(covFin_contrib,h1csLL_tmp,
		  "h2covFin_from"+label,"cCovFin_from_"+label,
		  PlotCovCorrOpt_t(),&h2cov_contrib);
      TH1D *h1_dCS_contrib= uncFromCov(covFin_contrib,
				       "h1_dCS_from_"+label,
				       h1deltacsLL_tmp,NULL,0);
      if (!h2cov_contrib || !h1_dCS_contrib) return;
      hsVec[iFlag-1].SetStyle(h1_dCS_contrib);
      plotHistoAuto(h1_dCS_contrib,"canvContrUnc",1,1,"LPE",label);
    }
  }

  if (1) {
    if (!eeCovS.PrintZRangeUnc("ee",h1csEE_tmp,0) ||
	!mmCovS.PrintZRangeUnc("mm",h1csMM_tmp,0) ||
	!llCovS.PrintZRangeUnc("ll",h1csLL_tmp,0)) {
      std::cout << "error\n";
    }

    if (0) {
      std::cout << "\n EE    MM    Combi-XSect\n";
      const TH1D *h1def=h1csEE_tmp;
      for (int ibin=1; ibin<=h1def->GetNbinsX(); ibin++) {
	std::cout << "ibin=" << ibin << " " << h1def->GetBinLowEdge(ibin)
		  << " -- "
		  << (h1def->GetBinLowEdge(ibin)+h1def->GetBinWidth(ibin))
		  << "  "
		  << h1csEE_tmp->GetBinContent(ibin) << "   "
		  << h1csMM_tmp->GetBinContent(ibin) << "   "
		  << h1csLL_tmp->GetBinContent(ibin) << "\n";
      }
    }
  }




  //TMatrixD sumContrChk( sumContributedCov, TMatrixD::kMinus, *blue->getCov() );
  //sumContrChk.Print();


  if (printCanvases || (saveTag.Length()>0)) {
    TString outFName="dyll-combi-" + fileTag + ".root";
    TFile fout(outFName,"recreate");
    blue->write(fout,fileTag);
    SaveCanvases("ALL","dyll-combi-"+fileTag,&fout);
    //h2cov_fromYield->Write();
    //h1_dCS_fromYield->Write();
    writeTimeTag(&fout);
    fout.Close();
    std::cout << "file <" << fout.GetName() << "> created\n";
  }
}

// -------------------------------------------------------

int adjustMMUnc(const finfomap_t &mmCovFNames,
		std::vector<TMatrixD> &mmCovV,
		CovStruct_t &mmCovS, int plotChangedCov,
		int plotCmp)
{
  finfomap_t mmUncFiles, mmUncStat, mmUncSyst;
  fvarweightmap_t relUncW;

  TH1D* h1csEE_loc=loadHisto(eeCSFName, eeCSH1Name, "h1csEE_loc", 1,h1dummy);
  TH1D* h1csMM_loc=loadHisto(mmCSFName, mmCSH1Name, "h1csMM_loc", 1,h1dummy);
  if (!h1csEE_loc || !h1csMM_loc) return 0;

  TString inpPath="./";
  TString inpFile=inpPath + "DYmm_ROOTFile_Input-20170207.root";

  TH1D *h1csMM= cloneHisto(h1csMM_loc,"h1csMM","h1csMM");
  if (0) {
    h1csMM= loadHisto(inpFile, "h_DiffXSec", "h1csMM_KL", 1,h1dummy);
    plotHisto(h1csMM_loc,"cCSmm",1,1,"LPE","old #mu#mu cs");
    plotHistoSame(h1csMM,"cCSmm","LPE","new #mu#mu cs");
    printRatio(h1csMM,h1csMM_loc);
    return 0;
  }

  mmUncFiles[_varYield]= inpFile;
  mmUncStat [_varYield]= "h_RelStatUnc";
  mmUncSyst [_varYield]= "";
  relUncW   [_varYield]= (useUncorrYieldUnc) ? 1. : 0.;
  mmUncFiles[_varBkg]= inpFile;
  mmUncStat [_varBkg]= "h_RelUnc_Stat_BkgEst";
  mmUncSyst [_varBkg]= "h_RelUnc_Syst_BkgEst";
  relUncW   [_varBkg]= 1;
  mmUncFiles[_varDetRes]= inpFile;
  mmUncStat [_varDetRes]= "h_RelUnc_Stat_DetRes";
  mmUncSyst [_varDetRes]= "h_RelUnc_Syst_DetRes";
  relUncW   [_varDetRes]= 1;
  mmUncFiles[_varEff]= inpFile;
  mmUncStat [_varEff]= "";
  mmUncSyst [_varEff]= "";
  relUncW   [_varEff]= 0;
  mmUncFiles[_varRhoFile]= inpFile;
  mmUncStat [_varRhoFile]= "h_RelUnc_Stat_EffSF";
  mmUncSyst [_varRhoFile]= "h_RelUnc_Syst_EffSF";
  relUncW   [_varRhoFile]= 1;
  mmUncFiles[_varFSRRes]= inpFile;
  mmUncStat [_varFSRRes]= "h_RelUnc_Stat_FSR";
  mmUncSyst [_varFSRRes]= "h_RelUnc_Syst_FSR";
  relUncW   [_varFSRRes]= 1;
  mmUncFiles[_varAcc]= inpFile;
  mmUncStat [_varAcc]= "";
  mmUncSyst [_varAcc]= "h_RelUnc_Syst_Acc";
  relUncW   [_varAcc]= 1;

  const TH1D *h1central= (1) ? h1csMM : NULL;
  int change_covV=1+(plotChangedCov!=0) ? 1:0;
  int res=adjustUnc( 0,
		     mmCovFNames, mmCovV,
		     mmUncFiles,mmUncStat,mmUncSyst,relUncW,
		     h1csMM,h1central,
		     plotCmp,change_covV,
		     "_mm",
		     mmCovS);
  if (!res) {
    std::cout << "error in adjustMMUnc\n";
    return 0;
  }

  return 1;
}
// -----------------------------------------------------------------------

int adjustEEUnc(const finfomap_t &eeCovFNames,
		std::vector<TMatrixD> &eeCovV,
		CovStruct_t &eeCovS, int plotChangedCov,
		int plotCmp)
{
  finfomap_t eeUncFiles, eeUncStat, eeUncSyst;
  fvarweightmap_t relUncW;

  TH1D* h1csEE_loc=loadHisto(eeCSFName, eeCSH1Name, "h1csEE_loc", 1,h1dummy);
  TH1D* h1csMM_loc=loadHisto(mmCSFName, mmCSH1Name, "h1csMM_loc", 1,h1dummy);
  if (!h1csEE_loc || !h1csMM_loc) return 0;

  TString inpPath="./";
  TString inpFile=inpPath + "DYee_ROOTFile_Input-20170208.root";

  TH1D *h1csEE= cloneHisto(h1csEE_loc,"h1csEE","h1csEE");
  if (0) {
    h1csEE= loadHisto(inpFile, "h_DiffXSec", "h1csEE_RC", 1,h1dummy);
    plotHisto(h1csEE_loc,"cCSee",1,1,"LPE","old #it{ee} cs");
    plotHistoSame(h1csEE,"cCSee","LPE","new #it{ee} cs");
    printRatio(h1csEE,h1csEE_loc);
    return 0;
  }

  eeUncFiles[_varYield]= inpFile;
  eeUncStat [_varYield]= "h_RelStatUnc";
  eeUncSyst [_varYield]= "";
  relUncW   [_varYield]= (useUncorrYieldUnc) ? 0.01 : 0.;
  eeUncFiles[_varBkg]= inpFile;
  eeUncStat [_varBkg]= "h_RelUnc_Stat_BkgEst";
  eeUncSyst [_varBkg]= "h_RelUnc_Syst_BkgEst";
  relUncW   [_varBkg]= 0.01;
  eeUncFiles[_varDetRes]= inpFile;
  eeUncStat [_varDetRes]= "h_RelUnc_Stat_DetRes"; // this is actually systematic
  eeUncSyst [_varDetRes]= "h_RelUnc_Syst_DetRes";
  relUncW   [_varDetRes]= 0.01;
  eeUncFiles[_varEff]= inpFile;
  eeUncStat [_varEff]= "";
  eeUncSyst [_varEff]= "";
  relUncW   [_varEff]= 0;
  eeUncFiles[_varRhoFile]= inpFile;
  eeUncStat [_varRhoFile]= "h_RelUnc_Stat_EffSF";
  eeUncSyst [_varRhoFile]= "h_RelUnc_Syst_EffSF";
  relUncW   [_varRhoFile]= 0.01;
  eeUncFiles[_varFSRRes]= inpFile;
  eeUncStat [_varFSRRes]= "";
  eeUncSyst [_varFSRRes]= "h_RelUnc_Syst_FSR";
  relUncW   [_varFSRRes]= 0.01;
  eeUncFiles[_varAcc]= inpFile;
  eeUncStat [_varAcc]= "";
  eeUncSyst [_varAcc]= "h_RelUnc_Syst_Acc";
  relUncW   [_varAcc]= 0.01;

  const TH1D *h1central= (1) ? h1csEE : NULL;
  int change_covV=1 + (plotChangedCov!=0) ? 1:0;
  int res=adjustUnc( 1,
		     eeCovFNames, eeCovV,
		     eeUncFiles,eeUncStat,eeUncSyst,relUncW,
		     h1csEE,h1central,
		     plotCmp,change_covV,
		     "_ee",
		     eeCovS);

  if (!res) {
    std::cout << "error in adjustEEUnc\n";
    return 0;
  }
  return 1;
}
// -----------------------------------------------------------------------

int loadUncData(int isEE,
		const finfomap_t &fNames,
		const finfomap_t &h1names,
		const fvarweightmap_t &relUncW,
		const TH1D *h1binning,
		const TH1D *h1centralVal, // if the uncertainty is relative
		std::vector<TH1D*> &h1uncV,
		TString tag)
{
  TString lepTag=(isEE) ? "ee" : "mm";

  std::cout << "\n";
  if (h1centralVal) printHisto(h1centralVal);


  for (finfomap_t::const_iterator it= fNames.begin(); it!=fNames.end(); it++) {
    TString fname= it->second;
    TString h1newName= "h1unc_" + lepTag + "_" + variedVarName(it->first) + tag;
    TH1D *h1unc= cloneHisto(h1binning, h1newName,h1newName);
    h1unc->Reset();
    h1uncV.push_back(h1unc);
    if (fname.Length()==0) continue;

    finfomap_t::const_iterator itHName= h1names.find(it->first);
    if (itHName==h1names.end()) {
      std::cout << "could not find the histo name for "
		<< variedVarName(it->first) << "\n";
      return 0;
    }
    if (itHName->second.Length()==0) continue;

    fvarweightmap_t::const_iterator itRelW= relUncW.find(it->first);
    if (itRelW==relUncW.end()) {
      std::cout << "coult not find the relUncW for "
		<< variedVarName(it->first) << "\n";
      return 0;
    }

    TH1D *h1= loadHisto(fname,itHName->second,itHName->second + "_tmp",h1dummy);
    if (!h1) return 0;
    std::cout << "got " << h1newName << "\n";
    //printHisto(h1);
    //std::cout << "\n\n" << std::endl;
    if (itRelW->second!=0.) {
      if (!copyContents(h1unc, h1)) return 0;
    }
    if ((h1centralVal!=NULL) && (itRelW->second!=0.)) {
      std::cout << "scaling uncertainty for " << variedVarName(it->first) << "\n";
      h1unc->Multiply(h1centralVal);
      h1unc->Scale(itRelW->second);
    }
  }
  return 1;
}

// -------------------------------------------------------

int aggregateUnc(const finfomap_t &fnames,
		 const std::vector<TMatrixD> &covStatV,
		 const std::vector<TH1D*> &h1SystV,
		 CovStruct_t &covS)
{

  covS.Zero();

  int idx= findVaried(fnames,_varYield);
  if (idx==-1) {
    std::cout << "failed to find varYield\n";
    return 0;
  }
  covS.covStat_Yield= covStatV[idx];

  idx=0;
  for (finfomap_t::const_iterator it=fnames.begin(); it!=fnames.end();
       it++, idx++) {
    if (it->first == _varYield) continue;
    covS.covStat_nonYield += covStatV[idx];
  }

  if (fnames.size() != h1SystV.size()) {
    std::cout << "aggregateUnc: fnames.size=" << fnames.size() << ", "
	      << h1SystV.size() << "\n";
    return 0;
  }

  idx=0;
  for (finfomap_t::const_iterator it=fnames.begin(); it!=fnames.end();
       it++, idx++) {
    if (h1SystV[idx]->Integral()==0) {
      std::cout << "aggregate systematic error: skipping "
		<< variedVarName(it->first) << " (no input)\n";
      continue;
    }
    if (it->first == _varAcc) {
      covS.covSyst_Acc = convert2cov1D(h1SystV[idx],0);
    }
    else {
      TMatrixD tmpM= convert2cov1D(h1SystV[idx],0);
      std::cout << "h1SystV[idx].size=" << h1SystV[idx]->GetNbinsX()
		<< ", tmpM.size=" << tmpM.GetNrows() << " x " << tmpM.GetNcols()
		<< ", covS.covSyst_nonAcc.size="
		<< covS.covSyst_nonAcc.GetNrows() << " x "
		<< covS.covSyst_nonAcc.GetNcols() << "\n";
      covS.covSyst_nonAcc += convert2cov1D(h1SystV[idx],0);
    }
  }

  return 1;
}

// -------------------------------------------------------

int adjustUnc(int isEE,
	      const finfomap_t &covFNames,
	      std::vector<TMatrixD> &covV,
	      const finfomap_t &uncFiles,
	      const finfomap_t &uncStat,
	      const finfomap_t &uncSyst,
	      const fvarweightmap_t &relUncW,
	      const TH1D *h1binning,
	      const TH1D *h1central,
	      int plotCmp, // 1 - plot uncertainties, 2 - and covariances
	      int change_covV, // 1 - change values, 2 - verify change by plot
	      TString tag,
	      CovStruct_t &covS)
{
  if (uncFiles.size()!=covV.size()) {
    std::cout << "sizes are different!\n";
  }

  std::vector<TString> labelV;
  for (finfomap_t::const_iterator it= covFNames.begin();
       it!=covFNames.end(); it++) {
    TVaried_t var= it->first;
    if (uncFiles.find(var)==uncFiles.end()) {
      std::cout << "mmUncFiles does not contain " << variedVarName(var) << "\n";
      return 0;
    }
    labelV.push_back(variedVarName(var) + tag);
  }

  std::vector<TH1D*> h1uncStatV, h1uncSystV;
  if (!loadUncData(isEE,uncFiles,uncStat,relUncW,
		   h1binning,h1central,h1uncStatV,tag+"_stat")) {
    std::cout << "failed to load statistical component\n";
    return 0;
  }
  if (!loadUncData(isEE,uncFiles,uncSyst,relUncW,
		   h1binning,h1central,h1uncSystV,tag+"_syst")) {
    std::cout << "failed to load systematic component\n";
    return 0;
  }

  // correct for the input: Stat_DetRes in the electron channel is a part of
  // Syst_DetRes
  if (isEE) {
    int idx= findVaried(covFNames,_varDetRes);
    if (idx==-1) {
      std::cout << "failed to find detRes (needed for the electron channel)\n";
      return 0;
    }
    std::cout << "idx(_varDetRes)=" << idx << "\n";
    //printHisto(h1uncStatV[idx]);
    //printHisto(h1uncSystV[idx]);
    if (!addInQuadrature(h1uncSystV[idx],h1uncStatV[idx],1)) {
      return 0;
    }
    h1uncStatV[idx]->Reset();
    std::cout << "cleared h1detRes_ee statistical error from RC\n";
    //return 0;
  }

  // eliminate eff syst. uncertainty
  if (1 && change_covV) {
    int idx1= findVaried(covFNames,_varEff);
    int idx2= findVaried(covFNames,_varAcc);
    if ((idx1==-1) || (idx2==-1)) {
      std::cout << "failed to find eff or acc\n";
      return 0;
    }
    std::cout << "idx(_varEff)=" << idx1 << "\n";
    std::cout << "idx(_varAcc)=" << idx2 << "\n";
    if (change_covV) {
      covV[idx1].Zero();
      covV[idx2].Zero();
    }
    /*
    if (idx==0) { std::cout << "  --- idx=0, fix the code\n"; return 0; }
    if (!h1uncStatV[idx]) {
      h1uncStatV[idx]= cloneHisto(h1uncStatV[0],Form("h1eff_%d",isEE),"h1Eff");
      h1uncStatV[idx]->Reset();
    }
    */
  }

  int res=compareUnc(labelV,h1uncStatV,covV,plotCmp,change_covV,tag);
  if (res && plotCmp && (change_covV==2))
    res=compareUnc(labelV,h1uncStatV,covV,2,0,tag+"_adj");

  if (res) {
    res= aggregateUnc(covFNames,covV,h1uncSystV,covS);
  }

  if (!res) {
    std::cout << "error in adjustUnc\n";
    return 0;
  }

  return 1;
}

// -------------------------------------------------------

int compareUnc(const std::vector<TString> &labelV,
	       const std::vector<TH1D*> &h1V_target,
	       std::vector<TMatrixD> &covV,
	       int plotCmp,
	       int change_covV,
	       TString tag)
{
  for (unsigned int i=0; i<labelV.size(); i++) {
    TString h2Name=h1V_target[i]->GetName() + TString("_cov");
    TH2D *h2cov= convert2histo(covV[i],h1V_target[i],h2Name,h2Name,NULL);
    TH1D *h1uncFromCov= uncFromCov(h2cov,NULL,0);
    if (plotCmp) {
      TString canvName="cCmpUnc_" + labelV[i] + "_" + tag;
      TH1D *h1= cloneHisto(h1V_target[i],
			   h1V_target[i]->GetName() + TString("_plotClone")+tag,
			   h1V_target[i]->GetTitle());
      hsRed.SetStyle(h1);
      h1->GetYaxis()->SetTitleOffset(1.5);
      logAxis(h1,1+8,"M [GeV]","","");
      TCanvas *cx=plotHisto(h1,canvName,1,1,"LP",labelV[i] + " target");
      setLeftMargin(cx,0.15);
      moveLegend(cx,0.05,0.);
      plotHistoSame(h1uncFromCov,canvName,"LP",labelV[i] + " cov");
    }
    if (plotCmp==2) {
      TString canvName2= "cCov_" + labelV[i] + "_" + tag;
      PlotCovCorrOpt_t opt;
      plotCovCorr(covV[i],h1V_target[1],"h2cov_"+labelV[i]+"_"+tag,
		  canvName2,opt);
    }

    if (change_covV) {
      if (h1V_target[i]->Integral()!=double(0)) {
	std::cout << "adjusting covV of " << labelV[i] << "\n";
	if (!changeCov(covV[i],h1V_target[i])) return 0;
      }
      else {
	std::cout << "no adjustment of covV for " << labelV[i] << "\n";
      }
    }
  }
  return 1;
}

// -------------------------------------------------------
// -------------------------------------------------------

TMatrixD constructEMAccCov(const CovStruct_t &eeCovS, const CovStruct_t &mmCovS,
			   double corrFactor)
{
  TMatrixD emCovTot(TMatrixD::kZero, eeCovS.covSyst_Acc);
  for (int ir=0; ir<eeCovS.covSyst_Acc.GetNrows(); ir++) {
    for (int ic=0; ic<mmCovS.covSyst_Acc.GetNcols(); ic++) {
      emCovTot(ir,ic) = corrFactor *
	sqrt(eeCovS.covSyst_Acc(ir,ic)) *
	sqrt(mmCovS.covSyst_Acc(ir,ic));
    }
  }
  return emCovTot;
}

// -------------------------------------------------------

TMatrixD constructMeasCov(const CovStruct_t &eeCovS, const CovStruct_t &mmCovS,
			  CovStruct_t::TCovPart_t covFlag, double accCorrFactor,
			  BLUEResult_t *blue)
{
  TMatrixD measCov_ee(blue->measCovMatrix(eeCovS.GetPart(covFlag),1));
  TMatrixD measCov_mm(blue->measCovMatrix(mmCovS.GetPart(covFlag),0));
  TMatrixD measCov( measCov_ee, TMatrixD::kPlus, measCov_mm );
  if (covFlag==CovStruct_t::_varCov_systAcc) {
    TMatrixD emCov( constructEMAccCov(eeCovS,mmCovS,accCorrFactor) );
    TMatrixD measEMCov( blue->measCovMatrix(emCov,2) );
    measCov += measEMCov;
  }
  return measCov;
}

// -------------------------------------------------------
