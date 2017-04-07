#ifndef CovStruct_H
#define CovStruct_H

//#include "inputs.h"
//#include "Blue.h"
#include "analyseBLUEResult.h"
#include "DYbinning.h"

// -------------------------------------------------------
// -------------------------------------------------------

struct CovStruct_t {
  typedef enum { _varCov_none=0, _varCov_statYield, _varCov_statNonYield,
		 _varCov_systAcc, _varCov_systNonAcc } TCovPart_t;

  int plotStatYieldOnly;
  TMatrixD covStat_Yield, covStat_nonYield;
  TMatrixD covSyst_Acc, covSyst_nonAcc;
public:
  CovStruct_t(const TMatrixD &M) :
    plotStatYieldOnly(0),
    covStat_Yield(TMatrixD::kZero,M), covStat_nonYield(TMatrixD::kZero,M),
    covSyst_Acc(TMatrixD::kZero,M), covSyst_nonAcc(TMatrixD::kZero,M) {}

  CovStruct_t(const TMatrixD &set_covStatYield,
	      const TMatrixD &set_covStatNonYield,
	      const TMatrixD &set_covSystAcc,
	      const TMatrixD &set_covSystNonAcc) :
    plotStatYieldOnly(0),
    covStat_Yield(set_covStatYield),
    covStat_nonYield(set_covStatNonYield),
    covSyst_Acc(set_covSystAcc),
    covSyst_nonAcc(set_covSystNonAcc)
  {}

  CovStruct_t(const CovStruct_t &s) :
    plotStatYieldOnly(s.plotStatYieldOnly),
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

  void Scale(double x)
  {
    covStat_Yield*=x; covStat_nonYield*=x;
    covSyst_Acc*=x; covSyst_nonAcc*=x;
  }

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

  // index value will be for the array, not for histogram bin
  int PrintZRangeUnc(TString label,
		     int idxMin, int idxMax_p1,
		     const TH1D *h1cs_perMBW,
		     int printNumbers=0) const
  {
    std::cout << "ZmassRange: " << idxMin << " .. " << idxMax_p1 << "\n";
    std::vector<double> zCov;
    TH1D *h1cs_abs= timesMassBinWidth(h1cs_perMBW);
    if (!this->ZrangeCov(idxMin,idxMax_p1,h1cs_abs, zCov,printNumbers)) {
      std::cout << "CovStruct_t::PrintZRangeUnc: failed to get the ZrangeCov\n";
      return 0;
    }
    std::cout << "integrated uncertainties:\n";
    if (1 && printNumbers) {
      printHisto(h1cs_abs);
      std::cout << " adding numbers:\n";
      for (int ibin=idxMin+1; ibin<idxMax_p1; ibin++) {
	std::cout << "ibin=" << ibin << "  " << h1cs_abs->GetBinContent(ibin)
		  << "\n";
      }
    }
    std::cout << "  - " << label << ": "
	      << h1cs_abs->Integral(idxMin+1,idxMax_p1) << " ";
    for (int ii=0; ii<4; ii++) {
      std::cout << "+- " << sqrt(zCov[ii]) << " (" << this->GetPartName(ii+1) << ")";
    }
    std::cout << "\n";
    return 1;
  }

  // automatic determination
  int PrintZRangeUnc(TString label, const TH1D *h1cs_perMBW,
		     int printNumbers=0) const
  {
    int idxMin=-1, idxMax_p1=-1;
    // index value will be for the array, not for histogram bin
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
      std::cout << "+- " << sqrt(zCov[ii]) << " (" << this->GetPartName(ii+1) << ")";
    }
    std::cout << "\n";
    return 1;
  }

  void Plot(TString tag, TH1D *h1binning, const PlotCovCorrOpt_t &opt) const {
    plotCovCorr(covStat_Yield,h1binning,
		"h2CovStatYield_"+tag,"c2CovStatYield_"+tag,opt);
    if (!plotStatYieldOnly) {
      plotCovCorr(covStat_nonYield,h1binning,
		  "h2CovStatNonYield_"+tag,"c2CovStatNonYield_"+tag,opt);
      plotCovCorr(covSyst_Acc,h1binning,
		  "h2CovSystAcc_"+tag,"c2CovSystAcc_"+tag,opt);
      plotCovCorr(covSyst_nonAcc,h1binning,
		  "h2CovSystNonAcc_"+tag,"c2CovSystNonAcc_"+tag,opt);
    }
  }

  void PlotUnc(TString tag, TH1D *h1binning, TH1D *h1centralVal=NULL,
	       double yrangeMin=0, double yrangeMax=0,
	       int logX=1, int logY=1) {
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
    plotHisto(h1statYield,cName,logX,logY,"LP","stat Yield");
    if (!plotStatYieldOnly) {
      plotHistoSame(h1statNonYield,cName,"LP","stat nonYield");
      plotHistoSame(h1systAcc,cName,"LP","syst Acc");
      plotHistoSame(h1systNonAcc,cName,"LP","syst nonAcc");
    }
    plotHistoSame(h1tot,cName,"LP","total");
  }

};

// -------------------------------------------------------
// -------------------------------------------------------

TMatrixD constructEMAccCov(const CovStruct_t &eeCovS, const CovStruct_t &mmCovS,
			   double corrFactor=1.);
TMatrixD constructMeasCov(const CovStruct_t &eeCovS, const CovStruct_t &mmCovS,
			  CovStruct_t::TCovPart_t covFlag, double accCorrFactor,
			  const BLUEResult_t *blue);

CovStruct_t constructCombCov(const CovStruct_t &eeCovS,
			     const CovStruct_t &mmCovS,
			     const BLUEResult_t *blue,
			     int corrAccCase,
			     int &flag);

int aggregateUnc(const finfomap_t &fnames, // needed for first value
		 const std::vector<TMatrixD> &covStatV,
		 const std::vector<TH1D*> &h1SystV,
		 CovStruct_t &covM);

int adjustUnc(TString lepTag,
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
	      CovStruct_t &covM,
	      std::string specialArgs);

CovStruct_t* detailedCovPlots(const BLUEResult_t *blue,
			      int corrCase, // flag
			      const TH1D *h1csEE,
			      const TH1D *h1csMM,
			      const CovStruct_t &eeCovS,
			      const CovStruct_t &mmCovS,
			      std::string showCanvs,
			      PlotCovCorrOpt_t *optCC_user=NULL);

// -------------------------------------------------------
// -------------------------------------------------------

#endif
