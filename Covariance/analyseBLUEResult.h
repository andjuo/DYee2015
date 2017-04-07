#ifndef analyseBLUEResult_H
#define analyseBLUEResult_H

#include "inputs.h"
#include "crossSection.h"
#include "Blue.h"

extern TString eeCSFName;
extern TString mmCSFName;
extern TString theoryCSFName;
extern TString eeCSHisto1DName;
extern TString mmCSHisto1DName;
extern TString theoryCSHisto1DName;
extern int theoryCSHisto1D_perMBW; // whether theory is per mass bin width

// -------------------------------------------------------
// -------------------------------------------------------

int loadUncData(TString lepTag,
		const finfomap_t &fNames,
		const finfomap_t &h1names,
		const fvarweightmap_t &relUncW,
		const TH1D *h1binning,
		const TH1D *h1centralVal, // if the uncertainty is relative
		std::vector<TH1D*> &h1uncV,
		TString tag);

int loadCovData(TString lepTag,
		const finfomap_t &covFNames,
		std::vector<TH2D*> &covAsH2D, std::vector<TMatrixD> &covV,
		TMatrixD &eeCovStat, TMatrixD &eeCovSyst);

int addLumiCov(TMatrixD &cov, double lumiRelUnc, const TH1D *h1cs);
int addLumiCov(TMatrixD &eeCov, TMatrixD &mmCov, TMatrixD &emCov,
	       double lumiRelUnc, const TH1D *h1csEE, const TH1D *h1csMM);

int compareUncAndScale(const std::vector<TString> &labelV,
		       const std::vector<TH1D*> &h1V_target,
		       std::vector<TMatrixD> &covV,
		       int plotCmp,
		       int change_covV,
		       TString tag);

BLUEResult_t* combineData(const TMatrixD *covEE_inp,
			  const TMatrixD *covMM_inp,
			  const TMatrixD *covEM_inp,
			  TString outputFileTag, TString plotTag,
			  int printCanvases, std::string showCanvases="ALL",
			  double internal_scale=1.,
			  PlotCovCorrOpt_t *ccOpt_user=NULL);

// -------------------------------------------------------
// -------------------------------------------------------

#endif
