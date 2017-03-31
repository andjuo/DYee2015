#ifndef analyseBLUEResult_H
#define analyseBLUEResult_H

#include "inputs.h"
#include "crossSection.h"
#include "Blue.h"

// -------------------------------------------------------
// -------------------------------------------------------

int loadCovData(TString lepTag,
		const finfomap_t &covFNames,
		std::vector<TH2D*> &covAsH2D, std::vector<TMatrixD> &covV,
		TMatrixD &eeCovStat, TMatrixD &eeCovSyst);
int addLumiCov(TMatrixD &cov, double lumiRelUnc, const TH1D *h1cs);
int addLumiCov(TMatrixD &eeCov, TMatrixD &mmCov, TMatrixD &emCov,
	       double lumiRelUnc, const TH1D *h1csEE, const TH1D *h1csMM);


BLUEResult_t* combineData(const TMatrixD *covEE_inp,
			  const TMatrixD *covMM_inp,
			  const TMatrixD *covEM_inp,
			  TString outputFileTag, TString plotTag,
			  int printCanvases, std::string showCanvases="ALL");

// -------------------------------------------------------
// -------------------------------------------------------

#endif
