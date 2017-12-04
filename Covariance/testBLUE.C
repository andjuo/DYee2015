#include "inputs.h"
#include "crossSection.h"
#include "Blue.h"
#include "analyseBLUEResult.h"

void testBLUE(int includeLumiUnc=1)
{

  TString fname="dyll-combi-_corr_wLumi_inpYieldUnc.root";
  TString subdir="BLUE_inp/";
  TString fileTag="_testBLUEresult";
  TString plotTag="_testBLUE";
  int printCanvases=0;


  std::string showCombinationCanvases;
  if (0) {
    showCombinationCanvases="ALL";
  }
  else {
    showCombinationCanvases+= " showCombiCS";
    showCombinationCanvases+= " showCombiCSCheck";
    showCombinationCanvases+= " showCombiCSUnc";
    showCombinationCanvases+= " showCSCheck";
  }


  TFile fin(fname);
  if (!fin.IsOpen()) {
    std::cout << "failed to open the file <" << fin.GetName() << ">\n";
    return;
  }

  TH1D *h1csEE_loc= loadHisto(fin,subdir+"h1csEE","h1csEE_loc",1,h1dummy);
  TH1D *h1csMM_loc= loadHisto(fin,subdir+"h1csMM","h1csMM_loc",1,h1dummy);
  TH1D *h1csLL_orig= loadHisto(fin,subdir+"h1csLL","h1csLL_orig",1,h1dummy);
  TMatrixD *covEE= (TMatrixD*)fin.Get(subdir+"covTotEE_noLumi");
  TMatrixD *covMM= (TMatrixD*)fin.Get(subdir+"covTotMM_noLumi");
  TMatrixD *covEM= (TMatrixD*)fin.Get(subdir+"covTotEM_noLumi");
  TMatrixD* covTotLumi= (TMatrixD*)fin.Get(subdir+"covTotLumi");
  fin.Close();

  if (!covEE || !covMM || !covEM || !covTotLumi) {
    std::cout << "cov matrix is not loaded\n";
    nullPtr(covEE,"covEE");
    nullPtr(covMM,"covMM");
    nullPtr(covEM,"covEM");
    nullPtr(covTotLumi,"covTotLumi");
    return;
  }

  if (0) {
    printMDim("covEE",*covEE);
    printMDim("covMM",*covMM);
    printMDim("covEM",*covEM);
    printMDim("covTotLumi",*covTotLumi);
    std::cout << "\n";
    return;
  }

  if (includeLumiUnc) {
    if (!covTotLumi) {
      std::cout << "includeLumiUnc=1, but covTotLumi is null\n";
      return;
    }

    int nBins=43;
    for (int ir=0; ir<nBins; ir++) {
      for (int ic=0; ic<nBins; ic++) {
	(*covEE)(ir,ic) += (*covTotLumi)(ir,ic);
      }
    }
    for (int ir=0; ir<nBins; ir++) {
      for (int ic=0; ic<nBins; ic++) {
	(*covMM)(ir,ic) += (*covTotLumi)(ir+nBins,ic+nBins);
      }
    }
    for (int ir=0; ir<nBins; ir++) {
      for (int ic=0; ic<nBins; ic++) {
	(*covEM)(ir,ic) += (*covTotLumi)(ir,ic+nBins);
      }
    }
  }

  eeCSFName=fname;
  eeCSHisto1DName= subdir+"h1csEE";
  mmCSFName=fname;
  mmCSHisto1DName= subdir+"h1csMM";
  theoryCSFName=fname;
  theoryCSHisto1D_perMBW=1;
  theoryCSHisto1DName= "h1csTheoryRed";
  theoryCSHisto1DName="h1Combi";
  theoryCSHisto1DName= subdir+"h1csLL";

  PlotCovCorrOpt_t ccOpt;
  double internalScale=1.;
  BLUEResult_t *blue=
    combineData(covEE,covMM,covEM,fileTag,plotTag,printCanvases,
		showCombinationCanvases,internalScale,&ccOpt);
		
  //printHisto(h1csLL_orig);


  if (showCombinationCanvases.find("showCSCheck")!=std::string::npos) {
    std::cout << "\ncompare CSs\n";
    if (0) {
      h1csEE_loc->SetMarkerStyle(kMultiply);
      h1csMM_loc->SetMarkerStyle(kMultiply);
      plotHistoSame(h1csEE_loc,"cCombiCS","LPE","ee_loc");
      plotHistoSame(h1csMM_loc,"cCombiCS","LPE","#mu#mu_loc");
    }
    else {
      h1csLL_orig->SetMarkerStyle(24);
      plotHistoSame(h1csLL_orig,"cCombiCS","LPE","\\ell\\ell\\ orig");
      TH1D *h1csLL_origErr= errorAsCentral(h1csLL_orig,0);
      plotHistoSame(h1csLL_origErr,"cErr","LP","\\ell\ell\\ orig");
    }
  }
}
