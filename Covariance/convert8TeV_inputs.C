#include "inputs.h"
#include "helper_8TeV.h"

// ---------------------------------------------------------------

void convert8TeV_inputs(int is2D=1, int saveFile=0)
{
  std::cout << "is2D=" << is2D << "\n";
  TString path="materials_SMP14003/";
  TString tag="NONE";

  HistoStyle_t hsEE(hsGreen);
  HistoStyle_t hsMM(hsColor46);
  HistoStyle_t hsLL(hsBlue);
  hsMM.markerStyle= 21;

  switch(is2D) {
  case 0: tag="1D"; break;
  case 1: tag="2D_mdf"; break;
  case 2: tag="2D"; break;
  default: std::cout << "bad is2D value\n"; return;
  }

  path += tag + "/";

  TString fname_xsEE1D= "ee_results_2014_abs_adhoc_PI_Bayes.dat";
  TString fname_xsMM1D= "mumu_results_2014_abs_PI_Bayes.dat";
  TString fname_xsLL1D= "output_data.txt";
  TString fname_covEE1D= "covForBlue-1D-absolute_EE_Bayes.dat";
  TString fname_covMM1D= "covForBlue-1D-absolute_MuMu.dat";
  TString fname_covLL1D= "output_covMat.txt";

  TString fname_xsEE2D= "ee_results_2D_2014_abs_adhoc_PI_Bayes.dat";
  TString fname_xsMM2D= "mumu_results_2D_2014_abs_PI_Bayes.dat";
  TString fname_xsLL2D= "output_data.txt";
  TString fname_covEE2D= "covForBlue-2D-absolute_EE_Bayes.dat";
  TString fname_covMM2D= "covForBlue-2D-absolute_MuMu.dat";
  TString fname_covLL2D= "output_covMat.txt";

  TH1D *h1ee=NULL, *h1mm=NULL, *h1ll=NULL;
  TH2D *h2covEE=NULL, *h2covMM=NULL, *h2covLL=NULL;

  // This is 1D case
  if (is2D==0) {

    if (1) {
      h1ee= loadASfile_asHisto1D(path + fname_xsEE1D,"h1ee",1);
      h1mm= loadASfile_asHisto1D(path + fname_xsMM1D,"h1mm",1);
      h1ll= loadASfile_asHisto1D(path + fname_xsLL1D,"h1ll",2);

      hsEE.SetStyle(h1ee);
      hsMM.SetStyle(h1mm);
      hsLL.SetStyle(h1ll);

      plotHisto(h1ee,"cXS1D",1,1,"LPE1","ee");
      plotHistoSame(h1mm,"cXS1D","LPE1","#mu#mu");
      plotHistoSame(h1ll,"cXS1D","LPE1","#it{ll}");
    }

    TH2D *h2covEE_abs= loadBLUECov(path + fname_covEE1D,
				   "h2covEE_abs", "h2covEE_abs",
				   _nMassBins2012,_massBinLimits2012,0);
    h2covEE= perMassBinWidth(h2covEE_abs,0);
    h2covMM= loadBLUECov(path + fname_covMM1D,
			 "h2covMM","h2covMM",
			 _nMassBins2012,_massBinLimits2012,0);
    h2covLL= loadBLUECov(path + fname_covLL1D,
			 "h2covLL","h2covLL",
			 _nMassBins2012,_massBinLimits2012,2);

    PlotCovCorrOpt_t optCC;
    plotCovCorr(h2covEE,"cCovEE",optCC);
    plotCovCorr(h2covMM,"cCovMM",optCC);
    plotCovCorr(h2covLL,"cCovLL",optCC);

    TH1D *h1accRelUnc= loadHisto(path + "syst_acc.root","invm_FEWZ41",
				 "h1accUncRel",h1dummy);
    if (!h1accRelUnc) return;
    // for some reason Alexey reduced the uncertainty by sqrt(1.4)
    h1accRelUnc->Scale(1/sqrt(1.4));

    if (1) {
      TH1D *h1uncEE_cov= uncFromCov(h2covEE);
      TH1D *h1uncMM_cov= uncFromCov(h2covMM);
      TH1D *h1uncLL_cov= uncFromCov(h2covLL);
      TH1D *h1uncLL= (h1ll) ? errorAsCentral(h1ll,0) : NULL;
      TH1D *h1ll_wAccUnc=
	(h1ll) ? addUncertainty("h1ll_wAccUnc",h1ll,h1accRelUnc,h1ll) : NULL;
      TH1D *h1uncLL_wAccUnc=errorAsCentral(h1ll_wAccUnc,0);
      hsEE.SetStyle(h1uncEE_cov);
      hsMM.SetStyle(h1uncMM_cov);
      hsLL.SetStyle(h1uncLL_cov);
      if (h1uncLL) hsViolet.SetStyle(h1uncLL);
      plotHisto(h1uncEE_cov,"cXS1Dunc",1,1,"LP","ee unc(cov)");
      plotHistoSame(h1uncMM_cov,"cXS1Dunc","LP","#mu#mu unc(cov)");
      plotHistoSame(h1uncLL_cov,"cXS1Dunc","LP","#it{ll} unc(cov)");
      //if (h1uncLL) plotHistoSame(h1uncLL,"cXS1Dunc","LP","#it{ll} unc(BLUE)");
      if (0 && h1uncLL_wAccUnc) {
	hsColor6(h1uncLL_wAccUnc);
	plotHistoSame(h1uncLL_wAccUnc,"cXS1Dunc","LP",
		      "#it{ll} unc(BLUE)+accUnc");
      }
      if (0) {
	TH1D *h1xsec= loadHisto(path+"rshape_comb_41bin.root","hxsec",
				"h1xsec_orig",h1dummy);
	hsRed.SetStyle(h1xsec);
	TH1D *h1uncOrig= errorAsCentral(h1xsec);
	plotHistoSame(h1uncOrig,"cXS1Dunc","LP","#it{ll} (full orig)");
      }
    }
  }

  // This is 2D case
  if (is2D==1) {
    if (1) {
      h1ee= loadASfile_asOneHisto2D(path + fname_xsEE2D,"h1ee",0);
      h1mm= loadASfile_asOneHisto2D(path + fname_xsMM2D,"h1mm",0);
      h1ll= loadASfile_asOneHisto2D(path + fname_xsLL2D,"h1ll",2);

      hsEE.SetStyle(h1ee);
      hsMM.SetStyle(h1mm);
      hsLL.SetStyle(h1ll);

      plotHisto(h1ee,"cXS1D",0,0,"LPE1","ee");
      plotHistoSame(h1mm,"cXS1D","LPE1","#mu#mu");
      plotHistoSame(h1ll,"cXS1D","LPE1","#it{ll}");
    }

    TH2D *h2covEE_abs= loadBLUECov(path + fname_covEE2D,
				   "h2covEE_abs", "h2covEE_abs",
				   _nFlatBins2D,_flatBins2D,0);
    h2covEE= perMassRapidityBinWidth(h2covEE_abs,0);
    h2covMM= loadBLUECov(path + fname_covMM2D,
			 "h2covMM","h2covMM",
			 _nFlatBins2D,_flatBins2D,0);
    h2covLL= loadBLUECov(path + fname_covLL1D,
			 "h2covLL","h2covLL",
			 _nFlatBins2D,_flatBins2D,2);

    PlotCovCorrOpt_t optCC;
    optCC.setLogScale(0,0);
    plotCovCorr(h2covEE,"cCovEE",optCC);
    plotCovCorr(h2covMM,"cCovMM",optCC);
    plotCovCorr(h2covLL,"cCovLL",optCC);

    if (1) {
      TH1D *h1uncEE_cov= uncFromCov(h2covEE);
      TH1D *h1uncMM_cov= uncFromCov(h2covMM);
      TH1D *h1uncLL_cov= uncFromCov(h2covLL);
      TH1D *h1uncLL= (h1ll) ? errorAsCentral(h1ll,0) : NULL;
      hsEE.SetStyle(h1uncEE_cov);
      hsMM.SetStyle(h1uncMM_cov);
      hsLL.SetStyle(h1uncLL_cov);
      if (h1uncLL) hsViolet.SetStyle(h1uncLL);
      plotHisto(h1uncEE_cov,"cXS1Dunc",0,1,"LP","ee unc(cov)");
      plotHistoSame(h1uncMM_cov,"cXS1Dunc","LP","#mu#mu unc(cov)");
      plotHistoSame(h1uncLL_cov,"cXS1Dunc","LP","#it{ll} unc(cov)");
      if (h1uncLL) plotHistoSame(h1uncLL,"cXS1Dunc","LP","#it{ll} unc(BLUE)");
      if (0) {
	TH1D *h1xsec= loadHisto(path+"absex_DET2D_comb.root","hxsec",
				"h1xsec_orig",h1dummy);
	hsRed.SetStyle(h1xsec);
	TH1D *h1uncOrig= errorAsCentral(h1xsec);
	plotHistoSame(h1uncOrig,"cXS1Dunc","LP","#it{ll} (full orig)");
      }
    }
  }


  if (saveFile) {
    if (!h1ee || !h1mm || !h1ll || !h2covEE || !h2covMM || !h2covLL) {
      std::cout << "at least one needed pointer is null\n";
      return;
    }

    h1ee->SetStats(0);
    h1mm->SetStats(0);
    h1ll->SetStats(0);
    h2covEE->SetStats(0);
    h2covMM->SetStats(0);
    h2covLL->SetStats(0);

    TString foutName="out_dy8TeV_1D.root";
    foutName.ReplaceAll("_1D","_"+tag);
    TFile fout(foutName,"recreate");
    h1ee->Write("h1ee_orig");
    h1mm->Write("h1mm_orig");
    h1ll->Write("h1ll_orig");
    h2covEE->Write("h2eeCov_orig");
    h2covMM->Write("h2mmCov_orig");
    h2covLL->Write("h2llCov_orig");
    writeMsg("packed materials_SMP14003/1D");
    writeTimeTag();
    fout.Close();
    std::cout << "file <" << fout.GetName() << "> saved\n";
  }
}

// ---------------------------------------------------------------
