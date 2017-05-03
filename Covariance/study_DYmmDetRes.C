#include "inputs.h"
#include "../RooUnfold/src/RooUnfoldResponse.h"
#include "../RooUnfold/src/RooUnfoldBayes.h"
#include "../RooUnfold/src/RooUnfoldInvert.h"
#include "crossSection.h"

//#include "study_FSRResp.C" // for case definitions

typedef enum { _mu76Xfsr=0,
	       _mu76XfsrKL=1,
	       _mu76XdetRes4p2=2,
	       _mu76XdetRes4p3=3,
	       _mu76XdetRes4p2KL=4,
	       _mu76XdetRes4p3KL=5,
	       _elRidhiFSR=6,
	       _elRidhiDetRes=7,
	       _elV2SkimFSR=8,
	       _elV2SkimDetRes=9,
	       _elV3FSR=10,
	       _elV3DetRes=11
} TStudyCase_t;

typedef enum { _unfoldBayes=0, _unfoldInvert } TUnfoldMethod_t;


void study_DYmmDetRes(int nSample, int studyCase=_mu76XdetRes4p3)
{
  closeCanvases(5);

  int baselineTest=1;
  int iLepton=1; // 0 - el, 1 - mu, 2 - lepton

  // main result
  TString fnameMain="cs_DYmm_13TeV_76X_cs.root";
  TString nameh1obsData;
  TString nameh1bkgEst;
  TString nameh1preUnfData="h1PostFSR";
  TString nameh1unfData="h1PreFSR";

  // Response matrix
  TString fnameResp="dymm_test_Mu76X.root";
  TString nameResp="rooUnf_fsrResp";
  TString fnameMC=fnameResp;
  TString nameh1measMC="h1_postFsr_Mweighted";
  TString nameh1trueMC="h1_preFsr_Mweighted";
  TString nameh2migrationMatrix=""; // define only if nameResp is empty
  int hasOverflows=0;
  int nIters=20;
  double mcScale=1.;
  TUnfoldMethod_t unfoldMethod=_unfoldBayes;

  // Covariance estimation due to randomized response matrix
  TString covFName="cov_mumu_varFSRRes_10000.root";
  TString nameh2cov="h2Cov";

  // output
  TString outDirTag="_mu76Xfsr";


  if (studyCase==_mu76Xfsr) {
    // all is set
  }
  else if ((studyCase==_mu76XdetRes4p2) || (studyCase==_mu76XdetRes4p3)) {
    iLepton=1;
    baselineTest=1;
    TString hltV= (studyCase==_mu76XdetRes4p2) ? "4p2" : "4p3";
    TString runPeriod= (studyCase==_mu76XdetRes4p2) ? "A" : "B";
    std::cout << " studyCase=" << hltV << "(runPeriod=" << runPeriod << ")\n";
    TString dataPath="./";
    fnameMain=dataPath + "cs_DYmm_13TeVMuApproved_cs" + runPeriod + ".root";
    nameh1preUnfData="h1Signal";
    nameh1unfData="h1Unf";
    nameh1obsData="h1Yield";
    nameh1bkgEst="h1Bkg";
    fnameResp="dymm_test_RECO_MuApproved.root";
    nameResp="rooUnf_detResRespPU";
    fnameMC=fnameResp;
    nameh1measMC="h1recoSel_MWPU";
    nameh1trueMC="h1_postFsrInAccSel_MWPU";
    nameh2migrationMatrix="";
    double lumiTot=2832.673;
    double lumiA= 865.919;
    double lumiB= lumiTot-lumiA;
    mcScale= (studyCase==_mu76XdetRes4p2) ? lumiA/lumiTot : lumiB/lumiTot;
    hasOverflows=0;
    nIters=17;
    covFName="cov_mumu_varDetRes_10000.root";
    nameh2cov="h2Cov";
    outDirTag="_mu76XdetRes"+hltV;
    unfoldMethod=_unfoldInvert;
  }

  TH1D *h1obsData= loadHisto(fnameMain,nameh1obsData,"h1obsData",h1dummy);
  TH1D *h1bkgEst= loadHisto(fnameMain,nameh1bkgEst,"h1bkgEst",h1dummy);
  TH1D *h1measData_file= loadHisto(fnameMain,nameh1preUnfData,
				   "h1measData_file",h1dummy);
  TH1D *h1unfData_file= loadHisto(fnameMain,nameh1unfData,
				  "h1unfData_file",h1dummy);
  if (!h1obsData || !h1bkgEst || !h1measData_file || !h1unfData_file) {
    std::cout << "failed to load histos from fnameMain=<" << fnameMain << ">\n";
    return;
  }
  if (0) {
    printRatio(h1measData_file,h1unfData_file);
    return;
  }

  if (0) {
    TH1D *h1measData_wBkgErr= cloneHisto(h1obsData,
			    "h1measData_wBkgErr","h1measData_wBkgErr");
    h1measData_wBkgErr->Add(h1bkgEst,-1);
    TH1D *h1measData_noBkgErr= cloneHisto(h1obsData,
				  "h1measData_noBkgErr","h1measData_noBkgErr");
    removeError(h1bkgEst);
    h1measData_noBkgErr->Add(h1bkgEst,-1);

    hsGreen.SetStyle(h1measData_noBkgErr);
    hsBlue.SetStyle(h1measData_wBkgErr);
    hsRed.SetStyle(h1measData_file);
    h1measData_file->SetLineWidth(4);

    plotHisto(h1measData_file,"cSigErrCheck",1,1,"LPE1","sig.data (file)");
    plotHistoSame(h1measData_wBkgErr,"cSigErrCheck","LPE1","sig (w.bkg err)");
    plotHistoSame(h1measData_noBkgErr,"cSigErrCheck","LPE1","sig (no bkg err)");

    TH1D *h1measData_unc= errorAsCentral(h1measData_file);
    TH1D *h1measData_wBkgErr_unc= errorAsCentral(h1measData_wBkgErr);
    TH1D *h1measData_noBkgErr_unc= errorAsCentral(h1measData_noBkgErr);
    plotHisto(h1measData_unc,"cSigErrCheck_Unc",1,1,"LP","sig.data (file)");
    plotHistoSame(h1measData_wBkgErr_unc,"cSigErrCheck_Unc","LP","sig (w.bkg err)");
    plotHistoSame(h1measData_noBkgErr_unc,"cSigErrCheck_Unc","LP","sig (no bkg err)");

    printRatio(h1measData_file,h1measData_wBkgErr);
    printRatio(h1measData_file,h1measData_noBkgErr);
    return;
  }

  removeError(h1bkgEst);
  TH1D *h1measData= cloneHisto(h1obsData,"h1measData","h1measData");
  h1measData->Add(h1bkgEst,-1);

  if (0) {
    hsRed.SetStyle(h1measData_file);
    h1measData_file->SetLineWidth(4);
    plotHisto(h1measData_file,"cChkMeas2",1,1,"LPE1","h1measData_file");
    plotHistoSame(h1measData,"cChkMeas2","LPE1","h1measData (no bkg.err)");
    return;
  }


  // Load unfolding matrix
  RooUnfoldResponse *rooResp= NULL;
  TH1D* h1measMC= loadHisto(fnameMC,nameh1measMC,"h1measMC",h1dummy);
  TH1D* h1trueMC= loadHisto(fnameMC,nameh1trueMC,"h1trueMC",h1dummy);
  if (!h1measMC || !h1trueMC) {
    std::cout << "failed to load histos from fnameMC=<" << fnameMC << ">\n";
    return;
  }
  if (mcScale!=double(1.0)) {
    h1measMC->Scale(mcScale);
    h1trueMC->Scale(mcScale);
  }
  if (nameResp.Length()) {
    rooResp= loadRooUnfoldResponse(fnameResp,nameResp,"DRR");
  }
  else {
    std::cout << "nameResp is not provided\n";
    return;
  }
  rooResp->UseOverflow(false);

  // Check MC unfolding on MC
  RooUnfold *fsrBayesMC=NULL;
  if (unfoldMethod==_unfoldBayes) {
    fsrBayesMC= new RooUnfoldBayes( rooResp, h1measMC, nIters, false, "fsrBayesMC" );
  }
  else {
    fsrBayesMC= new RooUnfoldInvert( rooResp, h1measMC, "fsrBayesMC" );
  }
  TH1D *h1unfBayesMC= cloneHisto(fsrBayesMC->Hreco(),"h1UnfMC",
				 TString("UnfMC;") +
				 h1trueMC->GetXaxis()->GetTitle() + TString(";")
				 + "unfolded", h1trueMC);

  histoStyle(h1measMC, kBlue,24);
  histoStyle(h1trueMC, kRed,26);
  histoStyle(h1unfBayesMC, kGreen, 5);

  plotHisto(h1measMC, "cCmpMC", 1,1, "LPE1", "measMC");
  plotHistoSame(h1trueMC, "cCmpMC", "LPE", "trueMC");
  plotHistoSame(h1unfBayesMC, "cCmpMC", "LPE", "bayesUnfMC");
  printRatio(h1trueMC,h1unfBayesMC);
  //return;

  // do unfolding on our h1measData
  TH1D *h1unfData=NULL;
  if (1) {
    RooUnfold *fsrBayes_tmp=NULL;
    if (unfoldMethod==_unfoldBayes) {
      fsrBayes_tmp= new RooUnfoldBayes( rooResp, h1measData, nIters, false,
					"fsrBayes_tmp" );
    }
    else {
      fsrBayes_tmp= new RooUnfoldInvert( rooResp, h1measData, "fsrBayes_tmp" );
    }
    h1unfData= cloneHisto(fsrBayes_tmp->Hreco(),"h1UnfData",
			  TString("UnfData;") +
			  h1measData->GetXaxis()->GetTitle() + TString(";")
			  + "unfolded", h1measData);
  }

  // verify the unfolding
  if (0) {
    hsRed.SetStyle(h1unfData_file);
    plotHistoError(h1unfData_file,"cChkUnf2err",1,1,"LP","h1unfData_file");
    plotHistoErrorSame(h1unfData,"cChkUnf2err","LP","h1unfData");

    plotHisto(h1unfData_file,"cChkUnf2",1,1,"LPE1","h1unfData_file");
    plotHistoSame(h1unfData,"cChkUnf2","LPE1","h1unfData");

    printRatio(h1unfData_file,h1unfData);
    return;
  }

  TString axisLabel_unf= niceMassAxisLabel(iLepton,"unf");
  TString cov_axisLabel_unf= axisLabel_unf + TString(";") + axisLabel_unf;
 
   TH2D *h2unfCov=NULL; // covariance from my toys
  if (baselineTest) {
    std::vector<TH1D*> vecPreFSR;
    int nonNegative=0;
    for (int i=0; i<nSample; i++) {
      TString nameH1Var= TString("hVar") + Form("%d",i);
      TH1D *h1_in= cloneHisto(h1obsData,nameH1Var,nameH1Var);
      randomizeWithinErr(h1obsData, h1_in, nonNegative);
      h1_in->Add(h1bkgEst,-1);
      RooUnfold *fsrBayes= NULL;
      if (unfoldMethod==_unfoldBayes) {
	fsrBayes= new RooUnfoldBayes( rooResp, h1_in, nIters, false, "bayes_"+TString(Form("%d",i)));
      }
      else {
	fsrBayes= new RooUnfoldInvert( rooResp, h1_in, "invert_"+TString(Form("%d",i)) );
      }
      TString outTag=TString("h1_out_") + Form("%d",i);
      TH1D *h1_out= cloneHisto(fsrBayes->Hreco(),outTag,outTag, h1unfData_file);
      vecPreFSR.push_back(h1_out);
      delete fsrBayes;
      delete h1_in;
    }
    TH1D *h1avgUnf=NULL;
    deriveCovariance(vecPreFSR,"toys","toys",&h1avgUnf,&h2unfCov);
    //TString axisLabel_unf= niceMassAxisLabel(iLepton,"unf");
      //"M_{\\ell\\!\\ell\\!,\\,\\text{unf}}\\text{ [GeV]}";
    h2unfCov->GetXaxis()->SetTitle(axisLabel_unf);
    h2unfCov->GetYaxis()->SetTitle(axisLabel_unf);
    h2unfCov->GetYaxis()->SetTitleOffset(1.8);

    if (0 && (vecPreFSR.size()>0)) {
      plotHisto(vecPreFSR[0],"cToy",1,1,"LPE1");
      for (unsigned int i=1; (i<50) && (i<vecPreFSR.size()); i++) {
	plotHistoSame(vecPreFSR[i],"cToy","LPE");
      }
    }
  }

  PlotCovCorrOpt_t ccOpt;
  ccOpt.yTitleOffset= 1.8;
  if (h2unfCov) plotCovCorr(h2unfCov,"cCovUnf",ccOpt,NULL);


  // chi^2 calculation
  TVectorD vTrueMC( convert2vec(h1trueMC) );
  TVectorD vMeasMC( convert2vec(h1measMC) );
  //TVectorD vMeas( convert2vec(h1measData_file) );
  //TVectorD vUnf( convert2vec(h1unfData_file) );

  TMatrixD mSmear( convert2mat( rooResp->Hresponse() ) );
  for (int ic=0; ic<mSmear.GetNcols(); ic++) {
    double sum=0;
    for (int ir=0; ir<mSmear.GetNrows(); ir++) {
      sum += mSmear(ir,ic);
    }
    for (int ir=0; ir<mSmear.GetNrows(); ir++) {
      mSmear(ir,ic) /= sum;
    }
  }

  if (0) {
    TMatrixD mSmear2(rooResp->Mresponse());
    TMatrixD mSmearDiff(mSmear);
    mSmearDiff-= mSmear2;
    mSmearDiff.Print();
    return;
  }
  else {
    mSmear = rooResp->Mresponse();
  }

  TVectorD vSmearCheck= mSmear * vTrueMC;
  printRatio("vMeasMC", vMeasMC, "vSmearCheck", vSmearCheck);
  //return;

  TMatrixD mSmearInv(TMatrixD::kInverted,mSmear);
  if (1) {
    TVectorD vUnsmearCheck= mSmearInv * vMeasMC;
    printRatio("vTrueMC",vTrueMC, "vUnsmearCheck",vUnsmearCheck);
    if (0) {
      TMatrixD unityTest( mSmearInv * mSmear );
      unityTest.Print("unityTest");
    }
  }

  if (0) {
    TMatrixD mResp ( rooResp->Mresponse() );
    TMatrixD diff ( mResp - mSmear );
    diff.Print();
  }

  if (!h2unfCov) {
    std::cout << "h2unfCov is null\n";
    return;
  }

  // check properties of mSmear matrix
  if (1) {
    std::cout << "\n\tREPLACING mSmear by mSmearInv\n\n"; mSmear=mSmearInv;
    std::cout << "check diagonal-dominant properties of mSmear matrix (row)\n";
    int isDiagDom=1;
    double smallestRatio=100.;
    for (int ir=0; ir<mSmear.GetNrows(); ir++) {
      double sumAbsNonDiag=0, sumRow=0;
      for (int ic=0; ic<mSmear.GetNcols(); ic++) {
	if (ir!=ic) sumAbsNonDiag+= fabs(mSmear(ir,ic));
	sumRow+= mSmear(ir,ic);
      }
      double dg= fabs(mSmear(ir,ir));
      std::cout << "ir=" << ir << "  ";
      std::cout << "dgEntry=" << mSmear(ir,ir) << ", sumRow=" << sumRow << "  ";
      std::cout << dg << " " << sumAbsNonDiag << " ";
      if (dg>sumAbsNonDiag) {
	std::cout << "ok  " << dg/sumAbsNonDiag << "\n";
	if (smallestRatio> dg/sumAbsNonDiag) smallestRatio=dg/sumAbsNonDiag;
      }
      else {
	std::cout << "!!\n";
	isDiagDom=0;
      }
    }
    std::cout << "isDiagDom=" << isDiagDom << "\n";
    if (isDiagDom) std::cout << " -- smallest ratio=" << smallestRatio << "\n";
    plotCovCorr(mSmear,h1measData_file,"h2mSmear","cMSmear");
    return;
  }
  if (1) {
    std::cout << "\n\tREPLACING mSmear by mSmearInv\n\n"; mSmear=mSmearInv;
    std::cout << "check diagonal-dominant properties of mSmear matrix (column)\n";
    int isDiagDom=1;
    double smallestRatio=100.;
    for (int ic=0; ic<mSmear.GetNcols(); ic++) {
      double sumAbsNonDiag=0, sumCol=0;
      for (int ir=0; ir<mSmear.GetNrows(); ir++) {
	if (ir!=ic) sumAbsNonDiag+= fabs(mSmear(ir,ic));
	sumCol+= mSmear(ir,ic);
      }
      double dg= fabs(mSmear(ic,ic));
      std::cout << "ic=" << ic << "  ";
      std::cout << "dgEntry=" << mSmear(ic,ic) << ", sumCol=" << sumCol << "  ";
      std::cout << dg << " " << sumAbsNonDiag << " ";
      if (dg>sumAbsNonDiag) {
	std::cout << "ok  " << dg/sumAbsNonDiag << "\n";
	if (smallestRatio> dg/sumAbsNonDiag) smallestRatio=dg/sumAbsNonDiag;
      }
      else {
	std::cout << "!!\n";
	isDiagDom=0;
      }
    }
    std::cout << "isDiagDom=" << isDiagDom << "\n";
    if (isDiagDom) std::cout << " -- smallest ratio=" << smallestRatio << "\n";

    plotCovCorr(mSmear,h1measData_file,"h2mSmear","cMSmear");
    return;
  }

  TH1D *h1uncFromCov= uncFromCov(h2unfCov);
  TH1D *h1unfCheck= cloneHisto(h1measData_file,"h1unfCheck","h1unfCheck");

  TH1D *h1uncFromProp= cloneHisto(h1uncFromCov,"h1uncFromProp","h1unfFromProp");
  h1uncFromProp->Scale(1e-4);
  h1uncFromProp->Reset();
  TH1D *h1uncSqrtNu= cloneHisto(h1uncFromCov,"h1uncSqrtNu","h1uncSqrtNu");
  h1uncSqrtNu->Scale(1e-4);
  h1uncSqrtNu->Reset();

  for (int ibin=1; ibin<=h1uncFromProp->GetNbinsX(); ibin++) {
    double sumUnf=0;
    double sumUncSqr=0;
    for (int jbin=1; jbin<=h1uncFromProp->GetNbinsX(); jbin++) {
      double tij= mSmearInv(ibin-1,jbin-1);
      double unc= h1measData->GetBinError(jbin);
      sumUncSqr+= tij*tij*unc*unc;
      sumUnf+= tij * h1measData->GetBinContent(jbin);
    }
    h1unfCheck->SetBinContent(ibin, sumUnf);
    h1uncFromProp->SetBinContent(ibin, sqrt(sumUncSqr));
    h1uncSqrtNu->SetBinContent(ibin, sqrt(h1unfData->GetBinContent(ibin)));
  }

  hsColor6.SetStyle(h1uncFromCov);
  hsViolet.SetStyle(h1uncFromProp);
  hsViolet.SetStyle(h1unfCheck);
  hsGreen.SetStyle(h1uncSqrtNu);
  hsRed.SetStyle(h1measData_file);

  printHisto(h1uncFromCov);
  printHisto(h1uncFromProp);

  plotHisto(h1unfData,"cUnfCheck",1,1,"LPE","unf data");
  //plotHistoAuto(h1unfData_file,"cUnfCheck",1,1,"LPE","unf data (file)");
  plotHistoSame(h1unfCheck,"cUnfCheck","LPE","unf check");
  hsRed.SetStyle(h1measData);
  plotHistoSame(h1measData,"cUnfCheck","LPE","signal data");
  //plotHistoSame(h1measData_file,"cUnfCheck","LPE","signal data");

  plotHisto(h1uncFromCov,"cUncertaintyCheck",1,1,"LPE","unc from cov");
  plotHistoSame(h1uncFromProp,"cUncertaintyCheck","LPE","unc by prop.");
  plotHistoSame(h1uncSqrtNu,"cUncertaintyCheck","LPE","sqrt(N_(u))");
  if (0 && h1obsData) {
    TH1D *h1uncSqrtNobs= cloneHisto(h1obsData,"h1uncSqrtNobs","h1uncSqrtNobs");
    for (int ibin=1; ibin<=h1obsData->GetNbinsX(); ibin++) {
      h1uncSqrtNobs->SetBinContent(ibin, sqrt(h1obsData->GetBinContent(ibin)));
    }
    hsColor46.SetStyle(h1uncSqrtNobs);
    plotHistoSame(h1uncSqrtNobs,"cUncertaintyCheck","LPE","sqrt(N_(obs))");
  }
  if (0 && h1unfData) {
    TH1D *h1unfUnc= errorAsCentral(h1unfData);
    hsRed.SetStyle(h1unfUnc);
    plotHistoSame(h1unfUnc,"cUncertaintyCheck","LPE",
		  "unc from unfData");
  }
  if (0 && h1measData) {
    TH1D *h1unfUnc= errorAsCentral(h1measData);
    hsBlack.SetStyle(h1unfUnc);
    plotHistoSame(h1unfUnc,"cUncertaintyCheck","LPE","unc from sigData");
  }
}
