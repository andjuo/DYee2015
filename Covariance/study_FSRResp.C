// Small macro for testing purposes

#include "inputs.h"
#include "../RooUnfold/src/RooUnfoldResponse.h"
#include "../RooUnfold/src/RooUnfoldBayes.h"
#include "../RooUnfold/src/RooUnfoldInvert.h"
#include "crossSection.h"
//#include <TMatrixD.h>
#include <TDecompSVD.h>

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

// --------------------------------------------------------------------

void study_FSRResp(int nSample, int studyCase=_elV2SkimFSR, int doSave=0,
		   int nIters_user=-1)
{
  if (0) {
    TMatrixD A(5,5);
    for (int i=0; i<5; i++)
      for (int j=0; j<5; j++)
	A(i,j) = i*10+j;
    TMatrixD B(submatrix(A,1,A.GetNrows()-1));
    A.Print();
    B.Print();
    return;
  }

  int baselineTest=1;
  int iLepton=1; // 0 - el, 1 - mu, 2 - lepton

  // main result
  TString fnameMain="cs_DYmm_13TeV_76X_cs.root";
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
    fnameMain=dataPath + "cs_DYmm_13TeV_76X_cs" + runPeriod + ".root";
    nameh1preUnfData="h1Signal";
    nameh1unfData="h1Unf";
    fnameResp="dymm_test_RECO_Mu76X_withOverflow.root";
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
  }
  else if ((studyCase==_mu76XdetRes4p2KL) || (studyCase==_mu76XdetRes4p3KL)) {
    iLepton=1;
    baselineTest=1;
    TString hltV= (studyCase==_mu76XdetRes4p2KL) ? "4p2KL" : "4p3KL";
    TString runPeriod= (studyCase==_mu76XdetRes4p2KL) ? "A" : "B";
    std::cout << " studyCase=" << hltV << "(runPeriod=" << runPeriod << ")\n";
    TString dataPath="./";
    fnameMain=dataPath + "cs_DYmm_13TeV_76X_cs" + runPeriod + ".root";
    nameh1preUnfData="h1Signal";
    nameh1unfData="h1Unf";
    fnameResp= fnameMain;
    nameResp="detRes";
    fnameMC="dymm_test_RECO_Mu76X_withOverflow.root";
    nameh1measMC="h1recoSel_MWPU";
    nameh1trueMC="h1_postFsrInAccSel_MWPU";
    nameh2migrationMatrix="";
    double lumiTot=2832.673;
    double lumiA= 865.919;
    double lumiB= lumiTot-lumiA;
    mcScale= (studyCase==_mu76XdetRes4p2KL) ? lumiA/lumiTot : lumiB/lumiTot;
    hasOverflows=0;
    nIters=17;
    //mcScale=2832.673;
    covFName=""; //"cov_mumu_varDetRes_10000.root";
    nameh2cov="h2Cov";
    outDirTag="_mu76XdetRes"+hltV;
  }
  else if (studyCase==_mu76XfsrKL) {
    std::cout << " studyCase= _mu76XfsrKL\n";
    iLepton=1;
    fnameMain="cs_DYmm_13TeV_76X_cs.root";
    nameh1preUnfData="h1PostFSR";
    nameh1unfData="h1PreFSR";
    fnameResp= "cs_DYmm_13TeV_76X_csA.root";
    nameResp="FSRRes";
    fnameMC="dymm_test_Mu76X.root";
    nameh1measMC= "h1_postFsr_Mweighted";
    nameh1trueMC= "h1_preFsr_Mweighted";
    hasOverflows=0;
    nIters=17;
    //mcScale=2832.673;
    covFName="cov_mumu_varFSRRes_10000.root" ;
    nameh2cov="h2Cov";
    outDirTag="_mu76XfsrKL";
  }
  else if (studyCase==_elRidhiFSR) {
    iLepton=0;
    std::cout << " studyCase= _elRidhiFSR\n";
    std::cout << "\n\tdata is not real!!\n";
    baselineTest=0;
    TString dataPath="dir-DYee-Ridhi20160603/";
    fnameMain=dataPath + "DYEE_FSRCorr_PosttoPre.root";
    nameh1preUnfData="postFSR_Mass";
    nameh1unfData="preFSR_Mass";
    fnameResp=fnameMain;
    nameResp="";
    fnameMC=fnameMain;
    nameh1measMC="postFSR_Mass";
    nameh1trueMC="preFSR_Mass";
    nameh2migrationMatrix="RespMatrix_FSR";
    hasOverflows=0;
    covFName="";
    nameh2cov="";
    outDirTag="_elFSRRes";
  }
  else if (studyCase==_elRidhiDetRes) {
    iLepton=0;
    std::cout << " studyCase= _elRidhiDetRes\n";
    std::cout << "\n\tdata is not real!!\n";
    baselineTest=0;
    TString dataPath="dir-DYee-Ridhi20160603/";
    fnameMain=dataPath + "DYEE_DetRes_RecotoPost.root";
    nameh1preUnfData="reco_EEMass";
    nameh1unfData="gen_EEMass";
    fnameResp=fnameMain;
    nameResp="";
    fnameMC=fnameMain;
    nameh1measMC="reco_EEMass";
    nameh1trueMC="gen_EEMass";
    nameh2migrationMatrix="RespMatrix";
    hasOverflows=0;
    covFName="";
    nameh2cov="";
    outDirTag="_elDetRes";
  }
  else if (studyCase==_elV2SkimFSR) {
    iLepton=0;
    std::cout << "studyCase= _elV2SkimFSR\n";
    baselineTest=0;
    mcScale=2316.97;
    // measured postFSR and preFSR distributions
    fnameMain="cs_DYee_13TeV_El2.root";
    nameh1preUnfData="h1UnfRhoEffAcc";
    nameh1unfData="h1PreFSR";
    // FSR response matrix
    fnameResp=fnameMain;
    nameResp="FSRRes";
    TString simFName="dyee_test_El2skim.root";
    simFName="dyee_test_El2skim_excludeGap.root";
    if (0) { // closure
      mcScale=1.;
      fnameResp=simFName;
      nameResp="rooUnf_fsrResp";
    }
    nIters=15;
    // simulated distributions
    fnameMC=simFName;
    nameh1measMC="h1_postFsr_Mweighted";
    nameh1trueMC="h1_preFsr_Mweighted";
    nameh2migrationMatrix="";
    hasOverflows=0;
    // response matrix randomization
    covFName="";
    nameh2cov="";
    outDirTag="_elV2SkimFSR";
  }
  else if (studyCase==_elV2SkimDetRes) {
    iLepton=0;
    std::cout << "studyCase= _elV2SkimDetRes\n";
    baselineTest=0;
    mcScale=2316.97;
    // measured postFSR and preFSR distributions
    fnameMain="cs_DYee_13TeV_El2.root";
    nameh1preUnfData="h1Signal";
    nameh1unfData="h1Unf";
    // FSR response matrix
    fnameResp=fnameMain;
    nameResp="detRes";
    if (0) { // closure
      mcScale=1.;
      fnameResp="dyee_test_RECO_El2skim.root";
      nameResp="rooUnf_detResRespPU";
    }
    nIters=15;
    // simulated distributions
    fnameMC="dyee_test_RECO_El2skim.root";
    nameh1measMC="h1recoSel_MWPU";
    nameh1trueMC="h1_postFsrInAccSel_MWPU";
    nameh2migrationMatrix="";
    hasOverflows=0;
    // response matrix randomization
    covFName="";
    nameh2cov="";
    outDirTag="_elV2SkimDetRes";
  }
  else if (studyCase==_elV3FSR) {
    iLepton=0;
    std::cout << "studyCase= _elV3FSR\n";
    baselineTest=0;
    mcScale=2316.97;
    // measured postFSR and preFSR distributions
    fnameMain="cs_DYee_13TeV_El3.root";
    nameh1preUnfData="h1UnfRhoEffAcc";
    nameh1unfData="h1PreFSR";
    // FSR response matrix
    fnameResp=fnameMain;
    nameResp="FSRRes";
    TString simFName="dyee_test_dressed_El3.root";
    //simFName="dyee_test_El2skim_excludeGap.root";
    if (0) { // closure
      mcScale=1.;
      fnameResp=simFName;
      nameResp="rooUnf_fsrResp";
    }
    nIters=15;
    // simulated distributions
    fnameMC=simFName;
    nameh1measMC="h1_postFsr_Mweighted";
    nameh1trueMC="h1_preFsr_Mweighted";
    nameh2migrationMatrix="";
    hasOverflows=0;
    // response matrix randomization
    covFName="";
    nameh2cov="";
    outDirTag="_elV3FSR";
  }
  else if (studyCase==_elV3DetRes) {
    iLepton=0;
    std::cout << "studyCase= _elV3SkimDetRes\n";
    baselineTest=0;
    mcScale=2316.97;
    // measured postFSR and preFSR distributions
    fnameMain="cs_DYee_13TeV_El3.root";
    nameh1preUnfData="h1Signal";
    nameh1unfData="h1Unf";
    // FSR response matrix
    fnameResp=fnameMain;
    nameResp="detRes";
    if (0) { // closure
      mcScale=1.;
      fnameResp="dyee_test_dressed_El3.root";
      nameResp="rooUnf_detResRespPU";
    }
    nIters=15;
    // simulated distributions
    fnameMC="dyee_test_dressed_El3.root";
    nameh1measMC="h1recoSel_MWPU";
    nameh1trueMC="h1_postFsrInAccSel_MWPU";
    nameh2migrationMatrix="";
    hasOverflows=0;
    // response matrix randomization
    covFName="";
    nameh2cov="";
    outDirTag="_elV3DetRes";
  }

  if (nIters_user>0) {
    std::cout << "overriding the number of iterations to "<<nIters_user << "\n";
    nIters=nIters_user;
  }

  // Load data
  TH1D *h1measData= loadHisto(fnameMain,nameh1preUnfData,"h1measData",h1dummy);
  TH1D *h1unfData= loadHisto(fnameMain,nameh1unfData,"h1unfData",h1dummy);
  if (!h1measData || !h1unfData) {
    std::cout << "failed to load histos from fnameMain=<" << fnameMain << ">\n";
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
    if (nameh2migrationMatrix.Length()==0) {
      std::cout << "neither rooUnfoldResponse nor migration matrix name is given\n";
      return;
    }
    TH2D* h2mig= loadHisto(fnameMC,nameh2migrationMatrix,"h2mig",h2dummy);
    if (!h2mig) {
      std::cout << "failed to load the migration matrix\n";
      return;
    }
    rooResp= new RooUnfoldResponse(h1measMC,h1trueMC,h2mig,"DRR");
  }
  if (!rooResp) {
    std::cout << "failed to get RooUnfoldResponse\n";
    return;
  }
  rooResp->UseOverflow(false);

  // Check distributions of the response objects
  if (1) {
    TH1D *h1measResp= cloneHisto((TH1D*)rooResp->Hmeasured(),"h1measResp","measResp");
    TH1D *h1trueResp= cloneHisto((TH1D*)rooResp->Htruth(),"h1trueResp","trueResp");
    histoStyle(h1measResp,kBlue,5,1);
    histoStyle(h1trueResp,kBlue,5,1);
    plotHisto(h1measMC,"cMeasMCCheck",1,1,"LPE","h1measMC");
    plotHistoSame(h1measResp,"cMeasMCCheck","LPE","h1measResp");
    plotHisto(h1trueMC,"cTrueMCCheck",1,1,"LPE","h1trueMC");
    plotHistoSame(h1trueResp,"cTrueMCCheck","LPE","h1trueResp");
    std::cout << "integrals of measMC and measResp: " << h1measMC->Integral()
	      << " " << h1measResp->Integral() << "\n";
    printRatio(h1measMC,h1measResp);
    std::cout << "integrals of trueMC and trueResp: " << h1trueMC->Integral()
	      << " " << h1trueResp->Integral() << "\n";
    printRatio(h1trueMC,h1trueResp);
  }

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

  // Compare data unfolding
  RooUnfold *fsrBayesDataMy=NULL;
  if (unfoldMethod==_unfoldBayes) {
    fsrBayesDataMy= new RooUnfoldBayes( rooResp, h1measData, nIters, false, "fsrBayesDataMy" );
  }
  else {
    fsrBayesDataMy= new RooUnfoldInvert( rooResp, h1measData, "fsrBayesDataMy" );
  }
  TMatrixD mBayesMeasCovAll(fsrBayesDataMy->GetMeasuredCov());
  TMatrixD mBayesMeasCov(submatrix(mBayesMeasCovAll,hasOverflows,mBayesMeasCovAll.GetNrows()-hasOverflows));
  //TH2D *h2bayesMeasCov_byBinNum= new TH2D(mBayesMeasCov);
  //h2bayesMeasCov_byBinNum->SetDirectory(0);
  TString axisLabel_bmc= niceMassAxisLabel(iLepton,"preUnf");
  TH2D *h2bayesMeasCov= convert2histo(mBayesMeasCov,h1measData,"h2bayesMeasCov",
			      "bayesMeasCov;" + axisLabel_bmc + TString(";") + axisLabel_bmc);
  h2bayesMeasCov->GetYaxis()->SetTitleOffset(1.8);
  TMatrixD mBayesRecoCovAll(fsrBayesDataMy->Ereco(RooUnfold::kCovariance));
  TMatrixD mBayesRecoCov(submatrix(mBayesRecoCovAll,hasOverflows,mBayesRecoCovAll.GetNrows()-hasOverflows));
  //TH2D *h2bayesRecoCov_byBinNum= new TH2D(mBayesRecoCov);
  //h2bayesRecoCov_byBinNum->SetDirectory(0);
  TString axisLabel_unf= niceMassAxisLabel(iLepton,"unf");
  TString cov_axisLabel_unf= axisLabel_unf + TString(";") + axisLabel_unf;
  TH2D *h2bayesRecoCov= convert2histo(mBayesRecoCov,h1unfData,"h2bayesRecoCov",
				      "bayesRecoCov;" + cov_axisLabel_unf);
  h2bayesRecoCov->GetYaxis()->SetTitleOffset(1.8);

  TH1D *h1unfBayesMy= cloneHisto(fsrBayesDataMy->Hreco(), "h1unfBayesMy",
				 TString("UnfDataMy;") +
				 h1unfData->GetXaxis()->GetTitle() +
				 TString("; unfolded(my)"), h1unfData);

  histoStyle(h1unfBayesMy, kBlue, 25);
  histoStyle(h1unfData, kBlack, 5);
  histoStyle(h1measData, kGreen+1, 27);

  plotHisto(h1unfBayesMy, "cCmpData", 1,1, "LPE1", "myBayesUnfData");
  plotHistoSame(h1unfData, "cCmpData", "LPE", "unfData");
  plotHistoSame(h1measData, "cCmpData", "LPE", "measData");
  printRatio(h1unfBayesMy, h1unfData);
  //printHisto(h1unfBayesMy);
  //printHisto(h1unfData);

  TH2D *h2unfCov=NULL; // covariance from my toys
  if (baselineTest) {
    std::vector<TH1D*> vecPreFSR;
    int nonNegative=1;
    for (int i=0; i<nSample; i++) {
      TString nameH1Var= TString("hVar") + Form("%d",i);
      TH1D *h1_in= cloneHisto(h1measData,nameH1Var,nameH1Var);
      randomizeWithinErr(h1measData, h1_in, nonNegative);
      RooUnfold *fsrBayes= NULL;
      if (unfoldMethod==_unfoldBayes) {
	fsrBayes= new RooUnfoldBayes( rooResp, h1_in, nIters, false, "bayes_"+TString(Form("%d",i)));
      }
      else {
	fsrBayes= new RooUnfoldInvert( rooResp, h1_in, "invert_"+TString(Form("%d",i)) );
      }
      TString outTag=TString("h1_out_") + Form("%d",i);
      TH1D *h1_out= cloneHisto(fsrBayes->Hreco(),outTag,outTag, h1unfData);
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

  // compare covariances
  plotCovCorr(h2bayesMeasCov,"cCovBayesMeas",ccOpt,NULL);
  plotCovCorr(h2bayesRecoCov,"cCovBayesReco",ccOpt,NULL);
  if (h2unfCov) plotCovCorr(h2unfCov,"cCovUnf",ccOpt,NULL);

  // covariance of the cross sections due to the limited MC statistics
  // for unfolding
  TH2D* h2covRespMStat= NULL;
  if (covFName.Length()>0) {
    h2covRespMStat=loadHisto(covFName, nameh2cov, "h2covRespMStat", h2dummy);
    if (!h2covRespMStat) {
      std::cout << "failed to load histos from covFName=<" << covFName << ">\n";
      return;
    }
    TString axisLabel_mstat= niceMassAxisLabel(iLepton,"");
    h2covRespMStat->GetXaxis()->SetTitle(axisLabel_mstat);
    h2covRespMStat->GetYaxis()->SetTitle(axisLabel_mstat);
    h2covRespMStat->GetYaxis()->SetTitleOffset(1.8);

    if (0) {
      TMatrixD mCovRespMStat( convert2mat(h2covRespMStat ) );
      TH2D *h2Chk= new TH2D(mCovRespMStat);
      plotCovCorr(h2covRespMStat,"cCovRespMStat",ccOpt,NULL);
      plotCovCorr(h2Chk, "cCovChk",ccOpt,NULL);
      printHisto(h2Chk);
      mCovRespMStat.Print();
      return;
    }

    plotCovCorr(h2covRespMStat,"cCovRespMStat",ccOpt,NULL);
  }

  // compare uncertainties from the covariance
  TH1D* h1uncFromUnfCov= (h2unfCov) ? uncFromCov(h2unfCov,NULL) : NULL;
  TH1D* h1uncFromBayesMeasCov= cloneHisto(h1trueMC,"h1uncFromBayesMeasCov","h1uncFromBayesMeasCov");
  h1uncFromBayesMeasCov->SetStats(0);
  copyContents(h1uncFromBayesMeasCov, uncFromCov(h2bayesMeasCov,NULL));
  TH1D* h1uncFromBayesRecoCov= cloneHisto(h1trueMC,"h1uncFromBayesRecoCov","h1uncFromBayesRecoCov");
  copyContents(h1uncFromBayesRecoCov, uncFromCov(h2bayesRecoCov,NULL));
  histoStyle(h1uncFromBayesMeasCov, kBlue, kPlus);
  histoStyle(h1uncFromBayesRecoCov, kRed, 24,3);
  plotHisto(h1uncFromBayesMeasCov,"cCmpUncFromCov",1,1,"LPE", "bayes measCov");
  plotHistoSame(h1uncFromBayesRecoCov,"cCmpUncFromCov","LPE", "bayes recoCov");
  if (h1uncFromUnfCov) {
    histoStyle(h1uncFromUnfCov, kBlack, 5, 2);
    plotHistoSame(h1uncFromUnfCov,"cCmpUncFromCov", "LPE1","my Toys");
    if (0 && h2covRespMStat) {
      TH2D *h2toyCovTot= cloneHisto(h2covRespMStat,"h2toyCovTot","h2toyCovTot");
      h2toyCovTot->Add(h2unfCov);
      TH1D *h1uncFromToysTot= uncFromCov(h2toyCovTot,NULL);
      histoStyle(h1uncFromToysTot, kGreen+1, 27, 3);
      plotHistoSame(h1uncFromToysTot,"cCmpUncFromCov", "PE", "my Toys tot");
    }
  }

  if (0) {
    TH1D* h1measDataErr= errorAsCentral(h1measData);
    TH1D* h1unfDataErr= errorAsCentral(h1unfData);
    histoStyle(h1measDataErr, kGreen, 26, 3);
    histoStyle(h1unfDataErr , 9, 27, 3);
    plotHistoSame(h1measDataErr, "cCmpUncFromCov","PE", "measDataErr");
    plotHistoSame(h1unfDataErr, "cCmpUncFromCov","PE", "unfDataErr");
    printRatio(h1measDataErr, h1uncFromBayesMeasCov);
    printRatio(h1unfDataErr, h1uncFromBayesRecoCov);
    return;
  }
  
  //printRatio(h1uncFromBayesMeasCov,h1uncFromBayesRecoCov);

  plotHisto(h1unfData,"cUnfCheck",1,1,"LPE","unf Data");
  plotHistoSame(h1trueMC,"cUnfCheck","LPE1","true MC");
  printRatio(h1unfData,h1trueMC);

  // chi^2 calculation
  //printRatio(h1trueMC,h1unfData);
  TVectorD vTrueMC( convert2vec(h1trueMC) );
  TVectorD vMeasMC( convert2vec(h1measMC) );
  TVectorD vMeas( convert2vec(h1measData) );
  TVectorD vUnf( convert2vec(h1unfData) );

  double chi2_unfold=0;
  TMatrixD mCovEst (vTrueMC.GetNoElements(),vTrueMC.GetNoElements());
  mCovEst.Zero();
  if (h2unfCov) {
    mCovEst= convert2mat(h2unfCov);
    if (1) {
      TVectorD vDiff_trueMC_unf = vTrueMC - vUnf;
      std::cout << "vTrueMC - vUnf = ";  vDiff_trueMC_unf.Print();
      TMatrixD mCovEstInv(mCovEst);
      mCovEstInv.Invert(NULL);
      TVectorD mCovEstInv_times_vDiffMC = mCovEstInv * vDiff_trueMC_unf;
      double chi2= vDiff_trueMC_unf * mCovEstInv_times_vDiffMC;
      std::cout << "chi2=" << chi2 << "\n";
      chi2_unfold=chi2;
    }

    chi2_unfold= chi2estimate( vTrueMC, vUnf, mCovEst);
    std::cout << "chi2_unfold=" << chi2_unfold << " (using mCovEst -- my estimation)\n";
  }
  else {
    std::cout << "chi2_unfold using mCovEst is not available\n";
  }
  std::cout << "chi2_unfold=" << chi2estimate( vTrueMC, vUnf, mBayesRecoCov ) << " (using bayesRecoCov -- from RooUnfold)\n";

  // chi^2 bottom-line test
  TMatrixD mMeasCovInv(mBayesMeasCov);
  mMeasCovInv.Invert(NULL);
  TMatrixD migMatrix( convert2mat( rooResp->Hresponse() ) );
  if (0) {
    printHisto( rooResp->Hresponse() );
    migMatrix.Print();
  }
  //TMatrixD migMatrixT( TMatrixD::kTransposed, migMatrix );
  TMatrixD mSmear(migMatrix);
  if (0) {
    for (int ir=0; ir<mSmear.GetNrows(); ir++) {
      double sum=0;
      for (int ic=0; ic<mSmear.GetNcols(); ic++) {
	sum += mSmear(ir,ic);
      }
      for (int ic=0; ic<mSmear.GetNcols(); ic++) {
	mSmear(ir,ic) /= sum;
      }
    }
  }
  else { // correct way
    for (int ic=0; ic<mSmear.GetNcols(); ic++) {
      double sum=0;
      for (int ir=0; ir<mSmear.GetNrows(); ir++) {
	sum += mSmear(ir,ic);
      }
      for (int ir=0; ir<mSmear.GetNrows(); ir++) {
	mSmear(ir,ic) /= sum;
      }
    }
  }

  mSmear = rooResp->Mresponse();

  TVectorD vSmearCheck= mSmear * vTrueMC;
  printRatio("vMeasMC", vMeasMC, "vSmearCheck", vSmearCheck);
  if (1) {
    TMatrixD mSmearInv(TMatrixD::kInverted,mSmear);
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

  double chi2_bottomline= chi2estimate( vMeas, vSmearCheck, mBayesMeasCov );
  std::cout << "chi2_bottomline=" << chi2_bottomline << " (using vSmearCheck)\n";
  chi2_bottomline= chi2estimate( vMeas, vMeasMC, mBayesMeasCov );
  std::cout << "chi2_bottomline=" << chi2_bottomline << " (using vMeasMC)\n";

  TDecompSVD svd(mSmear);
  svd.Decompose();

  std::cout << "SVD info:\n";
  std::cout << " - decomposed = " << svd.kDecomposed << "\n";
  std::cout << " - singular = " << svd.kSingular << "\n";
  std::cout << " - condition number= " << svd.Condition() << "\n";

  TVectorD singVal( svd.GetSig() );
  //singVal.Print();
  std::cout << " svd.eigenvalues : " << singVal[0] << " .. " << singVal[singVal.GetNoElements()-1] << "\n";
  double lowBound= singVal[singVal.GetNoElements()-1];
  double condNumb= (lowBound<0) ? singVal[0]/1e-34 : singVal[0]/lowBound;
  std::cout << " condition number sigma_max/max(0,sigma_min)=" << condNumb << "\n";

  if (!doSave) { std::cout << "no saving of canvases\n"; return; }
  std::vector<TCanvas*> cV;
  if (findCanvases("",cV)) {
    TString mainTag="studyFSRResp" + outDirTag;
    mainTag.Append( (nSample<100) ? "_debug" : Form("_%d",nSample) );
    TString destDir="dir-canv-" + mainTag;
    TString outFName="outCanv_" + mainTag + TString(".root");
    TFile *fout= new TFile(outFName,"recreate");
    if (!fout->IsOpen()) {
      std::cout << "failed to create output file <" << fout->GetName() << ">\n";
    }
    else {
      SaveCanvases(cV,destDir,fout);
      TObjString timeTag(DayAndTimeTag(0));
      timeTag.Write("timeTag");
      fout->Close();
      std::cout << "file <" << fout->GetName() << "> created\n";
    }
  }

}
