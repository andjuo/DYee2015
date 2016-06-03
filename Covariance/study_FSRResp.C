// Small macro for testing purposes

#include "inputs.h"
#include "../RooUnfold/src/RooUnfoldResponse.h"
#include "../RooUnfold/src/RooUnfoldBayes.h"
#include "../RooUnfold/src/RooUnfoldInvert.h"
#include "crossSection.h"
//#include <TMatrixD.h>
#include <TDecompSVD.h>

void study_FSRResp(int nSample)
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

  // main result
  TString fnameMain="cs_DYmm_13TeV_76X_cs.root";
  TString nameh1preUnfData="h1PostFSR";
  TString nameh1unfData="h1PreFSR";

  // Response matrix
  TString fname="dymm_test_Mu76X.root";
  TString nameResp="rooUnf_fsrResp";
  TString nameh1meas="h1_postFsr_Mweighted";
  TString nameh1true="h1_preFsr_Mweighted";
  int nIters=10;

  // Covariance estimation
  TString covFName="cov_mumu_varFSRRes_10000.root";
  TString nameh2cov="h2cov_varFSRRes_76X";

  TH1D *h1measData= loadHisto(fnameMain,nameh1preUnfData,"h1measData",h1dummy);
  TH1D *h1unfData= loadHisto(fnameMain,nameh1unfData,"h1unfData",h1dummy);
  if (!h1measData || !h1unfData) {
    std::cout << "failed to load histos from fnameMain=<" << fnameMain << ">\n";
    return;
  }

  RooUnfoldResponse *rooResp= loadRooUnfoldResponse(fname,nameResp,"DRR");
  TH1D* h1measMC= loadHisto(fname,nameh1meas,"h1meas",h1dummy);
  TH1D* h1trueMC= loadHisto(fname,nameh1true,"h1true",h1dummy);
  if (!rooResp || !h1measMC || !h1trueMC) {
    std::cout << "failed to load histos from fname=<" << fname << ">\n";
    return;
  }

  TH2D* h2covRespMStat= loadHisto(covFName, nameh2cov, "h2covRespMStat", h2dummy);
  if (!h2covRespMStat) {
    std::cout << "failed to load histos from covFName=<" << covFName << ">\n";
    return;
  }
  //TH2D* h2corrRespMStat= cov2corr(h2covRespMStat);

  if (0) {
    TMatrixD mCovRespMStat( convert2mat(h2covRespMStat ) );
    TH2D *h2Chk= new TH2D(mCovRespMStat);
    plotCovCorr(h2covRespMStat,"cCovRespMStat",NULL,1);
    plotCovCorr(h2Chk, "cCovChk",NULL,1);
    printHisto(h2Chk);
    mCovRespMStat.Print();
    return;
  }

  // Check MC unfolding on MC
  RooUnfoldBayes fsrBayesMC( rooResp, h1measMC, nIters, false, "fsrBayesMC" );
  TH1D *h1unfBayesMC= cloneHisto(fsrBayesMC.Hreco(),"h1UnfMC",
				 TString("UnfMC;") +
				 h1trueMC->GetXaxis()->GetTitle() + TString(";")
				 + "unfolded", h1trueMC);

  histoStyle(h1measMC, kBlue,24);
  histoStyle(h1trueMC, kRed,4);

  plotHisto(h1measMC, "cCmpMC", 1,1, "LPE1");
  plotHistoSame(h1trueMC, "cCmpMC", "LPE");
  plotHistoSame(h1unfBayesMC, "cCmpMC", "LPE");
  printRatio(h1trueMC,h1unfBayesMC);


  // Compare data unfolding
  RooUnfoldBayes fsrBayesDataMy( rooResp, h1measData, nIters, false, "fsrBayesDataMy" );
  TMatrixD bayesMeasCovAll(fsrBayesDataMy.GetMeasuredCov());
  TMatrixD bayesMeasCov(submatrix(bayesMeasCovAll,1,bayesMeasCovAll.GetNrows()-1));
  TH2D *h2bayesMeasCov= new TH2D(bayesMeasCov);
  h2bayesMeasCov->SetDirectory(0);
  h2bayesMeasCov->SetBinContent(0,0, 0.);
  TMatrixD bayesRecoCovAll(fsrBayesDataMy.Ereco(RooUnfoldBayes::kCovariance));
  TMatrixD bayesRecoCov(submatrix(bayesRecoCovAll,1,bayesRecoCovAll.GetNrows()-1));
  TH2D *h2bayesRecoCov= new TH2D(bayesRecoCov);
  h2bayesRecoCov->SetDirectory(0);
  h2bayesRecoCov->SetBinContent(0,0, 0.);

  TH1D *h1unfBayesMy= cloneHisto(fsrBayesDataMy.Hreco(), "h1unfBayesMy",
				 TString("UnfDataMy;") +
				 h1unfData->GetXaxis()->GetTitle() +
				 TString("; unfolded(my)"), h1unfData);

  histoStyle(h1unfBayesMy, kBlue, 25);

  plotHisto(h1unfBayesMy, "cCmpData", 1,1, "LPE1");
  plotHistoSame(h1unfData, "cCmpData", "LPE");
  plotHistoSame(h1measData, "cCmpData", "LPE");
  printRatio(h1unfBayesMy, h1unfData);
  //printHisto(h1unfBayesMy);
  //printHisto(h1unfData);

  std::vector<TH1D*> vecPreFSR;
  int nonNegative=0;
  for (int i=0; i<nSample; i++) {
    TString nameH1Var= TString("hVar") + Form("%d",i);
    TH1D *h1_in= cloneHisto(h1measData,nameH1Var,nameH1Var);
    randomizeWithinErr(h1measData, h1_in, nonNegative);
    RooUnfoldBayes *fsrBayes= new RooUnfoldBayes( rooResp, h1_in, nIters, false, "bayes_"+TString(Form("%d",i)));
    TString outTag=TString("h1_out_") + Form("%d",i);
    TH1D *h1_out= cloneHisto(fsrBayes->Hreco(),outTag,outTag, h1unfData);
    vecPreFSR.push_back(h1_out);
    delete fsrBayes;
    delete h1_in;
  }
  TH1D *h1avgUnf=NULL;
  TH2D *h2unfCov=NULL;
  deriveCovariance(vecPreFSR,"h2cov","h2covUnf",&h1avgUnf,&h2unfCov);

  if (1 && (vecPreFSR.size()>0)) {
    plotHisto(vecPreFSR[0],"cToy",1,1,"LPE1");
    for (unsigned int i=1; (i<50) && (i<vecPreFSR.size()); i++) {
      plotHistoSame(vecPreFSR[i],"cToy","LPE");
    }
  }

  // compare covariances
  plotCovCorr(h2covRespMStat,"cCovRespMStat",NULL,1);
  plotCovCorr(h2bayesMeasCov,"cCovBayesMeas",NULL,1);
  plotCovCorr(h2bayesRecoCov,"cCovBayesReco",NULL,1);
  plotCovCorr(h2unfCov,"cCovUnf",NULL,1);

  // compare uncertainties from the covariance
  TH1D* h1uncFromUnfCov= uncFromCov(h2unfCov,NULL);
  TH1D* h1uncFromBayesMeasCov= cloneHisto(h1uncFromUnfCov,"h1uncFromBayesMeasCov","h1uncFromBayesMeasCov");
  copyContents(h1uncFromBayesMeasCov, uncFromCov(h2bayesMeasCov,NULL));
  TH1D* h1uncFromBayesRecoCov= cloneHisto(h1uncFromUnfCov,"h1uncFromBayesRecoCov","h1uncFromBayesRecoCov");
  copyContents(h1uncFromBayesRecoCov, uncFromCov(h2bayesRecoCov,NULL));
  histoStyle(h1uncFromBayesMeasCov, kBlue, kPlus);
  histoStyle(h1uncFromUnfCov, kBlack, 5, 2);
  histoStyle(h1uncFromBayesRecoCov, kRed, 24,3);
  plotHisto(h1uncFromBayesMeasCov,"cCmpUncFromCov",1,1,"LPE", "bayes measCov");
  plotHistoSame(h1uncFromBayesRecoCov,"cCmpUncFromCov","LPE", "bayes recoCov");
  plotHistoSame(h1uncFromUnfCov,"cCmpUncFromCov", "LPE1","my Toys");

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

  // chi^2 calculation
  //printRatio(h1trueMC,h1unfData);
  TVectorD vTrueMC( convert2vec(h1trueMC) );
  TVectorD vMeasMC( convert2vec(h1measMC) );
  TVectorD vMeas( convert2vec(h1measData) );
  TVectorD vUnf( convert2vec(h1unfData) );

  double chi2_unfold=0;
  TMatrixD mCovEst( convert2mat(h2unfCov) );
  if (0) {
    TVectorD vDiff_trueMC_unf = vTrueMC - vUnf;
    vDiff_trueMC_unf.Print();
    TMatrixD mCovEstInv(mCovEst);
    mCovEstInv.Invert(NULL);
    TVectorD mCovEstInv_times_vDiffMC = mCovEstInv * vDiff_trueMC_unf;
    double chi2= vDiff_trueMC_unf * mCovEstInv_times_vDiffMC;
    std::cout << "chi2=" << chi2 << "\n";
    chi2_unfold=chi2;
  }

  chi2_unfold= chi2estimate( vTrueMC, vUnf, mCovEst);
  std::cout << "chi2_unfold=" << chi2_unfold << " (using mCovEst -- my estimation)\n";
  std::cout << "chi2_unfold=" << chi2estimate( vTrueMC, vUnf, bayesRecoCov ) << " (using bayesRecoCov -- from RooUnfold)\n";

  // chi^2 bottom-line test
  TMatrixD mMeasCovInv(bayesMeasCov);
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

  TVectorD vSmearCheck= mSmear * vTrueMC;
  printRatio("vMeasMC", vMeasMC, "vSmearCheck", vSmearCheck);

  double chi2_bottomline= chi2estimate( vMeas, vSmearCheck, bayesMeasCov );
  std::cout << "chi2_bottomline=" << chi2_bottomline << " (using vSmearCheck)\n";
  chi2_bottomline= chi2estimate( vMeas, vMeasMC, bayesMeasCov );
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

}
