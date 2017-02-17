#include "crossSection.h"
//#include "studyCSSteps.h"

int distributionByXBin(const std::vector<TH1D*> &h1V,
		       const TH1D *h1base, int nSigma,
		       TString canvNameBase,
		       const TH1D *h1lowEdge=NULL,
		       int subtract=0,
		       double extraScale=1.);

// --------------------------------------------------------------
// --------------------------------------------------------------

void studyDYmmCScovYield(int nSample=10, int method=1,
			 int doSave=0, int plotSetA=-1)
{
  TVaried_t var=_varYield;
  closeCanvases(10);

  std::cout << "doSave=" << doSave << "\n";

  TVersion_t inpVer=_verMu1;
  TString inpVerTag="_v1";
  if (0) { inpVer=_verMu76X; inpVerTag="_76X"; }
  else if (1) { inpVer=_verMuApproved; inpVerTag=versionName(inpVer); }

  MuonCrossSection_t muCS("muCS",inpVerTag,0,0,inpVer);
  if (!muCS.load("cs_DYmm_13TeV" + inpVerTag, inpVerTag)) {
    std::cout << "loading failed\n";
    return;
  }

  // preserve original cross section
  TH1D *h1csOrig= cloneHisto(muCS.h1CS(),"h1csOrig","h1csOrig");
  //muCS.h1Theory(muCS.h1CS());
  //printHisto(muCS.h1Theory()); return;

  TString covFName="cov_mumu_" + variedVarName(var) + Form("_%d.root",nSample);
  if ((nSample==500) || (nSample==5000)) covFName.ReplaceAll(".root","_slim.root");
  TString cNameA="cVaried" + versionName(inpVer) + "a";
  TString cNameB="cVaried" + versionName(inpVer) + "b";
  if (covFName.Index("_slim")!=-1) {
    cNameA.ReplaceAll("cVaried","cVaried_"+variedVarName(var)+"_");
    cNameB.ReplaceAll("cVaried","cVaried_"+variedVarName(var)+"_");
  }
  TFile fin(covFName);
  TCanvas *canvA=(TCanvas*)fin.Get(cNameA);
  TCanvas *canvB=(TCanvas*)fin.Get(cNameB);
  fin.Close();
  if (!canvA || !canvB) {
    std::cout << "failed to get canvases\n";
    return;
  }
  canvA->Draw();
  canvB->Draw();

  std::vector<TH1D*> h1aV_all, h1aV;
  std::vector<TH1D*> h1bV_all, h1bV;
  getHistosFromCanvas(canvA,&h1aV_all,NULL);
  getHistosFromCanvas(canvB,&h1bV_all,NULL);
  std::cout << "loaded histograms: " << h1aV_all.size()
	    << " from CSa and " << h1bV_all.size() << " from CSb\n";
  if ((h1aV_all.size()==0) || (h1bV_all.size()==0)) {
    std::cout << "no histos\n";
    return;
  }

  int keepOnlyPositiveSignal=(method==3) ? 1 : 0;
  const TH1D* h1bkgA= muCS.csA().h1Bkg();
  const TH1D* h1bkgB= muCS.csB().h1Bkg();

  //std::cout << "h1YieldA last bin: "; printHisto(muCS.csA().h1Yield());
  if (1) {
    // Correct muon bkg last bin in csA
    std::cout << "\n\tcorrecting background in the last mass bin\n";
    muCS.editCSa().editH1Bkg()->SetBinContent(43,0);
    std::cout << "\n";
  }
  std::cout << "h1bkgA last bin: " << h1bkgA->GetBinContent(h1bkgA->GetNbinsX()) << "\n";

  h1aV.reserve(h1aV_all.size());
  for (unsigned int i=0; i<h1aV_all.size(); i++) {
    if (TString(h1aV_all[i]->GetName()).Index("h1varTmp")!=-1) {
      // debug printout
      if (0 && !allGE(h1aV_all[i],h1bkgA)) {
	TH1D *h1= cloneHisto(h1aV_all[i],h1aV_all[i]->GetName() + TString("_tt"),
			     "");
	h1->Add(h1bkgA,-1);
	for (int ibin=1; ibin<=h1->GetNbinsX(); ibin++) {
	  if (h1->GetBinContent(ibin)<0) {
	    std::cout << "histoA " << i << " is negative in "
		      << h1->GetBinLowEdge(ibin) << " .. "
		      << (h1->GetBinLowEdge(ibin)+h1->GetBinWidth(ibin)) << "\n";
	    break;
	  }
	}
      }
      // debug printout (end)

      if (!keepOnlyPositiveSignal ||
	  (keepOnlyPositiveSignal && allGE(h1aV_all[i],h1bkgA))) {
	//std::cout << "keep " << h1aV_all[i]->GetName() << "\n";
	h1aV.push_back(h1aV_all[i]);
      }
    }
  }
  h1bV.reserve(h1bV_all.size());
  for (unsigned int i=0; i<h1bV_all.size(); i++) {
    if (TString(h1bV_all[i]->GetName()).Index("h1varTmp")!=-1) {
      if (!keepOnlyPositiveSignal ||
	  (keepOnlyPositiveSignal && allGE(h1bV_all[i],h1bkgB))) {
	//std::cout << "keep " << h1bV_all[i]->GetName() << "\n";
	h1bV.push_back(h1bV_all[i]);
      }
    }
  }
  std::cout << "kept " << h1aV.size() << " and " << h1bV.size() << " histograms\n";
  if ((h1aV.size()==0) || (h1bV.size()==0)) {
    std::cout << "no histos\n";
    return;
  }

  TH1D *h1var= NULL;
  TH1D *h1bkg= NULL;
  if (plotSetA) {
    h1var=muCS.editCSa().getVariedHisto(var);
    h1bkg= cloneHisto(muCS.editCSa().h1Bkg(),"h1bkg_clone","h1bkg_clone");
  }
  else {
    h1var=muCS.editCSb().getVariedHisto(var);
    h1bkg= cloneHisto(muCS.editCSb().h1Bkg(),"h1bkg_clone","h1bkg_clone");
  }

  if (0) {
    plotHisto(h1var,"cVar",1,1,"LPE1",Form("orig yield (csA=%d)",plotSetA));
    plotHistoSame(h1aV[3],"cVar","3");
    return;
  }

  if (0) {
    int nSigma=10;
    std::vector<TH1D*> h1RndDistrInBinsV;
    for (int ibin=1; ibin<=h1var->GetNbinsX(); ibin++) {
      double m1= h1var->GetBinLowEdge(ibin);
      double m2= m1 + h1var->GetBinWidth(ibin);
      TString massRange=Form("%2.0lf-%2.0lf",m1,m2);
      double count= h1var->GetBinContent(ibin);
      double sigma= h1var->GetBinError(ibin);
      TH1D *h1= new TH1D("h1RndDistr_" + massRange,
		 "h1RndDistr_" + massRange + ";event count;presence count",
		 2*nSigma*sigma+1,count-nSigma*sigma,count+nSigma*sigma);
      h1->SetDirectory(0);
      for (unsigned int i=0; i<h1aV.size(); i++) {
	h1->Fill(h1aV[i]->GetBinContent(ibin));
      }
      h1->GetYaxis()->SetTitleOffset(1.5);
      TCanvas *cm=plotHisto(h1,"cMassRange_"+massRange,0,0,"LPE1");
      //cm->SetLeftMargin(0.2);
      TLine *line= new TLine(h1bkg->GetBinContent(ibin),
			     h1->GetMinimum(),
			     h1bkg->GetBinContent(ibin),
			     h1->GetMaximum());
      line->SetLineStyle(2);
      line->SetLineColor(kRed);
      line->Draw();
      cm->Modified();
      cm->Update();
    }
  }
  else if (plotSetA>=0) {
    int nSigma=10;
    int subtract=0;
    if (plotSetA==1) {
      distributionByXBin(h1aV,h1var,nSigma,"cMassRangeA_",
			 h1bkgA,subtract);
    }
    if (plotSetA==0) {
      distributionByXBin(h1bV,h1var,nSigma,"cMassRangeB_",
			 h1bkgB,subtract);
    }
  }

  if (method) std::cout << "method=" << method << "\n";

  int removeNegativeSignal= (method==2) ? 1:0;

  unsigned int imax=h1aV.size();
  if (imax>h1bV.size()) imax=h1bV.size();
  //imax=100;
  std::vector<TH1D*> h1csV;
#ifdef studyCSSteps_H
  StudyCSSteps_t csSteps;
#endif

  for (unsigned int i=0; i<imax; i++) {
    TH1D *h1cs= muCS.calcCrossSection(var,i,h1aV[i],h1bV[i],removeNegativeSignal);
    if (!h1cs) {
      std::cout << "got null h1cs\n";
      return;
    }
    copyStyle(h1cs,h1aV[i]);
    h1csV.push_back(h1cs);
#ifdef studyCSSteps_H
    if (!csSteps.addInfo(muCS,Form("_study%d",i),h1aV[i])) {
      std::cout << "csSteps.addInfo error\n";
      return;
    }
#endif
  }

  TString calcTag= "_allRnd";
  if (keepOnlyPositiveSignal) calcTag="_keepYieldWithPosSignal";
  if (removeNegativeSignal) calcTag="_nulifyNegSig";

  if (1) {
    TH1D *h1csForPlot= cloneHisto(h1csOrig,"h1csForPlot","h1csForPlot");
    for (int ibin=1; ibin<h1csForPlot->GetNbinsX(); ibin++) {
      h1csForPlot->SetBinError(ibin, 0.025*h1csForPlot->GetBinContent(ibin));
    }
    distributionByXBin(h1csV,h1csForPlot,100,"cCS"+calcTag,NULL,0,1e4);
  }

  if (1) {
    TString canvName="cs";
    TCanvas *cx=muCS.plotCrossSection(canvName,0,0);
    if (1)
    for (unsigned int i=0; i<h1csV.size(); i++) {
      plotHistoSame(h1csV[i],canvName,"LP","");
      if (i>=200) { std::cout << "stopping at " << i << "\n"; break; }
    }
    if(!cx) std::cout << "null canvas\n";
  }

#ifdef studyCSSteps_H
  {
    TH1D *h1YieldA_avg=NULL, *h1SigA_avg=NULL, *h1UnfA_avg=NULL;
    TH1D *h1UnfRhoCorrA_avg=NULL, *h1UnfRhoEffCorrA_avg=NULL;
    TH1D *h1UnfRhoEffAccCorrA_avg=NULL;
    TH1D *h1PostFsrA_avg=NULL, *h1PreFsrA_avg=NULL, *h1CSA_avg=NULL;
    TH1D *h1YieldB_avg=NULL, *h1SigB_avg=NULL, *h1UnfB_avg=NULL;
    TH1D *h1UnfRhoCorrB_avg=NULL, *h1UnfRhoEffCorrB_avg=NULL;
    TH1D *h1UnfRhoEffAccCorrB_avg=NULL;
    TH1D *h1PostFsrB_avg=NULL, *h1PreFsrB_avg=NULL, *h1CSB_avg=NULL;
    TH1D *h1PostFsr_avg=NULL, *h1PreFsr_avg=NULL, *h1CS_avg=NULL;
    deriveCovariance(csSteps.csA()->h1YieldV(),"_Yield4p2","yield 4p2",&h1YieldA_avg,NULL);
    deriveCovariance(csSteps.csA()->h1SigV(),"_Sig4p2","signal 4p2",&h1SigA_avg,NULL);
    deriveCovariance(csSteps.csA()->h1UnfV(),"_Unf4p2","unf 4p2",&h1UnfA_avg,NULL);
    deriveCovariance(csSteps.csA()->h1UnfRhoCorrV(),"_UnfRhoCorr4p2","unf/rho 4p2",&h1UnfRhoCorrA_avg,NULL);
    deriveCovariance(csSteps.csA()->h1UnfRhoEffCorrV(),"_UnfRhoEffCorr4p2","unf/rho/eff 4p2",&h1UnfRhoEffCorrA_avg,NULL);
    deriveCovariance(csSteps.csA()->h1UnfRhoEffAccCorrV(),"_UnfRhoEffCorrAcc4p2","unf/rho/eff/A 4p2",&h1UnfRhoEffAccCorrA_avg,NULL);
    deriveCovariance(csSteps.csA()->h1PostFsrV(),"_postFSR4p2","postFSR 4p2",&h1PostFsrA_avg,NULL);
    deriveCovariance(csSteps.csA()->h1PreFsrV(),"_preFSR4p2","preFSR 4p2",&h1PreFsrA_avg,NULL);
    deriveCovariance(csSteps.csA()->h1CSV(),"_cs4p2","cs 4p2",&h1CSA_avg,NULL);

    deriveCovariance(csSteps.csB()->h1YieldV(),"_Yield4p3","yield 4p3",&h1YieldB_avg,NULL);
    deriveCovariance(csSteps.csB()->h1SigV(),"_Sig4p3","signal 4p3",&h1SigB_avg,NULL);
    deriveCovariance(csSteps.csB()->h1UnfV(),"_Unf4p3","unf 4p3",&h1UnfB_avg,NULL);
    deriveCovariance(csSteps.csB()->h1UnfRhoCorrV(),"_UnfRhoCorr4p3","unf/rho 4p3",&h1UnfRhoCorrB_avg,NULL);
    deriveCovariance(csSteps.csB()->h1UnfRhoEffCorrV(),"_UnfRhoEffCorr4p3","unf/rho/eff 4p3",&h1UnfRhoEffCorrB_avg,NULL);
    deriveCovariance(csSteps.csB()->h1UnfRhoEffAccCorrV(),"_UnfRhoEffCorrAcc4p3","unf/rho/eff/A 4p3",&h1UnfRhoEffAccCorrB_avg,NULL);
    deriveCovariance(csSteps.csB()->h1PostFsrV(),"_postFSR4p3","postFSR 4p3",&h1PostFsrB_avg,NULL);
    deriveCovariance(csSteps.csB()->h1PreFsrV(),"_preFSR4p3","preFSR 4p3",&h1PreFsrB_avg,NULL);
    deriveCovariance(csSteps.csB()->h1CSV(),"_cs4p3","cs 4p3",&h1CSB_avg,NULL);

    deriveCovariance(csSteps.h1PostFsrV(),"_postFSR","postFSR",&h1PostFsr_avg,NULL);
    deriveCovariance(csSteps.h1PreFsrV(),"_preFSR","preFSR",&h1PreFsr_avg,NULL);
    deriveCovariance(csSteps.h1CSV(),"_cs","cs",&h1CS_avg,NULL);

    TFile fout("dymm_study_yieldErrPropagation.root","recreate");
    h1YieldA_avg->Write();
    h1SigA_avg->Write();
    h1UnfA_avg->Write();
    h1UnfRhoCorrA_avg->Write();
    //h1UnfRhoEffCorrA_avg->Write();
    h1UnfRhoEffAccCorrA_avg->Write();
    h1PostFsrA_avg->Write();
    //h1PreFsrA_avg->Write();
    //h1CSA_avg->Write();
    h1YieldB_avg->Write();
    h1SigB_avg->Write();
    h1UnfB_avg->Write();
    h1UnfRhoCorrB_avg->Write();
    //h1UnfRhoEffCorrB_avg->Write();
    h1UnfRhoEffAccCorrB_avg->Write();
    h1PostFsrB_avg->Write();
    //h1PreFsrB_avg->Write();
    //h1CSB_avg->Write();
    h1PostFsr_avg->Write();
    h1PreFsr_avg->Write();
    h1CS_avg->Write();
    writeTimeTag(&fout);
    fout.Close();
  }
#endif

  TH1D *h1avgCS=NULL;
  TH2D *h2cov=NULL, *h2corr=NULL;
  PlotCovCorrOpt_t ccOpt;
  deriveCovariance(h1csV,calcTag,"csCov_"+variedVarName(var),&h1avgCS,&h2cov);
  plotCovCorr(h2cov,"csCov",ccOpt,&h2corr);


  TString outFName=covFName;
  outFName.ReplaceAll(".root",calcTag + ".root");
  std::cout << "outFName=" << outFName << "\n";
  if (doSave) {
    std::cout << "should be saving\n";
    TFile fout(outFName,"RECREATE");
    if (!fout.IsOpen()) {
      std::cout << "failed to open the file <" << fout.GetName() << ">\n";
      return;
    }
    h2cov->Write();
    //h2cov->Write("h2Cov"); // needed for combination code
    h2corr->Write();
    SaveCanvases("ALL",
		 "dir-dymmCScov"+variedVarName(var)+Form("-%d",nSample)+calcTag,
		 &fout);
    writeTimeTag(&fout);
    fout.Close();
  }

  std::cout << "macro ran with inpVerTag=" << inpVerTag << "\n";
}


// --------------------------------------------------------------
// --------------------------------------------------------------

int distributionByXBin(const std::vector<TH1D*> &h1V,
		       const TH1D *h1base, int nSigma,
		       TString canvNameBase,
		       const TH1D *h1lowEdge,
		       int subtract,
		       double extraScale)
{
  if (subtract && !h1lowEdge) {
    std::cout << "distributionByXBin: subtract=1, but h1lowEdge is null\n";
    return 0;
  }
  TH1D *h1var=NULL;
  h1var= cloneHisto(h1base,h1base->GetName()+("_"+canvNameBase),
		    h1base->GetTitle());
  if (!h1base) {
    std::cout << "h1base is null\n";
    return 0;
  }
  if (subtract) h1var->Add(h1lowEdge,-1);

  std::vector<TH1D*> h1RndDistrInBinsV;
  for (int ibin=1; ibin<=h1var->GetNbinsX(); ibin++) {
    double m1= h1var->GetBinLowEdge(ibin);
    double m2= m1 + h1var->GetBinWidth(ibin);
    TString massRange=Form("%2.0lf-%2.0lf",m1,m2);
    double count= extraScale*h1var->GetBinContent(ibin);
    double sigma= extraScale*h1var->GetBinError(ibin);
    int nBins=2*nSigma+1;
    if (nBins>100) nBins=100;
    if (0) {
      std::cout << "ibin=" << ibin << ", count=" << count
		<< " , sigma=" << sigma << ", nBins=" << nBins << "\n";
    }
    TH1D *h1= new TH1D("h1RndDistr_" + massRange,
	       "h1RndDistr_" + massRange + ";event count;presence count",
	       nBins,count-nSigma*sigma,count+nSigma*sigma);
    if (!h1) {
      std::cout << "could not create h1RndDistr_\n";
      return 0;
    }
    h1->SetDirectory(0);
    for (unsigned int i=0; i<h1V.size(); i++) {
      double cnt=h1V[i]->GetBinContent(ibin);
      if (subtract) cnt -= h1lowEdge->GetBinContent(ibin);
      h1->Fill(extraScale*cnt);
    }
    h1->GetYaxis()->SetTitleOffset(1.5);
    TCanvas *cm=plotHisto(h1,canvNameBase+"_"+massRange,0,0,"LPE1");
    //cm->SetLeftMargin(0.2);
    if (h1lowEdge) {
      double x= extraScale*h1lowEdge->GetBinContent(ibin);
      if (subtract) x=0;
      TLine *line= new TLine(x,h1->GetMinimum(),x,h1->GetMaximum());
      line->SetLineStyle(2);
      line->SetLineColor(kRed);
      line->Draw();
      cm->Modified();
      cm->Update();
    }
  }
  return 1;
}

// --------------------------------------------------------------
// --------------------------------------------------------------

void work(TVersion_t inpVer,
	  MuonCrossSection_t &muCS, TVaried_t var, int nSample, int doSave)
{
  std::vector<TH1D*> rndCSVec, rndCSaVec, rndCSbVec;
  int res=1;
  int removeNegativeSignal=1;
  if (var!=_varRhoFile)
    res=muCS.sampleRndVec(var,nSample,removeNegativeSignal,
			  rndCSVec,&rndCSaVec,&rndCSbVec);
  else {
    TString inpVerTag=versionName(inpVer);
    TString loadFName=Form("dir-Rho%s/dymm_rhoRndVec_%s_%d.root",
			   inpVerTag.Data(),inpVerTag.Data(),
			   nSample);
    RndVecInfo_t info(loadFName,"h1rho_4p2_var","h1rho_4p3_var");
    res=muCS.sampleRndVec(var,nSample,info,removeNegativeSignal,
			  rndCSVec,&rndCSaVec,&rndCSbVec);
  }
  if (!res) return;
  TH1D* h1Avg=NULL;
  TH2D* h2Cov=NULL;
  if (!muCS.deriveCov(rndCSVec,var,&h1Avg,&h2Cov)) return;
  
  TCanvas *c= muCS.plotCrossSection("cs");
  if (!c) return;
  h1Avg->SetLineColor(kGreen+1);
  h1Avg->SetMarkerColor(kGreen+1);
  plotHistoSame(h1Avg,"cs","LPE","rnd avg");
  printRatio(muCS.h1CS(), h1Avg);

  TH2D* h2Corr=NULL;
  PlotCovCorrOpt_t ccOpt;
  TCanvas *cx= plotCovCorr(h2Cov,"ccov",ccOpt,&h2Corr);
  if (!cx) std::cout << "cx is null\n";

  TH1D *h1uncFromCov= uncFromCov(h2Cov);
  TH1D *h1relUncFromCov= uncFromCov(h2Cov,muCS.h1CS());

  if (!doSave) {
    std::cout << "Comparing uncertainty from covariance to the uncertainty "
	      << "in the provided cross-section\n";
    std::cout << "Expectation is that the covariance uncertainties will not "
	      << "exceed those from the cross-section\n";
    TH1D *h1CSUnc= errorAsCentral(muCS.h1CS(),0);
    TH1D *h1CSRelUnc= errorAsCentral(muCS.h1CS(),1);
    histoStyle(h1CSUnc,kBlue,24);
    histoStyle(h1CSRelUnc,kBlue,24);
    plotHisto(h1uncFromCov,"cUncCmp",1,1,"hist","uncFromCov");
    plotHistoSame(h1CSUnc,"cUncCmp","P","csUnc");

    plotHisto(h1relUncFromCov,"cRelUncCmp",1,1,"hist","rel.unc.from.cov");
    plotHistoSame(h1CSRelUnc,"cRelUncCmp","P","cs.rel.unc");
  }

  if (doSave) {
    TString fname="cov_mumu_" + variedVarName(var) + Form("_%d.root",nSample);
    if (doSave==2) fname.ReplaceAll(".root","_slim.root");
    if (var==_varFSRRes_Poisson) fname.ReplaceAll(".root","-Poisson.root");
    TFile fout(fname,"RECREATE");
    if (!fout.IsOpen()) {
      std::cout << "file <" << fname << "> could not be created\n";
      return;
    }
    std::vector<TString> dirs;
    std::vector<std::vector<TH1D*>*> hVs;
    if (doSave!=2) {
      dirs.push_back("rnd_CSa");
      hVs.push_back(&rndCSaVec);
      dirs.push_back("rnd_CSb");
      hVs.push_back(&rndCSbVec);
      dirs.push_back("rnd_CStot");
      hVs.push_back(&rndCSVec);
    }
    for (unsigned int iDir=0; iDir<dirs.size(); iDir++) {
      std::vector<TH1D*> *hV= hVs[iDir];
      if (hV->size()) {
	fout.mkdir(dirs[iDir]);
	fout.cd(dirs[iDir]);
	for (unsigned int ih=0; ih<hV->size(); ih++) {
	  hV->at(ih)->Write();
	}
	fout.cd("");
      }
    }
    if (muCS.h1CS()) muCS.h1CS()->Write("xsec");
    if (muCS.h1Theory()) muCS.h1Theory()->Write("xsec_KL");
    if (h1Avg) h1Avg->Write("h1xsecAvg");
    if (h2Cov) h2Cov->Write("h2Cov");
    if (h2Corr) h2Corr->Write("h2Corr");
    if (h1uncFromCov) h1uncFromCov->Write("h1uncFromCov");
    if (h1relUncFromCov) h1relUncFromCov->Write("h1relUncFromCov");
    if (c) c->Write();
    if (cx) cx->Write();
    std::vector<TCanvas*> canvV;
    TString canvList="muCS_a muCS_b";
    canvList.Append(" cVaried" + muCS.csA().tag());
    canvList.Append(" cVaried" + muCS.csB().tag());
    if (doSave==2) canvList="";
    std::cout << "canvList=" << canvList << "\n";
    if (findCanvases(canvList,canvV)) {
      for (unsigned int ic=0; ic<canvV.size(); ic++) {
	canvV[ic]->Write();
      }
    }
    TObjString timeTag(DayAndTimeTag(0));
    timeTag.Write("timeTag");
    fout.Close();
    std::cout << "file <" << fout.GetName() << "> created\n";
  }
}

// --------------------------------------------------------------


