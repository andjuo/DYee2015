#include "crossSection.h"

// --------------------------------------------------------------

void work(TVersion_t inpVer,
	  MuonCrossSection_t &muCS, TVaried_t var, int nSample, int doSave,
	  int nMaxSamples=-1);

// --------------------------------------------------------------

void studyDYmmCS(TVaried_t var= _varNone, int nSample=10, int doSave=0,
		 int nMaxSamples=-1)
{
  TVersion_t inpVer=_verMu1;
  TString inpVerTag="_v1";
  if (0) { inpVer=_verMu76X; inpVerTag="_76X"; }
  else if (0) { inpVer=_verMuApproved; inpVerTag=versionName(inpVer); }
  else if (1) { inpVer=_verMuMay2017; inpVerTag=versionName(inpVer); }

  MuonCrossSection_t muCS("muCS",inpVerTag,0,0,inpVer);
  if (!muCS.load("cs_DYmm_13TeV" + inpVerTag, inpVerTag)) {
    std::cout << "loading failed\n";
    return;
  }

  //TCanvas *ca=muCS.editCSa().plotCrossSection("csa");
  //TCanvas *cb=muCS.editCSb().plotCrossSection("csb");

  if (var==_varNone) {
    TCanvas *c= muCS.plotCrossSection("cs");
    if (!c) return;

    plotHisto(cloneHisto(muCS.h1Bkg(),"h1tmp1","tmp1"),"ctmp1",1,1,"hist",
	      "original bkg");

    std::vector<double> bkgW;
    bkgW.reserve(muCS.bkgWUnc().GetNoElements());
    for (int i=0; i<muCS.bkgWUnc().GetNoElements(); i++) {
      bkgW.push_back(1.);
    }
    if (!muCS.recalcBkg(&bkgW)) return;
    plotHistoSame(cloneHisto(muCS.h1Bkg(),"h1tmp2","tmp2"),"ctmp1","LPE",
		  "modified bkg");
 
    TCanvas *c2= muCS.plotCrossSection("cs2",1);
    if (!c2) return;

    printHisto(muCS.h1Bkg());
    printRatio(muCS.h1CS(),muCS.h1Theory());
  }
  else {
    std::cout << "perform study\n";
    work(inpVer,muCS,var,nSample,doSave,nMaxSamples);
  }

  std::cout << "macro ran with inpVerTag=" << inpVerTag << "\n";
}


// --------------------------------------------------------------
// --------------------------------------------------------------

void work(TVersion_t inpVer,
	  MuonCrossSection_t &muCS, TVaried_t var, int nSample, int doSave,
	  int nMaxSamples)
{
  std::vector<TH1D*> rndCSVec, rndCSaVec, rndCSbVec;
  int res=1;
  int check_KLRndEffMap=(1 && (nSample==500) && (nMaxSamples==-500)) ? 1:0;
  int removeNegativeSignal=1;
  if (var!=_varRhoFile) {
    std::cout << "nSample=" << nSample << ", removeNegativeSignal=" << removeNegativeSignal << "\n";
    res=muCS.sampleRndVec(var,nSample,removeNegativeSignal,
			  rndCSVec,&rndCSaVec,&rndCSbVec);
  }
  else {
    TString inpVerTag=versionName(inpVer);
    TString loadFName=Form("dir-Rho%s/dymm_rhoRndVec_%s_%d.root",
			   inpVerTag.Data(),inpVerTag.Data(),
			   nSample);
    // hack to verify KL rnd maps
    if (check_KLRndEffMap) {
      loadFName.ReplaceAll(".root","_KL.root");
    }
    RndVecInfo_t info(loadFName,"h1rho_4p2_var","h1rho_4p3_var",nMaxSamples);
    if ((nMaxSamples>0) && (nMaxSamples<nSample)) nSample=nMaxSamples;
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
    else if (doSave==3) fname.ReplaceAll(".root","_veryslim.root");
    if (var==_varFSRRes_Poisson) fname.ReplaceAll(".root","-Poisson.root");
    if (check_KLRndEffMap) fname.ReplaceAll(".root","_KL.root");
    TFile fout(fname,"RECREATE");
    if (!fout.IsOpen()) {
      std::cout << "file <" << fname << "> could not be created\n";
      return;
    }
    std::vector<TString> dirs;
    std::vector<std::vector<TH1D*>*> hVs;
    if (doSave!=3) {
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
    if (doSave!=3) {
      std::vector<TCanvas*> canvV;
      TString canvList="muCS_a muCS_b";
      canvList.Append("cVaried"+muCS.csA().tag());
      canvList.Append("cVaried"+muCS.csB().tag());
      if (doSave==2) canvList="";
      canvList.Append(" cVaried_" + variedVarName(var) + "_" + muCS.csA().tag());
      canvList.Append(" cVaried_" + variedVarName(var) + "_" + muCS.csB().tag());
      std::cout << "canvList=" << canvList << "\n";
      if (findCanvases(canvList,canvV)) {
	for (unsigned int ic=0; ic<canvV.size(); ic++) {
	  canvV[ic]->Write();
	}
      }
    }
    TObjString timeTag(DayAndTimeTag(0));
    timeTag.Write("timeTag");
    fout.Close();
    std::cout << "file <" << fout.GetName() << "> created\n";
  }
}

// --------------------------------------------------------------


