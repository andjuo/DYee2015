#include "crossSection.h"

// --------------------------------------------------------------

void work(TVersion_t inpVer,
	  CrossSection_t &muCS, TVaried_t var, int nSample, int doSave);

// --------------------------------------------------------------

void studyDYeeCS(TVaried_t var= _varNone, int nSample=10, int doSave=0)
{
  TVersion_t inpVer=_verEl3;
  TString inpVerTag=versionName(inpVer);

  CrossSection_t eeCS("elCS",inpVerTag,_csPreFsrFullSp,inpVer);
  if (!eeCS.load("cs_DYee_13TeV_" + inpVerTag + ".root", inpVerTag)) {
    std::cout << "loading failed\n";
    return;
  }

  if (var==_varNone) {
    TCanvas *c= eeCS.plotCrossSection("cs");
    if (!c) return;
  }
  else {
    std::cout << "perform study\n";
    work(inpVer,eeCS,var,nSample,doSave);
  }

  std::cout << "macro ran with inpVerTag=" << inpVerTag << "\n";
}


// --------------------------------------------------------------
// --------------------------------------------------------------

void work(TVersion_t inpVer,
	  CrossSection_t &eeCS, TVaried_t var, int nSample, int doSave)
{
  std::vector<TH1D*> rndCSVec;
  int res=1;
  if (var!=_varRhoFile)
    res=eeCS.sampleRndVec(var,nSample,rndCSVec);
  else {
    TString inpVerTag=versionName(inpVer);
    TString loadFName=Form("dir-Rho%s/dyee_rhoRndVec_%s_%d.root",
			   inpVerTag.Data(),inpVerTag.Data(),
			   nSample);
    RndVecInfo_t info(loadFName,"h1rho_var");
    res=eeCS.sampleRndVec(var,nSample,info, rndCSVec);
  }
  if (!res) return;
  TH1D* h1Avg=NULL;
  TH2D* h2Cov=NULL;
  if (!eeCS.deriveCov(rndCSVec,&h1Avg,&h2Cov)) return;
  
  TCanvas *c= eeCS.plotCrossSection("cs");
  if (!c) return;
  h1Avg->SetLineColor(kGreen+1);
  h1Avg->SetMarkerColor(kGreen+1);
  plotHistoSame(h1Avg,"cs","LPE","rnd avg");
  printRatio(eeCS.h1PreFsrCS(), h1Avg);

  TH2D* h2Corr=NULL;
  TCanvas *cx= plotCovCorr(h2Cov,"ccov",&h2Corr);
  if (!cx) std::cout << "cx is null\n";

  TH1D *h1uncFromCov= uncFromCov(h2Cov);
  TH1D *h1relUncFromCov= uncFromCov(h2Cov,eeCS.h1PreFsrCS());

  if (!doSave) {
    std::cout << "Comparing uncertainty from covariance to the uncertainty "
	      << "in the provided cross-section\n";
    std::cout << "Expectation is that the covariance uncertainties will not "
	      << "exceed those from the cross-section\n";
    TH1D *h1CSUnc= errorAsCentral(eeCS.h1PreFsrCS(),0);
    TH1D *h1CSRelUnc= errorAsCentral(eeCS.h1PreFsrCS(),1);
    histoStyle(h1CSUnc,kBlue,24);
    histoStyle(h1CSRelUnc,kBlue,24);
    plotHisto(h1uncFromCov,"cUncCmp",1,1,"hist","uncFromCov");
    plotHistoSame(h1CSUnc,"cUncCmp","P","csUnc");

    plotHisto(h1relUncFromCov,"cRelUncCmp",1,1,"hist","rel.unc.from.cov");
    plotHistoSame(h1CSRelUnc,"cRelUncCmp","P","cs.rel.unc");
  }

  if (doSave) {
    TString fname="cov_ee_" + variedVarName(var) + Form("_%d.root",nSample);
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
    if (eeCS.h1PreFsrCS()) eeCS.h1PreFsrCS()->Write("xsec");
    if (eeCS.h1Theory()) eeCS.h1Theory()->Write("xsec_RC");
    if (h1Avg) h1Avg->Write("h1xsecAvg");
    if (h2Cov) h2Cov->Write("h2Cov");
    if (h2Corr) h2Corr->Write("h2Corr");
    if (h1uncFromCov) h1uncFromCov->Write("h1uncFromCov");
    if (h1relUncFromCov) h1relUncFromCov->Write("h1relUncFromCov");
    if (c) c->Write();
    if (cx) cx->Write();
    std::vector<TCanvas*> canvV;
    TString canvList;
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


