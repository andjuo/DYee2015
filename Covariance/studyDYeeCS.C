#include "crossSection.h"

// --------------------------------------------------------------

void work(TVersion_t inpVer,
	  CrossSection_t &muCS, TVaried_t var, int nSample, int doSave);

void createSystFile(TVersion_t inpVer,
		    CrossSection_t &elCS, TVaried_t var, int doSave);

// --------------------------------------------------------------
// save=1 save all, 2 - save les

void studyDYeeCS(TVaried_t var= _varNone, int nSample=10, int doSave=0)
{
  TVersion_t inpVer=_verEl3;
  //inpVer=_verEl3mb41;
  //inpVer=_verEl3mb42;
  inpVer=_verElMay2017;
  inpVer=_verElMay2017false;
  TString inpVerTag=versionName(inpVer);

  CrossSection_t eeCS("elCS",inpVerTag,_csPreFsrFullSp,inpVer);
  if (!eeCS.load("cs_DYee_13TeV_" + inpVerTag + ".root", inpVerTag)) {
    std::cout << "loading failed\n";
    return;
  }

  if ((inpVer==_verEl3) || (inpVer==_verElMay2017) ||
      (inpVer==_verElMay2017false)) {
    eeCS.nItersDetRes(21);
    eeCS.nItersFSR(21);
    eeCS.calcCrossSection();
  }

  if (var==_varNone) {
    TCanvas *c= eeCS.plotCrossSection("cs");
    if (!c) return;
  }
  else if (var==_varRhoSyst) {
    createSystFile(inpVer,eeCS,var,doSave);
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
  int removeNegativeSignal=1;
  if ((var!=_varRhoFile) && (var!=_varRhoSystFile) &&(var!=_varRhoSystFileSymm))
    res=eeCS.sampleRndVec(var,nSample,removeNegativeSignal,rndCSVec);
  else if (var==_varRhoFile) {
    TString inpVerTag=versionName(inpVer);
    TString loadFName=Form("dir-Rho%s/dyee_rhoRndVec_%s_%d.root",
			   inpVerTag.Data(),inpVerTag.Data(),
			   nSample);
    RndVecInfo_t info(loadFName,"h1rho_var");
    res=eeCS.sampleRndVec(var,nSample,info,removeNegativeSignal, rndCSVec);
  }
  else if (var==_varRhoSystFile) {
    TString inpVerTag=versionName(inpVer);
    TString loadFName=Form("dir-RhoSyst%s/dyee_rhoRndSystVec_%s_%d.root",
			   inpVerTag.Data(),inpVerTag.Data(),
			   nSample);
    RndVecInfo_t info(loadFName,"rhoRndSyst/h1rho_rnd");
    res=eeCS.sampleRndVec(var,nSample,info,removeNegativeSignal, rndCSVec);
  }
  else if (var==_varRhoSystFileSymm) {
    TString inpVerTag=versionName(inpVer);
    TString loadFName=Form("dir-RhoSystSymm-%s/dyee_rhoRndSystVec_%s_%d.root",
			   inpVerTag.Data(),inpVerTag.Data(),
			   nSample);
    RndVecInfo_t info(loadFName,"rhoRndSyst/h1rho_rnd");
    res=eeCS.sampleRndVec(var,nSample,info,removeNegativeSignal, rndCSVec);
  }
  else {
    std::cout << "code error: this branch should not be reached\n";
    res=0;
  }
  if (!res) return;
  TH1D* h1Avg=NULL;
  TH2D* h2Cov=NULL;
  if (!eeCS.deriveCov(rndCSVec,&h1Avg,&h2Cov)) return;

  TString massLabel= niceMassAxisLabel(leptonIdx(inpVer),"",0);
  h2Cov->GetXaxis()->SetTitle(massLabel);
  h2Cov->GetYaxis()->SetTitle(massLabel);
  h2Cov->SetTitle("h2cov " + variedVarName(var));

  PlotCovCorrOpt_t ccOpt;
  ccOpt.yTitleOffset=1.8;

  
  TCanvas *c= eeCS.plotCrossSection("cs");
  if (!c) return;
  h1Avg->SetLineColor(kGreen+1);
  h1Avg->SetMarkerColor(kGreen+1);
  plotHistoSame(h1Avg,"cs","LPE","rnd avg");
  printRatio(eeCS.h1PreFsrCS(), h1Avg);

  TH2D* h2Corr=NULL;
  TCanvas *cx= plotCovCorr(h2Cov,"ccov",ccOpt,&h2Corr);
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
    TString fname="cov_ee_" + versionName(inpVer) + "_" +
      variedVarName(var) + Form("_%d.root",nSample);
    if ((inpVer==_verEl3mb41) || (inpVer==_verEl3mb42)) {
      fname.ReplaceAll("ee_",Form("ee%d_",DYtools::nMassBins));
    }
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
    canvList+=" cVaried_" + variedVarName(var) + "_" + versionName(inpVer);
    std::cout << "canvList=" << canvList << "\n";
    if (canvList.Length() && findCanvases(canvList,canvV)) {
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

void createSystFile(TVersion_t inpVer,
		    CrossSection_t &eeCS, TVaried_t var, int doSave)
{
  std::cout << "createSystFile for inpVer=" << versionName(inpVer) << "\n";

  if (var!=_varRhoSyst) {
    std::cout << "createSystFile is not ready for var="
	      << variedVarName(var) << "\n";
    return;
  }

  TCanvas *c= eeCS.plotCrossSection("cs");
  c->Update();
  TH1D *h1cs_file= cloneHisto(eeCS.h1PreFsrCS(),"h1cs_file","h1cs_file");

  std::vector<TH1D*> h1RhoV;
  TString fname="dyee-rho-syst2.root";
  TFile fin(fname);
  if (!fin.IsOpen()) {
    std::cout << "failed to open the file <" << fin.GetName() << ">\n";
    return;
  }
  for (TTnPSystType_t syst= _tnpSyst_none; syst!=_tnpSyst_last; next(syst)) {
    TString hname="h1rho_";
    if (syst==_tnpSyst_none) hname.Append("stat");
    else hname.Append(tnpSystName(syst,1));
    TString newHName=hname;
    if (syst==_tnpSyst_none) newHName.ReplaceAll("stat","orig");
    TH1D *h1= loadHisto(fin,hname,newHName,1,h1dummy);
    h1RhoV.push_back(h1);
  }
  fin.Close();

  TCanvas *cRho=NULL;
  if (1) {
    //plotHisto(eeCS.h1Rho(),"cRho",1,0,"LPE","cs version");
    h1RhoV[0]->GetYaxis()->SetTitleOffset(1.6);
    logAxis(h1RhoV[0],1,"M_{ee} [GeV]","","\\langle\\rho\\rangle\\text{ with alt.effs}");
    cRho=plotHisto(h1RhoV[0],"cRho",1,0,"LP","orig");
    setLeftRightMargins(cRho,0.14,0.05);
    for (unsigned int i=1; i<h1RhoV.size(); i++) {
      TString legStr=h1RhoV[i]->GetName();
      legStr.ReplaceAll("h1rho_","");
      plotHistoSame(h1RhoV[i],"cRho","LP",legStr);
    }
    moveLegend(cRho,0.45,0);
  }

  std::vector<TH1D*> h1rhoRelDiffV;
  int absValues=1;
  int relativeDiff=1;
  deriveRelSyst(NULL,h1RhoV,h1rhoRelDiffV,relativeDiff,absValues);
  if (0) {
    std::cout << " : " << h1RhoV.size() << ", " << h1rhoRelDiffV.size() << "\n";
    for (unsigned int i=0; i<h1rhoRelDiffV.size(); i++) {
      printRatio(h1RhoV[0], h1RhoV[i+1]);
      printHisto(h1rhoRelDiffV[i]);
    }
    return;
  }

  TCanvas *cRhoSyst=NULL;
  if (1) {
    TString title="uncertainty";
    if (relativeDiff) title.Prepend("relative ");
    //h1rhoRelDiffV[0]->SetTitle(title);
    h1rhoRelDiffV[0]->GetYaxis()->SetRangeUser(0,0.12);
    h1rhoRelDiffV[0]->GetYaxis()->SetTitleOffset(1.6);
    TString ytitle="\\delta\\langle\\rho\\rangle";
    if (relativeDiff) ytitle= "(" + ytitle + ")_{rel}";
    //h1rhoRelDiffV[0]->GetYaxis()->SetTitle(ytitle);
    logAxis(h1rhoRelDiffV[0],3,"M_{ee} [GeV]",ytitle,title);
    cRhoSyst=plotHisto(h1rhoRelDiffV[0],"cRhoRelDiff",1,0,"LP",
				tnpSystName(TTnPSystType_t(1),1));
    setLeftRightMargins(cRhoSyst,0.14,0.05);
    for (TTnPSystType_t syst=TTnPSystType_t(2); syst!=_tnpSyst_last; next(syst))
      {
	plotHistoSame(h1rhoRelDiffV[syst-1],"cRhoRelDiff","LP",tnpSystName(syst,1));
      }
    moveLegend(cRhoSyst,0.45,0.47);
  }

  std::vector<TH1D*> h1csV;
  eeCS.var(_varRho);
  for (unsigned int i=0; i<h1RhoV.size(); i++) {
    eeCS.h1Rho(h1RhoV[i]);
    TString hname="h1cs_" + tnpSystName(TTnPSystType_t(i));
    TH1D *h1cs= cloneHisto(eeCS.calcCrossSection(0),hname,hname);
    copyStyle(h1cs,h1RhoV[i]);
    h1csV.push_back(h1cs);
  }
  std::cout << "h1csV.size=" << h1csV.size() << "\n";

  if (1) {
    TCanvas *cRhoCS= plotHisto(h1cs_file,"cRhoCS",1,1,"LP","cs_file");
    for (unsigned int i=0; i<h1csV.size(); i++) {
      plotHistoSame(h1csV[i],"cRhoCS","LP",h1csV[i]->GetName());
    }
    cRhoCS->Update();
  }

  std::vector<TH1D*> h1csRelDiffV;
  deriveRelSyst(NULL,h1csV,h1csRelDiffV,relativeDiff,absValues);
  TCanvas *cCSRhoUnc=NULL;
  if (1) {
    logAxis(h1csRelDiffV[0],3,"M_{ee} [GeV]",
	    "(\\delta\\sigma)_{\\langle\\rho\\rangle\\,syst;rel}",
	    "relative uncertainty of the cross section");
    h1csRelDiffV[0]->GetYaxis()->SetRangeUser(0,0.15);
    h1csRelDiffV[0]->GetYaxis()->SetTitleOffset(1.6);
    cCSRhoUnc= plotHisto(h1csRelDiffV[0],"cCSRhoUnc",1,0,"LP",
			 tnpSystName(TTnPSystType_t(1),1));
    setLeftRightMargins(cCSRhoUnc,0.14,0.05);
    for (TTnPSystType_t syst=TTnPSystType_t(2); syst!=_tnpSyst_last; next(syst))
      {
	plotHistoSame(h1csRelDiffV[syst-1],"cCSRhoUnc","LP",tnpSystName(syst,1));
      }
    moveLegend(cCSRhoUnc,0.45,0.47);
  }

  if (1) {
    for (unsigned int i=0; i<h1csRelDiffV.size(); i++) {
      std::cout << "\n\ni=" << i << "\n";
      printRatio(h1csRelDiffV[i],h1rhoRelDiffV[i]);
    }
  }

  if (doSave) {
    TString fname="dyeeCS-rhoSyst2.root";
    TFile fout(fname,"RECREATE");
    if (!fout.IsOpen()) {
      std::cout << "failed to create a file <" << fout.GetName() << ">\n";
      return;
    }
    for (unsigned int ih=0; ih<h1RhoV.size(); ih++) {
      h1RhoV[ih]->Write();
    }
    for (unsigned int ih=0; ih<h1rhoRelDiffV.size(); ih++) {
      h1rhoRelDiffV[ih]->Write();
    }
    for (unsigned int ih=0; ih<h1csRelDiffV.size(); ih++) {
      h1csRelDiffV[ih]->Write();
    }
    if (cRho) cRho->Write();
    if (cRhoSyst) cRhoSyst->Write();
    if (cCSRhoUnc) cCSRhoUnc->Write();
    TObjString timeTag(DayAndTimeTag(0));
    timeTag.Write("timeTag");
    fout.Close();
    std::cout << "file <" << fout.GetName() << "> created\n";
  }


  return;
}

// --------------------------------------------------------------
