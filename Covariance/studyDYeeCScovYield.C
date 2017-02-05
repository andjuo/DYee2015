#include "crossSection.h"
#include "studyDYmmCScovYield.C"

// --------------------------------------------------------------
// --------------------------------------------------------------

void studyDYeeCScovYield(int nSample=10, int method=1, int doSave=0,
			 int plotMassCounts=0, TString extraFileTag="")
{
  TVaried_t var=_varYield;
  if (nSample<0) {
    nSample=-nSample;
    var=_varYieldPoisson;
  }
  closeCanvases(10);

  std::cout << "doSave=" << doSave << "\n";

  TVersion_t inpVer=_verEl3;
  inpVer=_verEl3mb41;
  inpVer=_verEl3mb42;
  TString inpVerTag=versionName(inpVer);

  CrossSection_t elCS("elCS",inpVerTag,_csPreFsrFullSp,inpVer);
  if (!elCS.load("cs_DYee_13TeV_" + inpVerTag+".root", inpVerTag)) {
    std::cout << "loading failed\n";
    return;
  }

  // preserve original cross section
  TH1D *h1csOrig= cloneHisto(elCS.h1PreFsrCS(),"h1csOrig","h1csOrig");
  //elCS.h1Theory(elCS.h1CS());
  //printHisto(elCS.h1Theory()); return;

  TString covFName="cov_ee_"+variedVarName(var) + Form("_%d_slim.root",nSample);
  if (inpVer==_verEl3mb41) covFName.ReplaceAll("ee_","ee41_");
  else if (inpVer==_verEl3mb42) covFName.ReplaceAll("ee_","ee42_");
  if (nSample!=2000) covFName.ReplaceAll("_slim","");
  TString cName="cVaried_" + variedVarName(var) + "_" + versionName(inpVer);
  TFile fin(covFName);
  if (!fin.IsOpen()) return;
  TCanvas *canv=(TCanvas*)fin.Get(cName);
  fin.Close();
  canv->Draw();

  std::vector<TH1D*> h1aV_all, h1aV;
  getHistosFromCanvas(canv,&h1aV_all,NULL);
  std::cout << "loaded histograms: " << h1aV_all.size() << " from CSa\n";
  if (h1aV_all.size()==0) {
    std::cout << "no histos\n";
    return;
  }

  int keepOnlyPositiveSignal=(method==3) ? 1 : 0;
  if (method==4) keepOnlyPositiveSignal=2;
  TH1D* h1bkg= cloneHisto(elCS.h1Bkg(),"h1bkg_code","h1bkg_code");

  /*
  if (1) {
    // Correct muon bkg last bin in csA
    std::cout << "\n\tcorrecting background in the last mass bin\n";
    elCS.editCSa().editH1Bkg()->SetBinContent(43,0);
    std::cout << "\n";
  }
  std::cout << "h1bkgA last bin: " << h1bkgA->GetBinContent(h1bkgA->GetNbinsX()) << "\n";
  */

  if (keepOnlyPositiveSignal==2) h1bkg->Scale(1.3);

  h1aV.reserve(h1aV_all.size());
  for (unsigned int i=0; i<h1aV_all.size(); i++) {
    if (TString(h1aV_all[i]->GetName()).Index("h1varTmp")!=-1) {
      // debug printout
      if (0 && !allGE(h1aV_all[i],h1bkg)) {
	TH1D *h1= cloneHisto(h1aV_all[i],h1aV_all[i]->GetName() + TString("_tt"),
			     "");
	h1->Add(h1bkg,-1);
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
	  (keepOnlyPositiveSignal && allGE(h1aV_all[i],h1bkg))) {
	//std::cout << "keep " << h1aV_all[i]->GetName() << "\n";
	h1aV.push_back(h1aV_all[i]);
      }
    }
  }
  std::cout << "kept " << h1aV.size() << " histograms\n";
  if (h1aV.size()==0) {
    std::cout << "no histos\n";
    return;
  }

  TH1D* h1var=elCS.getVariedHisto(var);

  if (0) {
    plotHisto(h1var,"cVar",1,1,"LPE1","orig yield");
    plotHistoSame(h1aV[3],"cVar","3");
    return;
  }

  if (plotMassCounts) {
    int nSigma=10;
    int subtract=0;
    distributionByXBin(h1aV,h1var,nSigma,"cMassRange_",
		       h1bkg,subtract);
  }

  if (method) std::cout << "method=" << method << "\n";

  int removeNegativeSignal= (method==2) ? 1:0;

  unsigned int imax=h1aV.size();
  //imax=100;
  std::vector<TH1D*> h1csV;
  TCanvas *cCSChk=NULL;
  if (0) cCSChk= plotHisto(h1csOrig,"cCSChk",1,1,"LP","cs Orig");
  for (unsigned int i=0; i<imax; i++) {
    TH1D *h1cs_tmp= elCS.calcCrossSection(var,i,h1aV[i],removeNegativeSignal);
    if (!h1cs_tmp) {
      std::cout << "got null h1cs_tmp\n";
      return;
    }
    TH1D *h1cs=cloneHisto(h1cs_tmp,Form("h1cs_%d",i),h1cs_tmp->GetTitle());
    copyStyle(h1cs,h1aV[i]);
    if (cCSChk && (i<100)) plotHistoSame(h1cs,"cCSChk",Form("i=%d",int(i)));
    h1csV.push_back(h1cs);
  }

  TString calcTag= "_allRnd";
  if (keepOnlyPositiveSignal) calcTag="_keepYieldWithPosSignal";
  if (keepOnlyPositiveSignal==2) calcTag.Append("2");
  if (removeNegativeSignal) calcTag="_nulifyNegSig";

  if (1) {
    if (!h1csOrig) {
      std::cout << "h1csOrig is null\n";
      return;
    }
    TH1D *h1csForPlot= cloneHisto(h1csOrig,"h1csForPlot","h1csForPlot");
    for (int ibin=1; ibin<h1csForPlot->GetNbinsX(); ibin++) {
      h1csForPlot->SetBinError(ibin, 0.025*h1csForPlot->GetBinContent(ibin));
    }
    //HERE("call distributionByXBin");
    double scale=1e4;
    if (keepOnlyPositiveSignal==2) scale=1e2;
    distributionByXBin(h1csV,h1csForPlot,100,"cCS"+calcTag,NULL,0,scale);
    //HERE("done distributionByXBin");
  }

  if (1) {
    TString canvName="cs";
    TCanvas *cx=elCS.plotCrossSection(canvName,0);
    if (1)
    for (unsigned int i=0; i<h1csV.size(); i++) {
      plotHistoSame(h1csV[i],canvName,"LP","");
      if (i>100) { std::cout << "stopping at 100\n"; break; }
    }
    if(!cx) std::cout << "null canvas\n";
  }

  TH1D *h1avgCS=NULL;
  TH2D *h2cov=NULL, *h2corr=NULL;
  PlotCovCorrOpt_t ccOpt;
  deriveCovariance(h1csV,calcTag,"csCov_"+variedVarName(var),&h1avgCS,&h2cov);
  plotCovCorr(h2cov,"csCov",ccOpt,&h2corr);

  TString outFName=covFName;
  if (extraFileTag) calcTag.Append(extraFileTag);
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
    //h2cov->Write("h2cov"); // needed for combination code
    h2corr->Write();
    TString outDir="dir-dyeeCScov"+variedVarName(var)+
      Form("-%d",nSample)+calcTag;
    if (inpVer==_verEl3mb41) outDir.ReplaceAll("dyee","dyee41");
    else if (inpVer==_verEl3mb42) outDir.ReplaceAll("dyee","dyee42");
    SaveCanvases("ALL",outDir,&fout);
    writeTimeTag(&fout);
    fout.Close();
  }
  else std::cout << "not saved\n";

  std::cout << "macro ran with inpVerTag=" << inpVerTag << "\n";
}


// --------------------------------------------------------------
// --------------------------------------------------------------
// --------------------------------------------------------------
// --------------------------------------------------------------
// --------------------------------------------------------------


