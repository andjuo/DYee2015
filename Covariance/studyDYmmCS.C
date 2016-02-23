#include "crossSection.h"

// --------------------------------------------------------------

void work(MuonCrossSection_t &muCS, TVaried_t var, int nSample, int doSave);

// --------------------------------------------------------------

void studyDYmmCS(TVaried_t var= _varNone, int nSample=10, int doSave=0)
{
  MuonCrossSection_t muCS("muCS","v1");
  if (!muCS.load("cs_DYmm_13TeV_v1","v1")) {
    std::cout << "loading failed\n";
    return;
  }

  //TCanvas *ca=muCS.editCSa().plotCrossSection("csa");
  //TCanvas *cb=muCS.editCSb().plotCrossSection("csb");

  if (var==_varNone) {
    TCanvas *c= muCS.plotCrossSection("cs");
    if (!c) return;

    plotHisto(cloneHisto(muCS.h1Bkg(),"h1tmp1","tmp1"),"ctmp1",1,1,"hist");

    std::vector<double> bkgW;
    bkgW.reserve(muCS.bkgWUnc().GetNoElements());
    for (int i=0; i<muCS.bkgWUnc().GetNoElements(); i++) {
      bkgW.push_back(2.);
    }
    if (!muCS.recalcBkg(bkgW)) return;
    plotHistoSame(cloneHisto(muCS.h1Bkg(),"h1tmp2","tmp2"),"ctmp1","LPE");
 
    TCanvas *c2= muCS.plotCrossSection("cs2",1);
    if (!c2) return;
  }
  else {
    std::cout << "perform study\n";
    work(muCS,var,nSample,doSave);
  }
 
}


// --------------------------------------------------------------
// --------------------------------------------------------------

void work(MuonCrossSection_t &muCS, TVaried_t var, int nSample, int doSave)
{
  std::vector<TH1D*> rndCSVec, rndCSaVec, rndCSbVec;
  int res=muCS.sampleRndVec(var,nSample,rndCSVec,&rndCSaVec,&rndCSbVec);
  if (!res) return;
  TH1D* h1Avg=NULL;
  TH2D* h2Cov=NULL;
  if (!muCS.deriveCov(rndCSVec,var,&h1Avg,&h2Cov)) return;
  
  TCanvas *c= muCS.plotCrossSection("cs");
  if (!c) return;
  h1Avg->SetLineColor(kGreen+1);
  h1Avg->SetMarkerColor(kGreen+1);
  plotHistoSame(h1Avg,"cs","LPE");

  TH2D* h2Corr=NULL;
  TCanvas *cx= plotCovCorr(h2Cov,"ccov",&h2Corr);
  if (!cx) std::cout << "cx is null\n";

  if (doSave) {
    TString fname="cov_mumu_" + variedVarName(var) + Form("_%d.root",nSample);
    TFile fout(fname,"RECREATE");
    if (!fout.IsOpen()) {
      std::cout << "file <" << fname << "> could not be created\n";
      return;
    }
    std::vector<TString> dirs;
    std::vector<std::vector<TH1D*>*> hVs;
    dirs.push_back("rnd_CSa");
    hVs.push_back(&rndCSaVec);
    dirs.push_back("rnd_CSb");
    hVs.push_back(&rndCSbVec);
    dirs.push_back("rnd_CStot");
    hVs.push_back(&rndCSVec);
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
    if (muCS.h1CS()) muCS.h1CS()->Write();
    if (muCS.h1Theory()) muCS.h1Theory()->Write();
    if (h1Avg) h1Avg->Write();
    if (h2Cov) h2Cov->Write();
    if (h2Corr) h2Corr->Write();
    if (c) c->Write();
    if (cx) cx->Write();
    std::vector<TCanvas*> canvV;
    TString canvList="muCS_a muCS_b";
    canvList.Append(" cVaried" + muCS.csA().tag());
    canvList.Append(" cVaried" + muCS.csB().tag());
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


