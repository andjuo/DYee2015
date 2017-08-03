#include "inputs.h"
#include "crossSection.h"
#include "DYmm13TeV_eff.h"

HistoStyle_t hsData(kRed,  5,1,1.0,2.);
HistoStyle_t hsDataBkgPdf(kBlue,24,1,0.8,2.);
HistoStyle_t hsDataSigPdf(kGreen+1,20,1,0.8,2.);
HistoStyle_t hsMC(kBlack,3,1,1.0,2.);
HistoStyle_t hsMCNLO(kBlue,24,1,0.8,2.);
HistoStyle_t hsDataTag(6 ,   26,1,0.8,2.);
HistoStyle_t hsDataZRange(46,   23,1,0.8,2.);
std::vector<HistoStyle_t*> hsV;

// ---------------------------------------------------------------------

void deriveSFUnc(const DYTnPEffColl_t &coll, const EventSpace_t &es,
		 int nToys, int rndProfile, int MConly=0, TString foutName="");

int createEffCollection_mm(TString effFileName,
			   TString effSystFileName,
			   const TString inpEffFileNameTags,
			   const TString includeEffKinds,
			   const TString includeEffSystSrc,
			   DYTnPEffColl_t &effColl, // create effColl
			   TString saveFileName="");

// ---------------------------------------------------------------------

void calcRhoRndVec_mm_syst2(int nToys=100, int rndProfileKind=0,
			    int limitToSrc=-1, int saveRndRhoVec=0,
			    int recreateCollection=0)
{
  closeCanvases(10);

  // prepare histo styles
  hsV.push_back(&hsMC);
  hsV.push_back(&hsDataBkgPdf);
  hsV.push_back(&hsDataSigPdf);
  hsV.push_back(&hsData);  // NLOvsLO
  hsV.push_back(&hsDataTag); // tag
  hsV.push_back(&hsDataZRange);
  hsV.push_back(&hsMCNLO);
  hsV.push_back(new HistoStyle_t(hsViolet));
  hsV.push_back(new HistoStyle_t(hsGreen));
  hsV.back()->markerStyle=27;

  TString fnameBase="dymm_effSyst_";
  int elChannelFlag=5;  // histoNames mm_v3
  int isElectronChannel=0;
  DYTnPEffColl_t coll(elChannelFlag);

  TString effFName="/mnt/sdb/andriusj/v20160527_1st_CovarianceMatrixInputs_76X/Input5/ROOTFile_TagProbeEfficiency_76X_v20160502.root";
  TString inpEffSystFName=fnameBase + "20170731.root";
  TString collFileName= fnameBase + "Coll.root";
  TString inputEffKinds="RecoID Iso HLTv4p2 HLTv4p3"; // recommended to match effKinds
  TString includeEffKinds= inputEffKinds;
  TString includeEffSystSrc="sgnChange bkgChange M60to130 M70to120 nBin30 nBin50 TagPt20 TagPt24"; // recommended to match effKinds
  std::vector<TString> sources;
  int MConly=0;
  addToVector(sources,includeEffSystSrc);
  if ((limitToSrc!=-1) && (limitToSrc!=1111)) {
    std::vector<TString> mm_effSystSrc;
    addToVector(mm_effSystSrc,includeEffSystSrc);
    collFileName.ReplaceAll(".root",
			    "-" + DYtools::effSystSrc_mm[limitToSrc] + ".root");
    includeEffSystSrc=DYtools::effSystSrc_mm[limitToSrc];
    if (includeEffSystSrc.Index("NLOvsLO")!=-1) MConly=1;
  }

  if (recreateCollection) {
    if (!createEffCollection_mm(effFName,inpEffSystFName,
				inputEffKinds,
				includeEffKinds,includeEffSystSrc,
				coll,
				collFileName)) {
      std::cout << "failure\n";
      return;
    }
    std::cout << "\n\tcollection recreated\n\n";
  }
  else {
    coll.setElChannel(isElectronChannel);
    TString tagList=". " + includeEffSystSrc;
    tagList.ReplaceAll(" "," _");
    if (!coll.load(collFileName,"",tagList)) {
      std::cout << "failed to load the collection\n";
      return;
    }
    std::cout << "\n\tcollection file <" << collFileName << "> loaded\n";
  }


  std::cout << "\ncollection lists numbers\n";
  coll.listNumbers();

  // test full error
  DYTnPEff_t *totUnc= coll.getTnPWithTotUnc();
  //totUnc->listNumbers();

  //
  // Prepare the event space
  std::cout << "\n  Prepare event space\n\n";
  EventSpace_t esPostFsr;
  TString esMainDirName="mainES";
  TVersion_t inpVersion= _verMuMay2017;
  TString dymmTestFNameBase="dymm_test_RECO_";
  TString fname= dymmTestFNameBase + versionName(inpVersion) + TString(".root");

  if (!esPostFsr.load(fname,esMainDirName)) return;
  std::cout << "loaded " << esMainDirName << " from file <" << fname << ">\n";


  // calculate the scale factors using the alternative efficiency values
  if (limitToSrc==1111) {
    std::vector<TH1D*> h1rhoSystV;
    for (int iSrc=-1; iSrc<coll.getSrcCount(); iSrc++) {
      TString tag=(iSrc==-1) ? "stat" : sources[iSrc];
      DYTnPEff_t *effSrc= (iSrc==-1) ? coll.getTnPWithStatUncPtr() :
	coll.randomizeByKind(int(DYTnPEffColl_t::_rnd_ampl),tag,iSrc,1,1);
      if (!effSrc) {
	std::cout << "null effSrc for iSrc=" << iSrc << "\n";
	return;
      }
      std::cout << "iSrc=" << iSrc << ", listing efficiencies\n";
      effSrc->listNumbers();
      //tnpEffReadyV[iSrc+1]->listNumbersBrief(1,5);
      std::cout << "numbers listed\n";
      TString hname="h1rho_" + tag;
      //TH1D *h1rho_test=
      //esPostFsr.calculateScaleFactor(*tnpEffReadyV[iSrc+1],0,hname+"_test",
      //			       hname);
      //h1rho_test->SetStats(0);
      TH1D *h1rho_src=
	esPostFsr.calculateScaleFactor(*effSrc,0,hname,
		       hname + ";M_{#mu#mu} [GeV];\\langle\\rho\\rangle");
      if (!h1rho_src) {
	std::cout << "null histo h1rho_src\n";
	return;
      }
      else {
	std::cout << "h1rho_src ok\n";
	printHisto(h1rho_src);
	//printRatio(h1rho_src,h1rho_test,0,0,1);
      }
      h1rho_src->GetYaxis()->SetRangeUser(0.84,0.98);
      h1rho_src->SetStats(0);
      if (iSrc+1>=int(hsV.size())) {
	std::cout << "hsV vector is too small\n";
	return;
      }
      hsV[iSrc+1]->SetStyle(h1rho_src);
      h1rhoSystV.push_back(h1rho_src);
      if (iSrc==-1) plotHisto(h1rho_src,"cRho_syst",1,0,"LP","central");
      else plotHistoSame(h1rho_src,"cRho_syst","LP",sources[iSrc]);
      //h1rho_test->SetMarkerStyle(20);
      //plotHistoSame(h1rho_test,"cRho_syst","LP","test");
    }

    if (1) {
      HERE("saving file");
      TFile fout("dymm-rho-syst2.root","recreate");
      for (unsigned int i=0; i<h1rhoSystV.size(); i++) {
	h1rhoSystV[i]->Write();
      }
      TObjString timeTag(DayAndTimeTag(0));
      timeTag.Write("timeTag");
      fout.Close();
      std::cout << "file <" << fout.GetName() << "> created\n";
    }

    HERE("stopping macro");
    return;
  }

  // calculate the rho
  int checkHLTv4p3=1;
  TString periodTag= (checkHLTv4p3) ? "_HLTv4p3" : "_HLTv4p2";
  TString checkCanvasName="cRho_totUnc_incorrect"+periodTag;
  TH1D *h1rho_totUnc_incorrect=
    esPostFsr.calculateScaleFactor(*totUnc,checkHLTv4p3,
				   "h1rho_totUnc_incorrect",
			   "rho new eff (MCstatUnc);M_{#mu#mu} [GeV];#rho");
  h1rho_totUnc_incorrect->GetYaxis()->SetRangeUser(0.9,1.0);
  h1rho_totUnc_incorrect->SetStats(0);
  h1rho_totUnc_incorrect->SetMarkerStyle(5);
  plotHisto(h1rho_totUnc_incorrect,checkCanvasName,1,0,"LPE1","test eff");


  // Comparison of avg.rho for debugging purposes
  if (1) {
    std::cout << "\n\ncomparison for debug\n";
    int loadElChannelFlag=isElectronChannel;
    //if (inpVersion == _verElMay2017false) loadElChannelFlag=3;
    DYTnPEff_t tnpEff(loadElChannelFlag);
    if (!tnpEff.load(fname)) return;
    tnpEff.listNumbers();

    if (0) {
      compareRanges(esPostFsr.h2EffBinDef(), tnpEff.h2Eff_RecoID_Data,1);
      compareRanges(esPostFsr.h2EffBinDef(), totUnc->h2Eff_RecoID_Data,1);
    }
    if (0) {
      totUnc->printEffRatios(tnpEff,0,1e-3);
      //for (unsigned int i=0; i<tnpEff.fullListSize(); i++) {
	//compareRanges(tnpEff.h2fullList(i), totUnc->h2fullList(i),1);
	//printRatio(tnpEff.h2fullList(i),totUnc->h2fullList(i),0,0,1e-3);
      //}
    }

    TH1D *h1rhoChk= esPostFsr.calculateScaleFactor(tnpEff,checkHLTv4p3,
						   "h1rhoChk",
					   "rho_chk;M_{#mu#mu} [GeV];#rho");
    HistoStyle_t(kRed,24,1).SetStyle(h1rhoChk);
    plotHistoSame(h1rhoChk,checkCanvasName,"LPE1","rho chk (eff.vMay2017)");

    TH1D *h1rho_KPL=NULL;
    TString kplLabel="";
    if (inpVersion==_verMuMay2017) {
      h1rho_KPL= loadAsymmGraphAsTH1D("/media/sf_CMSData/DY13TeV-CovInputs/v20170504_Input_Cov_mm/ROOTFile_Input6_CrossCheck.root",
				      "g_EffSF"+periodTag,
			      "h1rho_KPL"+periodTag,"h1rho_KPL"+periodTag,1);
      kplLabel="#rho (KPL" + periodTag + ")";
    }
    if (h1rho_KPL && 1) {
      HistoStyle_t(kGreen+1,7,1).SetStyle(h1rho_KPL);
      plotHistoSame(h1rho_KPL,checkCanvasName,"LPE1",kplLabel);
      //printRatio(h1rho_KPL,h1rhoChk);
    }
    //return;
  }

  TString foutName;
  if (saveRndRhoVec) {
    foutName="dir-RhoSyst" + versionName(inpVersion) + "/";
    gSystem->mkdir(foutName);
    foutName += "dymm_rhoRndSystVec_" + versionName(inpVersion);
    if ((limitToSrc!=-1) && (limitToSrc!=1111)) {
      foutName+="_" + DYtools::effSystSrc_mm[limitToSrc];
    }
    foutName += Form("_%d.root",nToys);
  }

  deriveSFUnc(coll,esPostFsr,nToys,rndProfileKind,MConly,foutName);

  return;
}

// ----------------------------------------------------------
// -------------------------------------------------------------

int createEffCollection_mm(TString effFileName,
			   TString effSystFileName,
			   const TString inpEffFileNameTags,
			   const TString includeEffKinds,
			   const TString includeEffSystSrc,
			   DYTnPEffColl_t &effColl,
			   TString saveFileName)
{
  std::cout << "entered createEffCollection_mm\n";

  int histoNamesFlag=4;
  int histoSystNamesFlag=5;

  DYTnPEff_t effStat(histoNamesFlag);
  if (!effStat.load(effFileName,"","")) {
    std::cout << "failed to load central values\n";
    return 0;
  }
  effStat.appendTag("_stat");
  effStat.listNumbers();

  std::vector<TString> kinds,sources;
  addToVector(kinds,inpEffFileNameTags);
  addToVector(sources,includeEffSystSrc);

  if (0) {
    // testing
    DYTnPEff_t effSyst(histoSystNamesFlag);
    if (!effSyst.load(effSystFileName,"sgnChange","")) {
      std::cout << "failed to load systematics\n";
      return 0;
    }
    effSyst.appendTag("_sigPdf");
    effSyst.listNumbers();
    return 0;
  } // end of testing

  std::vector<DYTnPEff_t*> tnpEffV; // contains histograms from files
  tnpEffV.reserve(sources.size());

  // load delta_eff histos
  TFile fin(effSystFileName);
  if (!fin.IsOpen()) {
    std::cout << "Failed to open the file <" << fin.GetName() << ">\n";
    return 0;
  }

  for (unsigned int iSrc=0; iSrc<sources.size(); iSrc++) {
    DYTnPEff_t *tnpEff=new DYTnPEff_t(histoSystNamesFlag);
    tnpEffV.push_back(tnpEff);

    if (!tnpEff->load(fin,sources[iSrc],"")) {
      std::cout << "failed to load systematics for source="
		<< sources[iSrc] << "\n";
      return 0;
    }
    tnpEff->appendTag("_" + sources[iSrc]+"_diff");
    tnpEff->listNumbers();
  }

  fin.Close();
  std::cout << "file <" << fin.GetName() << "> loaded\n";

  //std::cout << "calling assignStatErr\n";
  if (!tnpEffV[0]) {
    std::cout << "tnpEffV[0] is NULL\n";
    return 0;
  }
  // test numbers
  if (1) {
    std::cout << "\nList loaded numbers\n";
    for (unsigned int i=0; i<tnpEffV.size(); i++) {
      tnpEffV[i]->listNumbers();
    }
    std::cout << "\n -- listing done\n";
  }

  // collection
  if (!effColl.assignStatErr(effStat,"")) {
    std::cout << "failed to assign stat err\n";
    return 0;
  }

  std::cout << "adding systErrSources " << sources.size() << "\n";
  for (unsigned int iSrc=0; iSrc<sources.size(); iSrc++) {
    HERE("iSrc=%d",iSrc);
    DYTnPEff_t effSrc(effStat,"_");
    if (!effSrc.add(*tnpEffV[iSrc])) {
      std::cout << "cannot add syst.src " << sources[iSrc] << "\n";
      return 0;
    }
    if (1) {
      std::cout << "\n\tadded statEff\n";
      effSrc.listNumbers();
    }
    if (!effColl.addSystErrSource(effSrc,"_"+sources[iSrc])) {
      std::cout << "failed to assign syst err for iSrc=" << iSrc
		<< " (source name=" << sources[iSrc] << ")\n";
      return 0;
    }
  }

  if (saveFileName.Length()) {
    HERE("saving the collection\n");
    if (!effColl.save(saveFileName,"")) {
      std::cout << "failed to save the collection\n";
      return 0;
    }
    std::cout << "collection saved to file <" << saveFileName << ">\n";
  }

  return 1;
}


// -------------------------------------------------------------
// ----------------------------------------------------------

void deriveSFUnc(const DYTnPEffColl_t &coll, const EventSpace_t &es,
		 int nToys, int rndKind, int MConly, TString foutName)
{
  int ihChk=(MConly) ? 1 : 0; // ihChk=0 - data, 1 - MC (eg NLOvsLO)
  int iSrcChk=0;
  TH2D *h2chk=NULL;
  int doCheck=1;
  std::vector<TH1D*> h1v42V,h1v43V;
  h1v42V.reserve(nToys);
  h1v43V.reserve(nToys);

  TH1D* h1sf=es.calculateScaleFactor(*coll.getTnPWithTotUnc(),0,"h1sf_totUnc_4p2","h1sf_totUnc;M_{#mu#mu} [GeV];#rho");
  hsColor6.SetStyle(h1sf);
  plotHistoAuto(h1sf,"cRho_totUnc_incorrect_HLTv4p2",1,0,"LPE1","h1sf_totUnc");

  for (int itoy=0; itoy<nToys; itoy++) {
    TString tag=Form("_rnd%d",itoy);

    if (doCheck) h2chk= cloneHisto(coll.getTnPWithStatUnc().h2fullList(0),
				   "h2chk"+tag,"h2chk"+tag);
    DYTnPEff_t *effRnd= coll.randomizeByKind(rndKind,tag,-1,0,0,
				     (doCheck) ? &h2chk : NULL,ihChk,iSrcChk);

    if (h2chk) {
      if ((rndKind==0) || (rndKind==2)) {
	int iEta=0;
	TH1D *h1effSlice=
	  h2chk->ProjectionY(h2chk->GetName()+tag+Form("_iEta%d",iEta),
			     iEta+1,iEta+1);
	h1effSlice->SetStats(0);
	removeError(h1effSlice);
	logAxis(h1effSlice,1);
	h1effSlice->GetXaxis()->SetTitle("p_{T}");
	h1effSlice->GetXaxis()->SetRangeUser(0,200);
	h1effSlice->SetLineColor(itoy%10+1);
	//printHisto(h1effSlice);
	plotHistoAuto(h1effSlice,"cChkEta",1,0,"L");
      }
      if ((rndKind==0) || (rndKind==1)) {
	int iPt=0;
	TH1D *h1effSlice=
	  h2chk->ProjectionX(h2chk->GetName()+tag+Form("_iPt%d",iPt),
			     iPt+1,iPt+1);
	h1effSlice->SetStats(0);
	h1effSlice->GetXaxis()->SetTitle("#eta");
	removeError(h1effSlice);
	//logAxis(h1effSlice,1);
	//h1effSlice->GetXaxis()->SetRangeUser(0,200);
	h1effSlice->SetLineColor(itoy%10+1);
	plotHistoAuto(h1effSlice,"cChkPt",0,0,"L");
      }
    }

    if (!effRnd) {
      HERE("effRnd is null\n");
      return;
    }

    TH1D *h1rhoRnd42= es.calculateScaleFactor(*effRnd,0,"h1rho4p2"+tag,
				    "h1rho4p2"+tag+";M_{#mu#mu} [GeV];#rho");
    TH1D *h1rhoRnd43= es.calculateScaleFactor(*effRnd,1,"h1rho4p3"+tag,
				    "h1rho4p3"+tag+";M_{#mu#mu} [GeV];#rho");
    if (!h1rhoRnd42 || !h1rhoRnd43) {
      HERE("h1rhoRnd42 and/or h1rhoRnd43 is null\n");
      return;
    }

    removeError(h1rhoRnd42);
    h1rhoRnd42->SetStats(0);
    h1rhoRnd42->SetLineColor(itoy%10+1);
    logAxis(h1rhoRnd42,1);
    if (itoy<20) plotHistoAuto(h1rhoRnd42,"cRhoRnd42",1,0,"L");

    removeError(h1rhoRnd43);
    h1rhoRnd43->SetStats(0);
    h1rhoRnd43->SetLineColor(itoy%10+1);
    logAxis(h1rhoRnd43,1);
    if (itoy<20) plotHistoAuto(h1rhoRnd43,"cRhoRnd43",1,0,"L");

    h1v42V.push_back(h1rhoRnd42);
    h1v43V.push_back(h1rhoRnd43);
  }

  TH1D *h1avgSF42=NULL, *h1avgSF43=NULL;
  TH2D *h2cov42=NULL, *h2cov43=NULL;
  if (!deriveCovariance(h1v42V,
			"h1sf42_avg_corr",
			Form("sf HLTv4p2 from corr iKind=%d",rndKind),
			&h1avgSF42,&h2cov42) ||
      !deriveCovariance(h1v43V,
			"h1sf43_avg_corr",
			Form("sf_HLTv4p3 from corr iKind=%d",rndKind),
			&h1avgSF43,&h2cov43)) {
    std::cout << "failed to derive the covariance\n";
    return;
  }

  TCanvas *cCov42=NULL, *cCov43=NULL;
  if (h2cov42) {
    cCov42=plotCovCorr(h2cov42,"cCov_HLTv4p2");
  }
  if (h2cov43) {
    cCov43=plotCovCorr(h2cov43,"cCov_HLTv4p3");
  }

  TH1D* h1sf42_last=es.calculateScaleFactor(*coll.getTnPWithTotUnc(),0,
	    "h1sf42_totUnc_last","h1sf42_totUnc_last;M_{#mu#mu} [GeV];#rho");
  TH1D* h1sf43_last=es.calculateScaleFactor(*coll.getTnPWithTotUnc(),1,
	    "h1sf43_totUnc_last","h1sf43_totUnc_last;M_{#mu#mu} [GeV];#rho");
  hsColor46.SetStyle(h1sf42_last);
  hsBlue2.SetStyle(h1sf43_last);
  plotHistoAuto(h1avgSF42,"cAvgSF42",1,0,"LPE1",
		Form("from toy model %d",rndKind));
  plotHistoAuto(h1sf,"cAvgSF42",1,0,"LPE1","with tot unc");
  plotHistoAuto(h1sf42_last,"cAvgSF42",1,0,"LPE1","with tot unc (last)");

  plotHistoAuto(h1avgSF43,"cAvgSF43",1,0,"LPE1",
		Form("from toy model %d",rndKind));
  plotHistoAuto(h1sf43_last,"cAvgSF43",1,0,"LPE1","with tot unc");

  if (foutName.Length()) {
    TFile fout(foutName,"recreate");
    if (!fout.IsOpen()) {
      std::cout << "failed to create a file <" << fout.GetName() << ">\n";
      return;
    }
    h1sf->Write(); h1sf->SetDirectory(0);
    h1avgSF42->Write(); h1avgSF42->SetDirectory(0);
    h1avgSF43->Write(); h1avgSF43->SetDirectory(0);
    if (h2cov42) h2cov42->Write();
    if (h2cov43) h2cov43->Write();
    fout.mkdir("rhoRndSyst");
    fout.cd("rhoRndSyst");
    writeHistos(h1v42V,NULL); // do not provide file ptr
    writeHistos(h1v43V,NULL); // do not provide file ptr
    fout.cd();
    if (cCov42) cCov42->Write();
    if (cCov43) cCov43->Write();
    writeTimeTag(&fout);
    fout.Close();
    std::cout << "file <" << fout.GetName() << "> saved\n";
  }

}

// ---------------------------------------------------------------------
