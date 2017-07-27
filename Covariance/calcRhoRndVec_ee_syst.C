#include "inputs.h"
#include "crossSection.h"
#include "DYmm13TeV_eff.h"

const int nPt=5;
const double Ele_pt[nPt+1] = { 10.,20.,30.,40.,50.,2000. };
const int nEta=10;
const double SC_eta[nEta+1] = { -2.5,-2.0,-1.566,-1.4442,-0.8,0,0.8,1.4442,1.566,2.0,2.5 };

HistoStyle_t hsData(kRed,  5,1,1.0,2.);
HistoStyle_t hsDataBkgPdf(kBlue,24,1,0.8,2.);
HistoStyle_t hsDataSigPdf(kGreen+1,20,1,0.8,2.);
HistoStyle_t hsMC(kBlack,3,1,1.0,2.);
HistoStyle_t hsMCNLO(kBlue,24,1,0.8,2.);
HistoStyle_t hsDataTag(6 ,   26,1,0.8,2.);
HistoStyle_t hsDataZRange(46,   23,1,0.8,2.);
std::vector<HistoStyle_t*> hsV;

// ---------------------------------------------------------------------

std::vector<TH1D*>* randomizeSF(const DYTnPEff_t &tnpEff,
				const EventSpace_t &es,
				int nToys, TString nameBase);
std::vector<TH1D*>* randomizeSF_Uncorr(const DYTnPEffColl_t &eff,
				       int rndSrcIdx,
				       const EventSpace_t &es,
				       int nToys, TString nameBase);

void deriveSFUnc(const DYTnPEffColl_t &coll, const EventSpace_t &es,
		 int nToys, int testCase);
void deriveSFUnc_Uncorr(const DYTnPEffColl_t &coll, const EventSpace_t &es,
			int nToys, int testCase);

// ---------------------------------------------------------------------

void calcRhoRndVec_ee_syst(int nToys=100, int testCase=0, int uncorrFlag=0,
			   int limitToSrc=-1)
{

  TString fnameBase="dyee_effSyst_";
  const int nKinds=3;
  const TString kinds[nKinds] = { "RECO", "ID", "HLT" };
  const int nSrc=4;
  const TString sources[nSrc] = { "bkgPdf", "sigPdf", "NLOvsLO", "tag" };

  hsV.push_back(&hsMC);
  hsV.push_back(&hsDataBkgPdf);
  hsV.push_back(&hsDataSigPdf);
  hsV.push_back(&hsData);  // NLOvsLO
  hsV.push_back(&hsDataTag); // tag
  hsV.push_back(&hsDataZRange);

  TH2D *h2bin=new TH2D("h2bin","h2bin;#eta;p_{T} [GeV]",
		       nEta,SC_eta,nPt,Ele_pt);
  h2bin->SetStats(0);
  h2bin->SetDirectory(0);
  h2bin->Sumw2();
  for (int ibin=1; ibin<=h2bin->GetNbinsX(); ibin++) {
    for (int jbin=1; jbin<=h2bin->GetNbinsY(); jbin++) {
      h2bin->SetBinContent(ibin,jbin, 1.);
      h2bin->SetBinError  (ibin,jbin, 0.);
    }
  }
  TH2D *h2binZero= cloneHisto(h2bin,"h2binZero","h2binZero");
  h2binZero->Reset();

  int elChannelFlag=2;
  std::vector<DYTnPEff_t*> tnpEffV;
  tnpEffV.reserve(nSrc+1);
  for (int iSrc=-1; iSrc<nSrc; iSrc++) {
    tnpEffV.push_back(new DYTnPEff_t(elChannelFlag));
    for (int iKind=0; iKind<nKinds+1; iKind++) {
      for (int isMC=0; isMC<2; isMC++) {
	TString hName= Form("h2Eff_kind%d_isMC%d_iSrc%d",
			    iKind,isMC,iSrc);
	tnpEffV.back()->setHisto(iKind,isMC,
				 cloneHisto(h2bin,hName,hName));
      }
    }
  }

  std::vector<DYTnPEff_t*> tnpEffReadyV;
  tnpEffReadyV.reserve(nSrc+1);
  for (int iSrc=-1; iSrc<nSrc; iSrc++) {
    tnpEffReadyV.push_back(new DYTnPEff_t(elChannelFlag));
    for (int iKind=0; iKind<nKinds+1; iKind++) {
      for (int isMC=0; isMC<2; isMC++) {
	TString hName= Form("h2EffReady_kind%d_isMC%d_iSrc%d",
			    iKind,isMC,iSrc);
	tnpEffReadyV.back()->setHisto(iKind,isMC,
				      cloneHisto(h2bin,hName,hName));
      }
    }
  }

  // load efficiencies
  for (int iKind=0; iKind<nKinds; iKind++) {
    TString fname= fnameBase + kinds[iKind] + ".root";
    TFile fin(fname);
    if (!fin.IsOpen()) {
      std::cout << "failed to open the file <" << fin.GetName() << ">\n";
      return;
    }
    for (int isMC=0; isMC<2; isMC++) {
      TString histoNameBase=(isMC) ? "h2effMC" : "h2effData";
      TString histoNameMem= histoNameBase + "_" + kinds[iKind];
      TString histoNameFile= histoNameBase;
      if (kinds[iKind]=="ID") histoNameFile += "_inEta";
      TH2D *h2stat= loadHisto(fin, histoNameFile, histoNameMem + "_stat",
			      1,h2dummy);
      if (h2stat) {
	std::cout << "++ loaded histo " << h2stat->GetName() << "\n";
      }
      else return;
      tnpEffV[0]->setHisto(iKind,isMC,h2stat);
      tnpEffReadyV[0]->copyHisto(iKind,isMC,h2stat);

      for (int iSrc=0; iSrc<nSrc; iSrc++) {
	tnpEffV[iSrc+1]->copyHisto(iKind,isMC,h2stat);
	tnpEffReadyV[iSrc+1]->copyHisto(iKind,isMC,h2stat);
	// also include HLTv4p3
	if (iKind==nKinds-1) {
	  tnpEffV[iSrc+1]->copyHisto(iKind+1,isMC,h2stat);
	  tnpEffReadyV[iSrc+1]->copyHisto(iKind+1,isMC,h2stat);
	}
      }

      // load systematic errors
      for (int iSrc=0; iSrc<nSrc; iSrc++) {
	histoNameFile= histoNameBase + "_" + sources[iSrc];
	if (kinds[iKind]=="ID") histoNameFile += "_inEta";
	TString histoNameSyst= histoNameMem + "_syst_" + sources[iSrc];
	TH2D *h2syst= loadHisto(fin, histoNameFile, histoNameSyst, 0,h2dummy);
	if (h2syst) {
	  std::cout << "+ loaded " << histoNameSyst << "\n";
	  tnpEffV[iSrc+1]->setHisto(iKind,isMC,h2syst);
	  tnpEffReadyV[iSrc+1]->copyHisto(iKind,isMC,h2syst);
	}
	else {
	  std::cout << "- " << histoNameFile << " not found on file\n";
	  tnpEffV[iSrc+1]->setHisto(iKind,isMC,h2binZero);
	  tnpEffReadyV[iSrc+1]->copyHisto(iKind,isMC,
				  tnpEffReadyV[0]->getHisto(iKind,isMC));
	}
      }
    }
    fin.Close();
    std::cout << "file <" << fin.GetName() << "> loaded\n";
  }

  // test numbers
  if (0) {
    std::cout << "\nList loaded numbers\n";
    for (unsigned int i=0; i<tnpEffV.size(); i++) {
      tnpEffV[i]->listNumbers();
    }
    std::cout << "\n -- listing done\n";
  }

  // collection
  DYTnPEffColl_t coll(elChannelFlag);
  if (!coll.assignStatErr(*tnpEffV[0],"")) {
    std::cout << "failed to assign stat err\n";
    return;
  }
  for (int iSrc=0; iSrc<nSrc; iSrc++) {
    if ((limitToSrc!=-1) && (limitToSrc!=1111) && (iSrc!=limitToSrc)) {
      std::cout << " not including source " << sources[iSrc] << "\n";
      continue;
    }
    if (!coll.addSystErrSource(*tnpEffV[iSrc+1],"_"+sources[iSrc])) {
      std::cout << "failed to assign syst err for iSrc=" << iSrc << "\n";
      return;
    }
  }

  std::cout << "\ncollection lists numbers\n";
  coll.listNumbers(1);
  if (0) {
    TString foutNameColl= "dyee_effSyst_Coll.root";
    if ((limitToSrc!=-1) && (limitToSrc!=1111))
      foutNameColl.ReplaceAll(".root","-" + sources[limitToSrc] + ".root");
    if (!coll.save(foutNameColl,"")) {
      std::cout << "failed to save the collection\n";
      return;
    }
    std::cout << "collection saved to file <" << foutNameColl << ">\n";
    return;
  }


  // test full error
  DYTnPEff_t *totUnc= coll.getTnPWithTotUnc();
  //totUnc->listNumbers();
  //totUnc= coll.

  //
  // Prepare the event space
  std::cout << "\n  Prepare event space\n\n";
  EventSpace_t esPostFsr;
  TString esMainDirName="mainES";
  TVersion_t inpVersion= _verEl3;
  inpVersion=_verElMay2017false;
  TString dyeeTestFNameBase="dyee_test_dressed_";
  TString fname= dyeeTestFNameBase + versionName(inpVersion) + TString(".root");

  if (!esPostFsr.load(fname,esMainDirName)) return;
  std::cout << "loaded " << esMainDirName << " from file <" << fname << ">\n";

  // calculate the scale factors using the alternative efficiency values
  if (limitToSrc==1111) {
    std::vector<TH1D*> h1rhoSystV;
    for (int iSrc=-1; iSrc<nSrc; iSrc++) {
      TString tag=(iSrc==-1) ? "stat" : sources[iSrc];
      DYTnPEff_t *effSrc= (iSrc==-1) ?
	tnpEffV[0] : coll.randomizeByKind(int(DYTnPEffColl_t::_rnd_ampl),0,
					  tag,iSrc,1,1);
	//coll.getTnPShiftByUnc(tag,iSrc+1,1,1);
      if (!effSrc) {
	std::cout << "null effSrc for iSrc=" << iSrc << "\n";
	return;
      }
      std::cout << "iSrc=" << iSrc << ", listing efficiencies\n";
      effSrc->listNumbers();
      tnpEffReadyV[iSrc+1]->listNumbersBrief(1,5);
      std::cout << "numbers listed\n";
      TString hname="h1rho_" + tag;
      TH1D *h1rho_test=
	esPostFsr.calculateScaleFactor(*tnpEffReadyV[iSrc+1],0,hname+"_test",
				       hname);
      h1rho_test->SetStats(0);
      TH1D *h1rho_src=
	esPostFsr.calculateScaleFactor(*effSrc,0,hname,
			       hname + ";M_{ee} [GeV];\\langle\\rho\\rangle");
      if (!h1rho_src) {
	std::cout << "null histo h1rho_src\n";
	return;
      }
      else {
	std::cout << "h1rho_src ok\n";
	//printHisto(h1rho_src);
	printRatio(h1rho_src,h1rho_test,0,0,1);
      }
      h1rho_src->GetYaxis()->SetRangeUser(0.84,0.98);
      h1rho_src->SetStats(0);
      hsV[iSrc+1]->SetStyle(h1rho_src);
      h1rhoSystV.push_back(h1rho_src);
      if (iSrc==-1) plotHisto(h1rho_src,"cRho_syst",1,0,"LP","central");
      else plotHistoSame(h1rho_src,"cRho_syst","LP",sources[iSrc]);
      h1rho_test->SetMarkerStyle(20);
      //plotHistoSame(h1rho_test,"cRho_syst","LP","test");
    }

    if (1) {
      TFile fout("dyee-rho-syst.root","recreate");
      for (unsigned int i=0; i<h1rhoSystV.size(); i++) {
	h1rhoSystV[i]->Write();
      }
      TObjString timeTag(DayAndTimeTag(0));
      timeTag.Write("timeTag");
      fout.Close();
      std::cout << "file <" << fout.GetName() << "> created\n";
    }

    return;
  }

  // calculate the rho
  TH1D *h1rho_totUnc_incorrect=
    esPostFsr.calculateScaleFactor(*totUnc,0,"h1rho_totUnc_incorrect",
			   "rho new eff (MCstatUnc);M_{ee} [GeV];#rho");
  h1rho_totUnc_incorrect->GetYaxis()->SetRangeUser(0.9,1.0);
  h1rho_totUnc_incorrect->SetStats(0);
  plotHisto(h1rho_totUnc_incorrect,"cRho_totUnc_incorrect",1,0,"LPE1","Arun eff (MC stat unc)");

  if (1) {
    std::cout << "comparison for debug\n";
    DYTnPEff_t tnpEff(elChannelFlag);
    if (!tnpEff.load(fname)) return;
    //tnpEff.listNumbers();

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

    TH1D *h1rhoChk= esPostFsr.calculateScaleFactor(tnpEff,0,"h1rhoChk",
						   "rho_chk;M_{ee} [GeV];#rho");
    HistoStyle_t(kRed,20,1).SetStyle(h1rhoChk);
    plotHistoSame(h1rhoChk,"cRho_totUnc_incorrect","LPE1","rho chk (eff.v3)");

    if (1) {
      TH1D *h1rho_RC= loadHisto("/mnt/sdb/andriusj/v3_09092016_CovarianceMatrixInputs/ROOTFile_Input6_CrossCheck.root","h_EffSF","h1rho_RC",1,h1dummy);
      HistoStyle_t(kGreen+1,7,1).SetStyle(h1rho_RC);
      plotHistoSame(h1rho_RC,"cRho_totUnc_incorrect","LPE1","rho (Ridhi v3)");
      //printRatio(h1rho_RC,h1rhoChk);
    }
    return;
  }

  if (uncorrFlag) deriveSFUnc_Uncorr(coll,esPostFsr,nToys,testCase);
  else deriveSFUnc(coll,esPostFsr,nToys,testCase);

  return;
}

// ----------------------------------------------------------
// ----------------------------------------------------------

// ----------------------------------------------------------

std::vector<TH1D*>* randomizeSF_Uncorr(const DYTnPEffColl_t &coll,
				       int rndSrcIdx,
				       const EventSpace_t &es,
				       int nToys, TString nameBase)
{
  std::vector<TH1D*> *h1V= new std::vector<TH1D*>();
  h1V->reserve(nToys);
  int fibin1=1, fibin2=1;
  TH1D *h1follow= new TH1D("h1follow",Form("h1follow_%d_%d",fibin1,fibin2),
			   100,-1.,3.);
  h1follow->SetDirectory(0);
  const int debugMode=0; //_locdef_debugMode;
  for (int iToy=0; iToy<nToys; iToy++) {
    TString tag=Form("_srcIdx%d_rnd%d",rndSrcIdx,iToy);
    DYTnPEff_t *rndEff= coll.randomize(rndSrcIdx,tag);
    TH1D *h1toy= NULL;
    if (debugMode) { // debug mode
      h1toy=es.calculateScaleFactor_misc(*rndEff,
					 DYTnPEff_t::_misc_hlt4p2_leadLep,
					 "h1"+tag,"h1"+tag,
					 fibin1,fibin2,h1follow);
    }
    else { // realistic mode
      h1toy= es.calculateScaleFactor(*rndEff,0,
				     "h1"+tag,"h1"+tag);
    }
    h1V->push_back(h1toy);
  }
  if (h1follow && (h1follow->Integral()!=0))
    plotHisto(h1follow,"cFollowedBin",0,1,"LPE");
  return h1V;
}

// ----------------------------------------------------------

void deriveSFUnc_Uncorr(const DYTnPEffColl_t &coll, const EventSpace_t &es,
		 int nToys, int testCase)
{
  int rndSrcIdx=testCase;
  //if (testCase==-1111) rndSrcIdx=-1111;
  TString nameBase="_rnd";
  std::vector<TH1D*> *h1V = randomizeSF_Uncorr(coll,rndSrcIdx,es,nToys,nameBase);
  DYTnPEff_t *baseEff= coll.getTnPSource(rndSrcIdx);
  std::cout << "baseEff "; baseEff->listNumbers();
  std::cout << "calculateSF from baseEff\n";

  TH1D *h1sf= NULL;
  const int debugMode= 0; //_locdef_debugMode;
  if (debugMode) { // debug mode
    h1sf=es.calculateScaleFactor_misc(*baseEff,DYTnPEff_t::_misc_hlt4p2_leadLep,"h1sf_misc","h1sf_misc");
  }
  else {
    h1sf=es.calculateScaleFactor(*baseEff,0,"h1sf","h1sf");
  }
  TH1D *h1avgSF=NULL;
  TH2D *h2cov=NULL;
  if (!deriveCovariance(*h1V,"h1sf_avg",Form("sf from src=%d",rndSrcIdx),
			&h1avgSF,&h2cov)) {
    std::cout << "failed to derive the covariance\n";
    return;
  }
  printHisto(h1sf);
  printHisto(h1avgSF);
  plotHisto(h1avgSF,"cAvgSF_unCorr",1,0,"LPE1","avg");
  HistoStyle_t(kRed,24).SetStyle(h1sf);
  plotHistoSame(h1sf,"cAvgSF_unCorr","LPE1","sf");
  return;
  printHisto(h2cov);
  TH2D* h2corr=NULL;
  plotCovCorr(h2cov,"cSFCov",PlotCovCorrOpt_t(1,1,0),&h2corr);
  printHisto(h2corr);
}

// ----------------------------------------------------------

// ---------------------------------------------------------------------

std::vector<TH1D*>* randomizeSF(const DYTnPEff_t &tnpEff,
				const EventSpace_t &es,
				int nToys, TString nameBase)
{
  /*
  std::vector<TH1D*> *h1V= new std::vector<TH1D*>();
  h1V->reserve(nToys);
  const int debugMode=0; //_locdef_debugMode;
  for (int iToy=0; iToy<nToys; iToy++) {
    TString tag=Form("_srcIdx%d_rnd%d",rndSrcIdx,iToy);
    DYTnPEff_t *rndEff= coll.randomize(rndSrcIdx,tag);
    TH1D *h1toy= NULL;
    if (debugMode) { // debug mode
      h1toy=es.calculateScaleFactor_misc(*rndEff,
					 DYTnPEff_t::_misc_hlt4p2_leadLep,
					 "h1"+tag,"h1"+tag,
					 fibin1,fibin2,h1follow);
    }
    else { // realistic mode
      h1toy= es.calculateScaleFactor(*rndEff,0,
				     "h1"+tag,"h1"+tag);
    }
    h1V->push_back(h1toy);
  }
  if (h1follow && (h1follow->Integral()!=0))
    plotHisto(h1follow,"cFollowedBin",0,1,"LPE");
  return h1V;
  */
  return NULL;
}

// ---------------------------------------------------------------------

void deriveSFUnc(const DYTnPEffColl_t &coll, const EventSpace_t &es,
		 int nToys, int testCase)
{
  int ihChk=0;
  int iSrcChk=0;
  TH2D *h2chk=NULL;
  int doCheck=1;
  std::vector<TH1D*> h1V;
  h1V.reserve(nToys);

  for (int itoy=0; itoy<nToys; itoy++) {
    TString tag=Form("_rnd%d",itoy);

    if (doCheck) h2chk= cloneHisto(coll.getTnPWithStatUnc().h2fullList(0),
				   "h2chk"+tag,"h2chk"+tag);

    DYTnPEff_t *effRnd= coll.randomizeByKind(testCase,0,tag,-1,0,0,
					     &h2chk,ihChk,iSrcChk);

    if (h2chk) {
      if ((testCase==0) || (testCase==2)) {
	int iEta=0;
	TH1D *h1effSlice=
	  h2chk->ProjectionY(h2chk->GetName()+tag+Form("_iEta%d",iEta),
			     iEta+1,iEta+1);
	h1effSlice->SetStats(0);
	removeError(h1effSlice);
	logAxis(h1effSlice,1);
	h1effSlice->GetXaxis()->SetRangeUser(0,200);
	h1effSlice->SetLineColor(itoy%10+1);
	//printHisto(h1effSlice);
	plotHistoAuto(h1effSlice,"cChk",1,0,"L");
      }
      if ((testCase==0) || (testCase==1)) {
	int iPt=0;
	TH1D *h1effSlice=
	  h2chk->ProjectionX(h2chk->GetName()+tag+Form("_iPt%d",iPt),
			     iPt+1,iPt+1);
	h1effSlice->SetStats(0);
	removeError(h1effSlice);
	//logAxis(h1effSlice,1);
	//h1effSlice->GetXaxis()->SetRangeUser(0,200);
	h1effSlice->SetLineColor(itoy%10+1);
	plotHistoAuto(h1effSlice,"cChkPt",0,0,"L");
      }
    }

    TH1D *h1rhoRnd= es.calculateScaleFactor(*effRnd,0,"h1rho"+tag,
					    "h1rho"+tag+";M_{ee} [GeV];#rho");
    removeError(h1rhoRnd);
    h1rhoRnd->SetStats(0);
    h1rhoRnd->SetLineColor(itoy%10+1);
    logAxis(h1rhoRnd,1);
    if (itoy<20) plotHistoAuto(h1rhoRnd,"cRhoRnd",1,0,"L");

    h1V.push_back(h1rhoRnd);
  }

  TH1D *h1avgSF=NULL;
  TH2D *h2cov=NULL;
  if (!deriveCovariance(h1V,
			"h1sf_avg_corr",Form("sf from corr iKind=%d",testCase),
			&h1avgSF,&h2cov)) {
    std::cout << "failed to derive the covariance\n";
    return;
  }

  TH1D* h1sf=es.calculateScaleFactor(*coll.getTnPWithTotUnc(),0,"h1sf_totUnc","h1sf_totUnc;M_{ee} [GeV];#rho");
  plotHistoAuto(h1avgSF,"cAvgSF",1,0,"LPE1",Form("from toy model %d",testCase));
  plotHistoAuto(h1sf,"cAvgSF",1,0,"LPE1","with tot unc");
}

// ---------------------------------------------------------------------
