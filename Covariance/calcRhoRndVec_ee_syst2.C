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

void deriveSFUnc(const DYTnPEffColl_t &coll, const EventSpace_t &es,
		 int nToys, int testCase);

// ---------------------------------------------------------------------

void calcRhoRndVec_ee_syst2(int nToys=100, int testCase=0,
			    int limitToSrc=-1,
			    int recreateCollection=0)
{

  // prepare histo styles
  hsV.push_back(&hsMC);
  hsV.push_back(&hsDataBkgPdf);
  hsV.push_back(&hsDataSigPdf);
  hsV.push_back(&hsData);  // NLOvsLO
  hsV.push_back(&hsDataTag); // tag
  hsV.push_back(&hsDataZRange);

  TString fnameBase="dyee_effSyst_";
  int elChannelFlag=2;
  int isElectronChannel=1;
  DYTnPEffColl_t coll(elChannelFlag);

  TString collFileName= fnameBase + "Coll.root";
  TString inputEffKinds="RECO ID HLT"; // recommended to match effKinds
  TString includeEffKinds= inputEffKinds;
  TString includeEffSystSrc="bkgPdf sigPdf NLOvsLO tag";
  std::vector<TString> sources;
  addToVector(sources,includeEffSystSrc);
  if ((limitToSrc!=-1) && (limitToSrc!=1111)) {
    collFileName.ReplaceAll(".root",
			    "-" + DYtools::effSystSrc[limitToSrc] + ".root");
    includeEffSystSrc=DYtools::effSystSrc[limitToSrc];
  }

  if (recreateCollection) {
    if (!createEffCollection(fnameBase,inputEffKinds,
			     includeEffKinds,includeEffSystSrc,
			     nEta,SC_eta,nPt,Ele_pt,
			     coll,
			     collFileName,
			     1)) {
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
  coll.listNumbers(1);

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
  inpVersion= _verElMay2017false;
  TString dyeeTestFNameBase="dyee_test_dressed_";
  TString fname= dyeeTestFNameBase + versionName(inpVersion) + TString(".root");

  if (!esPostFsr.load(fname,esMainDirName)) return;
  std::cout << "loaded " << esMainDirName << " from file <" << fname << ">\n";


  // calculate the scale factors using the alternative efficiency values
  if (limitToSrc==1111) {
    std::vector<TH1D*> h1rhoSystV;
    for (int iSrc=-1; iSrc<coll.getSrcCount(); iSrc++) {
      TString tag=(iSrc==-1) ? "stat" : sources[iSrc];
      DYTnPEff_t *effSrc= (iSrc==-1) ? coll.getTnPWithStatUncPtr() :
	coll.randomizeByKind(int(DYTnPEffColl_t::_rnd_ampl),0,tag,iSrc,1,1);
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
			       hname + ";M_{ee} [GeV];\\langle\\rho\\rangle");
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
      hsV[iSrc+1]->SetStyle(h1rho_src);
      h1rhoSystV.push_back(h1rho_src);
      if (iSrc==-1) plotHisto(h1rho_src,"cRho_syst",1,0,"LP","central");
      else plotHistoSame(h1rho_src,"cRho_syst","LP",sources[iSrc]);
      //h1rho_test->SetMarkerStyle(20);
      //plotHistoSame(h1rho_test,"cRho_syst","LP","test");
    }

    if (1) {
      TFile fout("dyee-rho-syst2.root","recreate");
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
  h1rho_totUnc_incorrect->SetMarkerStyle(5);
  plotHisto(h1rho_totUnc_incorrect,"cRho_totUnc_incorrect",1,0,"LPE1","Arun eff (MC stat unc)");


  // Comparison of avg.rho for debugging purposes
  if (1) {
    std::cout << "\n\ncomparison for debug\n";
    int loadElChannelFlag=elChannelFlag;
    if (inpVersion == _verElMay2017false) loadElChannelFlag=3;
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

    TH1D *h1rhoChk= esPostFsr.calculateScaleFactor(tnpEff,0,"h1rhoChk",
						   "rho_chk;M_{ee} [GeV];#rho");
    HistoStyle_t(kRed,24,1).SetStyle(h1rhoChk);
    plotHistoSame(h1rhoChk,"cRho_totUnc_incorrect","LPE1","rho chk (eff.v3)");

    TH1D *h1rho_RC=NULL;
    TString rcLabel="";
    if (inpVersion==_verEl3) {
      h1rho_RC= loadHisto("/mnt/sdb/andriusj/v3_09092016_CovarianceMatrixInputs/ROOTFile_Input6_CrossCheck.root","h_EffSF","h1rho_RC",1,h1dummy);
      rcLabel="#rho (Ridhi v3)";
    }
    else if (inpVersion==_verElMay2017false) {
      h1rho_RC= loadHisto("/media/sf_CMSData/DY13TeV-CovInputs/v20170518_Input_Cov_ee/ROOTFile_Input6_CrossCheck.root","h_EffSF","h1rho_RC",1,h1dummy);
      rcLabel="#rho (Ridhi vMay20170518)";
    }
    if (h1rho_RC && 1) {
      HistoStyle_t(kGreen+1,7,1).SetStyle(h1rho_RC);
      plotHistoSame(h1rho_RC,"cRho_totUnc_incorrect","LPE1",rcLabel);
      //printRatio(h1rho_RC,h1rhoChk);
    }
    //return;
  }

  deriveSFUnc(coll,esPostFsr,nToys,testCase);

  return;
}

// ----------------------------------------------------------
// ----------------------------------------------------------

void deriveSFUnc(const DYTnPEffColl_t &coll, const EventSpace_t &es,
		 int nToys, int testCase)
{
  int ihChk=0;
  int iSrcChk=0;
  TH2D *h2chk=NULL;
  int doCheck=0;
  std::vector<TH1D*> h1V;
  h1V.reserve(nToys);

  for (int itoy=0; itoy<nToys; itoy++) {
    TString tag=Form("_rnd%d",itoy);

    if (doCheck) h2chk= cloneHisto(coll.getTnPWithStatUnc().h2fullList(0),
				   "h2chk"+tag,"h2chk"+tag);
    DYTnPEff_t *effRnd= coll.randomizeByKind(testCase,0,tag,-1,0,0,
				     (doCheck) ? &h2chk : NULL,ihChk,iSrcChk);

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

    if (!effRnd) {
      HERE("effRnd is null\n");
      return;
    }
    TH1D *h1rhoRnd= es.calculateScaleFactor(*effRnd,0,"h1rho"+tag,
					    "h1rho"+tag+";M_{ee} [GeV];#rho");
    if (!h1rhoRnd) {
      HERE("h1rhoRnd is null\n");
      return;
    }
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
  if (h2cov) {
    plotCovCorr(h2cov,"cCov");
  }

  TH1D* h1sf=es.calculateScaleFactor(*coll.getTnPWithTotUnc(),0,"h1sf_totUnc","h1sf_totUnc;M_{ee} [GeV];#rho");
  plotHistoAuto(h1avgSF,"cAvgSF",1,0,"LPE1",Form("from toy model %d",testCase));
  plotHistoAuto(h1sf,"cAvgSF",1,0,"LPE1","with tot unc");
}

// ---------------------------------------------------------------------
