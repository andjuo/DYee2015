#include "inputs.h"
#include <TLine.h>
#include <fstream>
#include "DYmm13TeV_eff.h"
//#include "plotDYeeESF.C"

// ----------------------------------------------------------------
// ----------------------------------------------------------------

void plotDYeeESF_coll(TString theSet="ID", std::string args="")
{

  int iChkSrc=0, testCase=0, demo=1, plotIdx=0, plotSame=0, save=0;
  int sfOnly=0, plotEffsAndSFs=0;
  setValue(iChkSrc,args,"chkSource=",0);
  setValue(testCase,args,"testCase=",0);
  setValue(demo,args,"demo=",1);
  setValue(plotIdx,args,"plotIdx=",0);
  setValue(plotSame,args,"plotSame=",1);
  setValue(save,args,"save=",0);
  setValue(sfOnly,args,"sfOnly=",0);
  setValue(plotEffsAndSFs,args,"plotEffsAndSFs=",0);


  int theSetIdx=-1;
  if (theSet=="RECO") theSetIdx=0;
  else if (theSet=="ID") theSetIdx=1;
  else if (theSet=="HLT") theSetIdx=2;
  else {
    std::cout << "unknown theSet=" << theSet << "\n";
    return;
  }

  const int nSrc=4;
  const TString sources[nSrc] = { "bkgPdf", "sigPdf", "NLOvsLO", "tag" };
  TString tagList=". ";
  for (int i=0; i<nSrc; i++) {
    tagList.Append("_"+sources[i]);
    if (i+1<nSrc) tagList.Append(" ");
  }
  if (iChkSrc>=nSrc) {
    std::cout << "error: iChkSrc>=nSrc\n";
    return;
  }

  DYTnPEffColl_t coll(0);
  TString effFileName="dyee_effSyst_Coll.root";
  TString activeSource="all";
  if (demo==1) {
    effFileName="dyee_effSyst_Coll-noSCGap.root";
    if (valueEquals(testCase,"12 13 14")) {
      activeSource="_bkgPdf";
      tagList=". " + activeSource;
    }
    else {
      effFileName="data_demo_TnPEffColl.root";
      tagList="_coll _collA";
    }
  }
  else {
    if ((iChkSrc<0) || (iChkSrc>=nSrc)) {
      std::cout << "bad value of iChkSrc=" << iChkSrc << "\n";
      return;
    }
    activeSource= "_" + sources[iChkSrc];
    tagList=". " + activeSource;
  }

  if (!coll.load(effFileName,"",tagList)) return;
  coll.listNumbers(0);

  HistoStyle_t hsData(kRed,  5,1,1.0,2.); //hsData.SetStyle(h2effData);
  HistoStyle_t hsDataBkgPdf(kBlue,24,1,0.8,2.); //hsDataBkgPdf.SetStyle(h2effData_bkgPdf);
  HistoStyle_t hsDataSigPdf(kGreen+1,20,1,0.8,2.); //hsDataSigPdf.SetStyle(h2effData_sigPdf);
  HistoStyle_t hsMC(kBlack,3,1,1.0,2.); //hsMC.SetStyle(h2effMC);
  HistoStyle_t hsMCNLO(kBlue,24,1,0.8,2.); //hsMCNLO.SetStyle(h2effMC_NLOvsLO);
  HistoStyle_t hsDataTag(6 ,   26,1,0.8,2.); //hsDataTag.SetStyle(h2effData_tag);
  HistoStyle_t hsDataZRange(46,   23,1,0.8,2.); //hsDataZRange.SetStyle(h2effData_Zrange);

  int allInOne=1;
  std::vector<HistoStyle_t> hsSystSrcV;
  if (allInOne) hsSystSrcV.push_back(hsData);
  hsSystSrcV.push_back(hsDataBkgPdf);
  hsSystSrcV.push_back(hsDataSigPdf);
  //hsSystSrcV.push_back(hsMCNLO);
  hsSystSrcV.push_back(hsDataTag);

  int excludeGap=1;

  std::vector<TH2D*> h2_data_var;
  std::vector<TString> label_data_var;
  h2_data_var.push_back(coll.getTnPWithStatUnc().h2vecPtr(theSetIdx)->at(0));
  label_data_var.push_back(theSet + ": data");

  if (testCase==0) {
    h2_data_var.push_back(coll.getTnPSystUncSource(iChkSrc)->h2vecPtr(theSetIdx)->at(0));
    label_data_var.push_back(theSet + ": data " + sources[iChkSrc]);
    plotEffs(h2_data_var,hsSystSrcV,label_data_var, "Data",theSet,0,1, excludeGap,0);
  }
  else if (testCase==1) {
    h2_data_var.push_back(coll.getTnPSystUncSource(iChkSrc)->h2vecPtr(theSetIdx)->at(0));
    label_data_var.push_back(theSet + ": data " + sources[iChkSrc]);
    //plotEffs(h2_data_var,hsSystSrcV,label_data_var, "Data",theSet,0,1, excludeGap);


    DYTnPEff_t *sigmaUp= coll.getTnPShiftByUnc("_shiftUp",iChkSrc,1,1);
    DYTnPEff_t *sigmaDown= coll.getTnPShiftByUnc("_shiftDown",iChkSrc,-1,-1);
    //sigmaUp->SetStyle(hsSystSrcV
    h2_data_var.push_back(sigmaUp->h2vecPtr(theSetIdx)->at(0));
    label_data_var.push_back(theSet + " data sigmaUp");
    h2_data_var.push_back(sigmaDown->h2vecPtr(theSetIdx)->at(0));
    label_data_var.push_back(theSet + " data sigmaDown");
    plotEffs(h2_data_var,hsSystSrcV,label_data_var, "Data",theSet,1,1, excludeGap,0);
    plotEffs(h2_data_var,hsSystSrcV,label_data_var, "Data",theSet,1,1, excludeGap,1);
  }
  else if (valueEquals(testCase,"2 3 4 12 13 14")) {
    if (testCase>10) testCase-=10;
    DYTnPEff_t *centralEff= new DYTnPEff_t(coll.getTnPWithStatUnc(),"_central");
    DYTnPEff_t *effAmplUU=
      coll.randomizeByKind( DYTnPEffColl_t::_rnd_ampl+(testCase-2),
			    0, "_ampl_UU", 1,1);
    DYTnPEff_t *effAmplDD=
      coll.randomizeByKind( DYTnPEffColl_t::_rnd_ampl+(testCase-2),
			    0, "_ampl_DD", -1,-1);
    DYTnPEff_t *effAmplUD=
      coll.randomizeByKind( DYTnPEffColl_t::_rnd_ampl+(testCase-2),
			    0, "_ampl_UD", 1,-1);
    DYTnPEff_t *effAmplDU=
      coll.randomizeByKind( DYTnPEffColl_t::_rnd_ampl+(testCase-2),
			    0, "_ampl_DU", -1,1);

    TString expl= "ampl ";
    if (testCase==2) expl="scale ";
    else if (testCase==3) expl="BvsE ";
    else if (testCase==4) expl="\"pT\" ";
    else expl="??? ";
    int skipGap=1;
    int plotVsEta=0; // flag value
    int plotVsPt=1;  // flag value

    if (!sfOnly) {
      if (plotSame) {
	centralEff->plotProfiles(plotIdx, skipGap,plotVsEta,
				 effAmplUU,"central",expl+"sigma UU",
				 effAmplDD,expl+"sigma DD");
	centralEff->plotProfiles(plotIdx, skipGap,plotVsPt,
				 effAmplUU,"central",expl+"sigma UU",
				 effAmplDD,expl+"sigma DD");
      }
      else {
	centralEff->plotProfiles(plotIdx, skipGap,plotVsEta,
				 effAmplUD,"central",expl+"sigma UD",
				 effAmplDU,expl+"sigma DU");
	centralEff->plotProfiles(plotIdx, skipGap,plotVsPt,
				 effAmplUD,"central",expl+"sigma UD",
				 effAmplDU,expl+"sigma DU");
      }
    }
    if (sfOnly!=-1) {
      if (plotSame) {
	centralEff->plotSFProfiles(plotIdx, skipGap,plotVsEta,
				   effAmplUU,"central",expl+"sigma UU",
				 effAmplDD,expl+"sigma DD");
	centralEff->plotSFProfiles(plotIdx, skipGap,plotVsPt,
				   effAmplUU,"central",expl+"sigma UU",
				   effAmplDD,expl+"sigma DD");
      }
      else {
	centralEff->plotSFProfiles(plotIdx, skipGap,plotVsEta,
				   effAmplUD,"central",expl+"sigma UD",
				   effAmplDU,expl+"sigma DU");
	centralEff->plotSFProfiles(plotIdx, skipGap,plotVsPt,
				   effAmplUD,"central",expl+"sigma UD",
				   effAmplDU,expl+"sigma DU");
      }
    }
  }
  else if (valueEquals(testCase,"102 103 104")) {
    if (testCase>100) testCase-=100;
    DYTnPEff_t *centralEff= new DYTnPEff_t(coll.getTnPWithStatUnc(),"_central");
    DYTnPEff_t *effAmplUU=
      coll.randomizeByKind( DYTnPEffColl_t::_rnd_ampl+(testCase-2),
			    0, "_ampl_UU", 1,1);
    DYTnPEff_t *effAmplDD=
      coll.randomizeByKind( DYTnPEffColl_t::_rnd_ampl+(testCase-2),
			    0, "_ampl_DD", -1,-1);
    DYTnPEff_t *effAmplUD=
      coll.randomizeByKind( DYTnPEffColl_t::_rnd_ampl+(testCase-2),
			    0, "_ampl_UD", 1,-1);
    DYTnPEff_t *effAmplDU=
      coll.randomizeByKind( DYTnPEffColl_t::_rnd_ampl+(testCase-2),
			    0, "_ampl_DU", -1,1);

    effAmplUU->listNumbers();

    TString expl= "ampl ";
    if (testCase==2) expl="scale ";
    else if (testCase==3) expl="BvsE ";
    else if (testCase==4) expl="pT-dep ";
    else expl="??? ";
    int skipGap=1;
    int plotVsEta=0; // flag value
    int plotVsPt=1;  // flag value

    if (plotEffsAndSFs) {
      if (!sfOnly) {
	if (plotSame) {
	  centralEff->plotProfiles(plotIdx, skipGap,plotVsEta,
				   effAmplUU,"central",expl+"sigma UU",
				   effAmplDD,expl+"sigma DD");
	  centralEff->plotProfiles(plotIdx, skipGap,plotVsPt,
				   effAmplUU,"central",expl+"sigma UU",
				   effAmplDD,expl+"sigma DD");
	}
	else {
	  centralEff->plotProfiles(plotIdx, skipGap,plotVsEta,
				   effAmplUD,"central",expl+"sigma UD",
				   effAmplDU,expl+"sigma DU");
	  centralEff->plotProfiles(plotIdx, skipGap,plotVsPt,
				   effAmplUD,"central",expl+"sigma UD",
				   effAmplDU,expl+"sigma DU");
	}
      }
      if (sfOnly!=-1) {
	if (plotSame) {
	  centralEff->plotSFProfiles(plotIdx, skipGap,plotVsEta,
				     effAmplUU,"central",expl+"sigma UU",
				     effAmplDD,expl+"sigma DD");
	  centralEff->plotSFProfiles(plotIdx, skipGap,plotVsPt,
				     effAmplUU,"central",expl+"sigma UU",
				     effAmplDD,expl+"sigma DD");
	}
	else {
	  centralEff->plotSFProfiles(plotIdx, skipGap,plotVsEta,
				     effAmplUD,"central",expl+"sigma UD",
				     effAmplDU,expl+"sigma DU");
	  centralEff->plotSFProfiles(plotIdx, skipGap,plotVsPt,
				     effAmplUD,"central",expl+"sigma UD",
				     effAmplDU,expl+"sigma DU");
	}
      }
    }

    std::cout << "\n\tPrepare event space\n\n";
    EventSpace_t esPostFsr;
    TString esMainDirName="mainES";
    TVersion_t inpVersion= _verEl3;
    TString dyeeTestFNameBase="dyee_test_dressed_";
    TString fname= dyeeTestFNameBase + versionName(inpVersion) + TString(".root");

    if (!esPostFsr.load(fname,esMainDirName)) return;
    std::cout << "loaded " << esMainDirName << " from file <" << fname << ">\n";

    TString axisTitles= ";M_{ee} [GeV];#rho";
    TH1D *h1rho_central=
      esPostFsr.calculateScaleFactor(*centralEff,0,"h1rho_central",
				     "rho nominal" + axisTitles);
    TH1D *h1rho_UU=
      esPostFsr.calculateScaleFactor(*effAmplUU,0,"h1rho_UU",
				     "rho UU" + axisTitles);
    TH1D *h1rho_DD=
      esPostFsr.calculateScaleFactor(*effAmplDD,0,"h1rho_DD",
				     "rho DD" + axisTitles);
    TH1D *h1rho_UD=
      esPostFsr.calculateScaleFactor(*effAmplUD,0,"h1rho_UD",
				     "rho UD" + axisTitles);
    TH1D *h1rho_DU=
      esPostFsr.calculateScaleFactor(*effAmplDU,0,"h1rho_DU",
				     "rho DU" + axisTitles);

    hsData.SetStyle(h1rho_central);
    hsDataBkgPdf.SetStyle(h1rho_UU);
    hsDataSigPdf.SetStyle(h1rho_DD);
    hsMC.SetStyle(h1rho_UD);
    hsDataZRange.SetStyle(h1rho_DU);

    std::vector<TH1D*> h1setA;
    h1setA.push_back(h1rho_central);
    h1setA.push_back(h1rho_UU);
    h1setA.push_back(h1rho_DD);
    h1setA.push_back(h1rho_UD);
    h1setA.push_back(h1rho_DU);
    double range_min=0, range_max=1.5;
    double leg_dx=0., leg_dy=0.;
    //setAutoRanges(h1setA,range_min,range_max,0,1);
    if (valueEquals(iChkSrc,"0 2 3")) {
      range_min=0.9; range_max=1.;
      leg_dx=0.45; leg_dy=0;
    }
    else if (iChkSrc==1) {
      range_min=0.8; range_max=1.1;
      leg_dx=0.45; leg_dy=-0.03;
    }
    h1rho_central->GetYaxis()->SetRangeUser(range_min,range_max);
    h1rho_central->SetTitle("Source: " + activeSource + ", variation: " + expl);
    //logAxis(h1rho_central,1);
    TString cName="cRho"+activeSource + "_" + expl;
    TCanvas *c=plotHisto(h1rho_central,cName,1,0,"LPE","nominal");
    plotHistoSame(h1rho_UU,cName,"L",activeSource + " " + expl + " UU");
    plotHistoSame(h1rho_DD,cName,"L",activeSource + " " + expl + " DD");
    plotHistoSame(h1rho_UD,cName,"L",activeSource + " " + expl + " UD");
    plotHistoSame(h1rho_DU,cName,"L",activeSource + " " + expl + " DU");
    moveLegend(c,leg_dx,leg_dy);
    if (save) SaveCanvas(c,cName,"dir-rhoSyst",1);
  }
}

// --------------------------------------------------------------
// --------------------------------------------------------------
