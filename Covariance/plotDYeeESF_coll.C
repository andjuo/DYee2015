#include "inputs.h"
#include <TLine.h>
#include <fstream>
#include "DYmm13TeV_eff.h"
#include "plotDYeeESF.C"

// ----------------------------------------------------------------
// ----------------------------------------------------------------

void plotDYeeESF_coll(TString theSet="ID", int iChkSrc=0,
		      int testCase=0, int save=0
		      )
{

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
  if (!coll.load("dyee_effSyst_Coll.root","",tagList)) return;
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
    plotEffs(h2_data_var,hsSystSrcV,label_data_var, "Data",theSet,0, excludeGap);
  }
  else if (testCase==1) {
    h2_data_var.push_back(coll.getTnPSystUncSource(iChkSrc)->h2vecPtr(theSetIdx)->at(0));
    label_data_var.push_back(theSet + ": data " + sources[iChkSrc]);
    //plotEffs(h2_data_var,hsSystSrcV,label_data_var, "Data",theSet,0, excludeGap);


    DYTnPEff_t *sigmaUp= coll.getTnPShiftByUnc("_shiftUp",iChkSrc,1,1);
    DYTnPEff_t *sigmaDown= coll.getTnPShiftByUnc("_shiftDown",iChkSrc,-1,-1);
    //sigmaUp->SetStyle(hsSystSrcV
    h2_data_var.push_back(sigmaUp->h2vecPtr(theSetIdx)->at(0));
    label_data_var.push_back(theSet + " data sigmaUp");
    h2_data_var.push_back(sigmaDown->h2vecPtr(theSetIdx)->at(0));
    label_data_var.push_back(theSet + " data sigmaDown");
    plotEffs(h2_data_var,hsSystSrcV,label_data_var, "Data",theSet,1, excludeGap);
    plotEffs(h2_data_var,hsSystSrcV,label_data_var, "Data",theSet,1, excludeGap, 1);
  }
}

// --------------------------------------------------------------
// --------------------------------------------------------------
