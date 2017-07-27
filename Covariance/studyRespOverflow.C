#include "inputs.h"
#include "crossSection.h"

// -----------------------------------------------------------

TH1D *unfold(RooUnfoldResponse *detRes, bool useOverflow, TString tag);
TH1D *getRatio(const TH1D *h1nom, const TH1D *h1denom);

// -----------------------------------------------------------

void studyRespOverflow()
{
  closeCanvases(2);

  RooUnfoldResponse *detRes_RC= loadRooUnfoldResponse("cs_DYee_13TeV_ElMay2017.root","detRes","detRes_RC",1);
  RooUnfoldResponse *detRes_AJ= loadRooUnfoldResponse("dyee_test_dressed_ElMay2017.root","rooUnf_detResRespPU","detRes_AJ",1);
  if (!detRes_RC || !detRes_AJ) return;


  TH1D *h1unf_RC_true= unfold(detRes_RC,true,"_RC_true");
  TH1D *h1unf_RC_false=unfold(detRes_RC,false,"_RC_false");
  TH1D *h1unf_AJ_true= unfold(detRes_AJ,true,"_AJ_true");
  TH1D *h1unf_AJ_false=unfold(detRes_AJ,false,"_AJ_false");

  TH1D *h1true_RC= cloneHisto((TH1D*)detRes_RC->Htruth(),"h1true_RC","h1true_RC");
  TH1D *h1true_AJ= cloneHisto((TH1D*)detRes_AJ->Htruth(),"h1true_AJ","h1true_AJ");

  hsBlack.SetStyle(h1true_RC);
  hsGreen.SetStyle(h1unf_RC_true);
  hsBlue.SetStyle(h1unf_RC_false);
  hsBlack.SetStyle(h1true_AJ);
  hsGreen.SetStyle(h1unf_AJ_true);
  hsBlue.SetStyle(h1unf_AJ_false);

  if (1) {
    TString cNameAJ="cUnf_AJ";
    plotHistoAuto(h1true_AJ,cNameAJ,1,1,"LPE","h1true AJ");
    plotHistoAuto(h1unf_AJ_true,cNameAJ,1,1,"LPE","h1unf true AJ");
    plotHistoAuto(h1unf_AJ_false,cNameAJ,1,1,"LPE","h1unf false AJ");
    TString cNameRC="cUnf_RC";
    plotHistoAuto(h1true_RC,cNameRC,1,1,"LPE","h1true RC");
    plotHistoAuto(h1unf_RC_true,cNameRC,1,1,"LPE","h1unf true RC");
    plotHistoAuto(h1unf_RC_false,cNameRC,1,1,"LPE","h1unf false RC");

    printRatio(h1unf_AJ_true,h1unf_AJ_false);
    printRatio(h1unf_RC_true,h1unf_RC_false);
  }

  if (1) {
    TString cNameAJ="cUnfRatio_AJ";
    TString cNameRC="cUnfRatio_RC";
    TH1D *h1unfRatio_AJ_true= getRatio(h1unf_AJ_true,h1true_AJ);
    TH1D *h1unfRatio_AJ_false=getRatio(h1unf_AJ_false,h1true_AJ);
    TH1D *h1unfRatio_RC_true= getRatio(h1unf_RC_true,h1true_RC);
    TH1D *h1unfRatio_RC_false=getRatio(h1unf_RC_false,h1true_RC);

    plotRatio(h1unfRatio_AJ_true,cNameAJ,1,0,"LPE","h1true AJ");
    plotHistoAuto(h1unfRatio_AJ_false,cNameAJ,1,0,"LPE","h1false AJ");
    plotRatio(h1unfRatio_RC_true,cNameRC,1,0,"LPE","h1true RC");
    plotHistoAuto(h1unfRatio_RC_false,cNameRC,1,0,"LPE","h1false RC");
  }

}

// -----------------------------------------------------------

TH1D *unfold(RooUnfoldResponse *detRes_inp, bool useOverflow, TString tag)
{
  RooUnfoldResponse *detRes= new RooUnfoldResponse(*detRes_inp);
  detRes->UseOverflow(useOverflow);
  TH1D *h1meas=cloneHisto((TH1D*)detRes->Hmeasured(),"h1meas"+tag,"h1meas"+tag);
  printHisto(h1meas);
  RooUnfoldBayes *bayes= new RooUnfoldBayes(detRes,h1meas,21,false);
  TH1D *h1unf= cloneHisto((TH1D*)bayes->Hreco(),"h1unf"+tag,"h1unf"+tag);
  printHisto(h1unf);
  delete bayes;
  return h1unf;
}

// -----------------------------------------------------------
// -----------------------------------------------------------

TH1D *getRatio(const TH1D *h1nom, const TH1D *h1denom)
{
  TString name= h1nom->GetName() + TString("__div__") +
    TString(h1denom->GetName());
  TH1D *h1= cloneHisto(h1nom,name,name);
  h1->Divide(h1denom);
  return h1;
}

// -----------------------------------------------------------
// -----------------------------------------------------------
