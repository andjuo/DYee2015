#include "inputs.h"

void plotEff()
{
  TString fname="/mnt/sdb/andriusj/v3_09092016_CovarianceMatrixInputs/ROOTFile_Input5_TagProbeEfficiency.root";
  if (1) {
    TH2D *h2reco=loadHisto(fname,"h_2D_Eff_RECO_Data","h2reco",1,h2dummy);
    plotHisto(h2reco,"cReco",0,1,"COLZ");
    //printHisto(h2reco);
    printHistoRange(h2reco);
  }
  if (1) {
    TH2D *h2id=loadHisto(fname,"h_2D_Eff_ID_Data","h2id",1,h2dummy);
    plotHisto(h2id,"cID",0,1,"COLZ");
    //printHisto(h2id);
    printHistoRange(h2id);
  }
  if (0) {
    TH2D *h2id=loadHisto(fname,"h_2D_Eff_ID_MC","h2idmc",1,h2dummy);
    plotHisto(h2id,"cID_MC",0,1,"COLZ");
    //printHisto(h2id);
    printHistoRange(h2id);
  }
  if (1) {
    TH2D *h2trg=loadHisto(fname,"h_2D_Eff_Trig_Data","h2trig",1,h2dummy);
    plotHisto(h2trg,"cHLT",0,1,"COLZ");
    //printHisto(h2trg);
    printHistoRange(h2trg);
  }
}
