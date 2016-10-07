#include "compareVersions.C"


void compareCorrs(int theCase=1, int noH2Ratios_user=0, TString showOnly="")
{
  closeCanvases();

  InfoBundle_t m,c1,c2;

  TString fname="dyee_test_dressed_El2skim3.root";
  TString fname_AK="/home/andriusj/DY13TeV/DYanalysis-20160817/ElectronNtupler/test/Analysis_Codes/AccEff/dyee_preFSR_forAccEff_v2.root";
  TString fname_AK_steps="/home/andriusj/DY13TeV/DYanalysis-20160817/ElectronNtupler/test/Analysis_Codes/AccEff/dyee_preFSR_forAccEff_v1steps.root";
  TString fname_RC="cs_DYee_13TeV_El3.root";
  TString fname_RC_unf="/mnt/sdc/andriusj/DY13TeV/EgammaWork-20160912/ElectronNtupler/work_76x/Unfolding/RespObj_detUnfolding.root";
  m.fileName= fname;
  c1.fileName= fname;

  if (theCase==1) {
    std::cout << "compare DYee corrections\n";
    m.fileName= fname;
    m.add(_iSignal,"h1recoSel_MWPU");
    m.add(_iUnf,"h1_postFsrInAccSel_MWPU");
    c1.fileName= fname;
    c1.add(_iSignal,"RRmeas+rooUnf_detResRespPU");
    c1.add(_iUnf,"RRtrue+rooUnf_detResRespPU");
    //c2.fileName= fname;
    c2.add(_iUnf, "h1_postFsrInAccSel_DieleOk_MWPU");
  }
  else if (theCase==10) {
    m.add(_iEffPass, "h1_postFsrInAccSel_DieleOk_MWPU" + lumiScale);
    m.add(_iEffPass, "h1_postFsrInAccSel_MWPU" + lumiScale);
    m.add(_iEffTot , "h1_postFsrInAcc_MW" + lumiScale);
    //m.add(_iAccPass, "h1_postFsrInAcc_MW" + lumiScale);
    //m.add(_iAcc, "h1Acc");
    c1.fileName= fname_AK;
    c1.add(_iEffPass,"h1_eff_sumPass");
    c1.add(_iEffTot, "h1_eff_sumTot");
    c1.add(_iAccPass,"h1_acc_sumPass");
    c1.add(_iAccTot, "h1_acc_sumTot");
    c1.add(_iEff, "h1eff");
    c1.add(_iAcc, "h1acc");
  }
  else if (theCase==11) {
    m.add(_iEff, "h1EffPU");
    m.add(_iEffAcc, "h1EffPUAcc");
    c2.fileName= fname;
    c2.add(_iEff,"h1_postFsrInAccSel_DieleOk_MWPU:ratio:h1_postFsrInAcc_MW");//RC
    c2.add(_iEff, "h1_postFsrInAccSel_MWPU:ratio:h1_postFsrInAcc_MW"); // MC
    c1.fileName= fname_RC;
    c1.add(_iEff, "h1Eff");
    c1.add(_iEffAcc,"h1EffAcc");
  }
  else if (theCase==12) {
    m.add(_iSignal,"h1recoSel_MWPU" + lumiScale);
    m.add(_iUnf,"h1_postFsrInAccSel_MWPU" + lumiScale); // no eleMatching
    c1.fileName= fname_RC;
    c1.add(_iSignal,"RRmeas+detRes"); // RC meas
    c1.add(_iUnf,"RRtrue+detRes"); // RC true
  }
  else if (theCase==121) {
    m.add(_iSignal,"h1recoSel_MWPU" + lumiScale);
    m.add(_iUnf,"h1_postFsrInAccSel_MWPU" + lumiScale); // no eleMatching
    c1.fileName= fname_RC_unf;
    c1.add(_iSignal,"RRmeas+Unfold_DetectorRes1"); // RC meas
    c1.add(_iUnf,"RRtrue+Unfold_DetectorRes1"); // RC true
  }
  else if (theCase==13) {
    m.add(_iEffPass,"h1_postFsrInAccSel_MWPU");
    c1.add(_iEffPass,"h1_postFsrInAccSel_MissAndMig_MWPU");
    c2.fileName= fname;
    c2.add(_iEffPass,"h1_postFsrInAccSel_DieleOk_MWPU");
  }
  else if (theCase==131) { // same as 13. changed order
    c1.add(_iEffPass,"h1_postFsrInAccSel_MWPU");
    m.add(_iEffPass,"h1_postFsrInAccSel_MissAndMig_MWPU");
    c2.fileName= fname;
    c2.add(_iEffPass,"h1_postFsrInAccSel_DieleOk_MWPU");
  }
  else if (theCase==132) { // similar to 13. Compare to AK_steps
    c1.add(_iEffPass,"h1_postFsrInAccSel_MWPU");
    m.add(_iEffPass,"h1_postFsrInAccSel_MissAndMig_MWPU");
    c2.fileName= fname_AK_steps;
    c2.add(_iEffPass,"h1_eff_sumPass_noTrigObjMatching" + lumiInvScale);
  }
  else if (theCase==21) {
    m.fileName= fname;
    m.add(_iUnf,"RRtrue+rooUnf_detResRespPU");
    c1.fileName= fname;
    c1.add(_iUnf,"h1_postFsrInAccSel_MissAndMig_MWPU");
  }
  else if (theCase==22) {
    m.add(_iEff,"h1EffPU");
    c1.add(_iEff,"h1_postFsrInAccSel_MissAndMig_MWPU:ratio:h1_postFsrInAcc_MW");
  }
  else {
    std::cout << "not ready for theCase=" << theCase << "\n";
  }

  drawComparison(m,c1,c2,noH2Ratios_user,showOnly);
}
