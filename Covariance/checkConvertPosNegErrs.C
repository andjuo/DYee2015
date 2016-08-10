#include "DYmm13TeV.h"
#include "DYmm13TeV_eff.h"
#include "DYbinning.h"
#include "inputs.h"

// --------------------------------------------------------------
// --------------------------------------------------------------

void checkConvertPosNegErrs()
{
  TString srcPath="/mnt/sdb/andriusj/v20160527_1st_CovarianceMatrixInputs_76X/";
  TString fnameEff= srcPath + "Input5/ROOTFile_TagProbeEfficiency_76X_v20160502.root";

  DYTnPEff_t tnpEff, tnpEff_avg, tnpEff_pos, tnpEff_neg;

  TFile finEff(fnameEff);
  if (!finEff.IsOpen()) {
    std::cout << "failed to open the file <" << fnameEff << ">\n";
    return;
  }
  if (!tnpEff.load_DYMM_TnP(finEff,1,"_KL")) { HERE("failed to load"); return; }
  if (!tnpEff_avg.load_DYMM_TnP_asymmEff(finEff,1,0,"_avg")) { HERE("failed to load avg"); return; }
  if (!tnpEff_pos.load_DYMM_TnP_asymmEff(finEff,1,1,"_pos")) { HERE("failed to load pos"); return; }
  if (!tnpEff_neg.load_DYMM_TnP_asymmEff(finEff,1,-1,"_neg")) { HERE("failed to load neg"); return; }
  finEff.Close();

  if (1) {
    HERE("\ncomparing ratios");

    tnpEff.printEffRatios(tnpEff_avg,1);
    HERE("\n\n");
    tnpEff.printEffRatios(tnpEff_pos,1);
    HERE("\n\n");
    tnpEff.printEffRatios(tnpEff_neg,1);
  }

  DYTnPEff_t tnpEff_asym;
  if (!tnpEff_asym.randomize(tnpEff_pos,tnpEff_neg,"_rndAsym1")) {
    std::cout << "randomization failed\n";
    return;
  }
  tnpEff.printEffRatios(tnpEff_asym,0);
}
