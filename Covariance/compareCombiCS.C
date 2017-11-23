#include "inputs.h"

void compareCombiCS()
{
  TString fname1="dyll-combi-_corr_wLumi_plotChCov_inpYieldUnc_noRhoSystCorr.root";
  TString label1="no #rho syst.corr";
  TString fname2="dyll-combi-_corr_wLumi_plotChCov_inpYieldUnc_wRhoSystCorr.root";
  TString label2="with #rho syst.corr";
  if (0) {
    fname2="dyll-combi-_corr_wLumi_plotChCov_inpYieldUnc-eeRhoSystSymm.root";
    label2="with #rho syst.corr (ee Symm)";
  }

  TH1D *h1cs1= loadHisto(fname1,"h1Combi","h1Combi"+label1,h1dummy);
  TH1D *h1cs2= loadHisto(fname2,"h1Combi","h1Combi"+label2,h1dummy);
  if (!h1cs1 || !h1cs2) return;

  printRatio(h1cs1,h1cs2,0,0,1e-3);
}
