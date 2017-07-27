#include "inputs.h"
#include "DYmm13TeV_eff.h"

void plotDYeeEff(int plotIdxOnly=-1, int skipGap=0, int plotVsPt=0)
{

  TVersion_t inpVersion=_verElMay2017;
  int hEffNamingVersion=3;
  TString fnameEff="/media/sf_CMSData/DY13TeV-CovInputs/v20170518_Input_Cov_ee/ROOTFile_Input5_TagProbeEfficiency.root";

  if (0) {
    inpVersion=_verEl3;
    hEffNamingVersion=2;
    fnameEff="/mnt/sdb/andriusj/v3_09092016_CovarianceMatrixInputs/ROOTFile_Input5_TagProbeEfficiency.root";
  }

  DYTnPEff_t eff(hEffNamingVersion);
  if (!eff.load(fnameEff,"","")) return;

  eff.plotProfiles(plotIdxOnly,skipGap,plotVsPt);

  std::cout << "\n\ninpVersionName=" << versionName(inpVersion) << "\n";

  const TH2D *h2def= eff.getHisto(0,0);
  for (int iEta=1; iEta<=h2def->GetNbinsX()+1; iEta++) {
    for (int iPt=1; iPt<=h2def->GetNbinsY()+1; iPt++) {
      int fi= DYtools::FlatIndex(h2def,iEta,iPt,1);
      double eta=h2def->GetXaxis()->GetBinCenter(iEta);
      double pt=h2def->GetYaxis()->GetBinCenter(iPt);
      int fi_chk= DYtools::FlatIndex(h2def,eta,pt,1);
      std::cout << "iEta=" << iEta << "(" << eta << "), iPt=" << iPt
		<< "(" << pt << "), fi=" << fi << " and " << fi_chk << "\n";
    }
  }
}
