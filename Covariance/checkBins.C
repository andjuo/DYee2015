#include "inputs.h"
#include "DYbinning.h"

void checkBins()
{
  TString srcPath="/media/ssd/v20160214_1st_CovarianceMatrixInputs/";
  srcPath="/mnt/sdb/andriusj/v20160214_1st_CovarianceMatrixInputs/";

  TString fname= srcPath + "Input5/ROOTFile_TagProbeEfficiency.root";
  TFile fin(fname);
  if (!fin.IsOpen()) {
    std::cout << "failed to open the file <" << fname << ">\n";
    return;
  }
  TH2D *h2= loadHisto(fin, "h_2D_Eff_Iso_MC", "h2",1,h2dummy);
  fin.Close();
  if (!h2) return;

  double v=0,dv=0;
  for (int ibin=0; ibin<=h2->GetNbinsX()+1; ibin++) {
    double eta= h2->GetXaxis()->GetBinCenter(ibin);
    for (int jbin=0; jbin<=h2->GetNbinsY()+1; jbin++) {
      double pt= h2->GetYaxis()->GetBinCenter(jbin);
      std::cout << "ibin=" << ibin << ", jbin=" << jbin << ", eta=" << eta
		<< ", pt=" << pt
		<< ", flatIdx=" << DYtools::FlatIndex(h2,eta,pt)
		<< "\n";
      if (DYtools::GetValue(h2, eta,pt, v,dv, 1.,0.)) {
	std::cout << " " << v << " +- " << dv << "\n";
      }
      else std::cout << " invalid values: " << v << " +- " << dv << "\n";
    }
  }
}
