#include "inputs.h"
#include "crossSection.h"

void checkMigMatrix(int fsrRes=0)
{
  TVersion_t inpVer=_verEl3;
  //inpVer=_verElMay2017false;
  TString inpVerTag=versionName(inpVer);

  TString eeCSFName="cs_DYee_13TeV_" + inpVerTag + ".root";
  TString fieldName=(fsrRes) ? "FSRRes" : "detRes";

  RooUnfoldResponse *r= loadRooUnfoldResponse(eeCSFName,fieldName,
					      fieldName,1);
  if (!r) return;

  TH2D *h2= (TH2D*)r->Hresponse();

  printHisto(h2);

  TH1D *h1diag= diagAsTH1D(h2);
  printHisto(h1diag);

  if (1) {
    std::cout << "\n\nLarge relative error\n";
    for (int ibin=1; ibin<=h2->GetNbinsX(); ibin++) {
      for (int jbin=1; jbin<=h2->GetNbinsY(); jbin++) {
	double v= h2->GetBinContent(ibin,jbin);
	double e= h2->GetBinError(ibin,jbin);
	if (e>v) {
	  std::cout << "ibin=" << ibin << ",jbin=" << jbin
		    << "  " << h2->GetXaxis()->GetBinLowEdge(ibin)
		    << "-" << (h2->GetXaxis()->GetBinLowEdge(ibin)+
			       h2->GetXaxis()->GetBinWidth(ibin))
		    << "  " << h2->GetYaxis()->GetBinLowEdge(jbin)
		    << "-" << (h2->GetYaxis()->GetBinLowEdge(jbin)+
			   h2->GetYaxis()->GetBinWidth(jbin))
		    << "  " << h2->GetBinContent(ibin,jbin) << " +- "
		    << h2->GetBinError(ibin,jbin) << "\n";
	}
      }
    }
  }
}
