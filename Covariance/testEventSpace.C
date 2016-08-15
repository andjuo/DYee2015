#include "DYmm13TeV_eff.h"

double *createArr(int nBins, double xmin, double xmax) {
  double *x= new double[nBins+1];
  double h=(xmax-xmin)/nBins;
  for (int i=0; i<=nBins; i++) x[i]= h*i + xmin;
  return x;
}

// -------------------------------------------------------------

void testEventSpace()
{
  TH2D *h2effBinDef= new TH2D("h2effBinDef","h2effBinDef;#eta;p_{T}",
			      10,createArr(10,-2.5,2.5),
			      10,createArr(10,0,500));
  for (int ibin=1; ibin<=h2effBinDef->GetNbinsX(); ibin++) {
    for (int jbin=1; jbin<=h2effBinDef->GetNbinsY(); jbin++) {
      h2effBinDef->SetBinContent(ibin,jbin, ibin+jbin);
    }
  }
  plotHisto(h2effBinDef,"cEffBinDef",0,0);
  EventSpace_t es("mainES",h2effBinDef);

  TLorentzVector v1,v2;
  v1.SetPtEtaPhiM(16.,-2.4, 1.2, 0.);
  v2.SetPtEtaPhiM(30.,2.4,  1.2, 0.);
  std::cout << "mass= " << (v1+v2).M() << "; w=1.5\n";
  es.fill(&v1,&v2, 1.5);
  v1.SetPtEtaPhiM(200,-0.2, 1.2, 0.);
  v2.SetPtEtaPhiM(50, 0.2, 1.2, 0.);
  std::cout << "mass= " << (v1+v2).M() << "; w=1.\n";
  es.fill(&v1,&v2, 1.);

  if (1) {
    v1.SetPtEtaPhiM(16.,-2.1, 1.2, 0.);
    v2.SetPtEtaPhiM(30.,2.1,  1.2, 0.);
    std::cout << "mass= " << (v1+v2).M() << "; w=1\n";
    es.fill(&v1,&v2, 1.0);
  }

  std::vector<TH2D*> h2ptSpace, h2etaSpace;
  std::vector<TH1D*> avgH1V=es.avgAxisValues(&h2ptSpace,&h2etaSpace);

  for (unsigned int im=0; im<h2ptSpace.size(); im++) {
    if (h2ptSpace[im]->Integral()!=0) {
      plotHisto(h2ptSpace[im],"cPtSpace"+DYtools::massStr(im),0,0);
    }
  }
  for (unsigned int im=0; im<h2etaSpace.size(); im++) {
    if (h2etaSpace[im]->Integral()!=0) {
      plotHisto(h2etaSpace[im],"cEtaSpace"+DYtools::massStr(im),0,0);
    }
  }

  histoStyle(avgH1V[0],kGreen+1,5);
  histoStyle(avgH1V[1],kGreen+1,5);
  histoStyle(avgH1V[2],kBlue,24);
  histoStyle(avgH1V[3],kBlue,24);
  plotHisto(avgH1V[1],"cAvgPt",1,0,"LPE","avgPt1");
  plotHistoSame(avgH1V[3],"cAvgPt","LPE","avgPt2");
  plotHisto(avgH1V[0],"cAvgEta",1,0,"LPE","avg|Eta1|");
  plotHistoSame(avgH1V[2],"cAvgEta","LPE","avg|Eta2|");

  for (unsigned int ih=0; ih<avgH1V.size(); ih++) {
    printHisto(avgH1V[ih]);
  }

}
