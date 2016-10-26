#include "inputs.h"

void simpleTest(int nToys)
{

  if (0) {
    TH1D *h1= new TH1D("h1","h1",1000,-5.,5.);
    TH1D *h1b= (TH1D*)h1->Clone("h1b");

    for (int i=0; i<nToys; i++) {
      double r= gRandom->Gaus(0.5,0.5) + gRandom->Gaus(0,0.5) + gRandom->Gaus(0,0.5);
      h1->Fill(r);
      double rb= gRandom->Gaus(0.5,0.5)*gRandom->Gaus(1,1)*gRandom->Gaus(1,1);
      h1b->Fill(rb);
    }
    //histoStyle(h1b,kRed);
    h1b->SetLineColor(kRed);
    plotHisto(h1,"cx",0,1,"LPE");
    plotHistoSame(h1b,"cx","LPE");
    std::cout << "h1 : " << h1->GetMean() << " +- " << h1->GetRMS() << "\n";
    std::cout << "h1b: " << h1b->GetMean() << " +- " << h1b->GetRMS() << "\n";
  }

  if (1) {
    TH2D *h2stat= new TH2D("h2stat","h2stat",1,0,1,2,0.,2.);
    h2stat->SetBinContent(1,1, 0.5);  h2stat->SetBinError(1,1, 0.5);
    h2stat->SetBinContent(1,2, 1.25); h2stat->SetBinError(1,2, 0.75);
    TH2D *h2systSrc= new TH2D("h2systSrc","h2systSrc",1,0,1, 2,0.,2.);
    h2systSrc->SetBinContent(1,1, 1.); h2systSrc->SetBinError(1,1,10000.);
    h2systSrc->SetBinContent(1,2, 1.); h2systSrc->SetBinError(1,2,10000.);
    TH2D* h2systRel= (TH2D*)h2systSrc->Clone("h2systRel");
    assignDiffAsUnc(h2systRel,h2stat,1);
    TH2D* h2syst= (TH2D*)h2systSrc->Clone("h2syst");
    assignDiffAsUnc(h2syst,h2stat,0);
    printHisto(h2stat);
    printHisto(h2systSrc);
    printHisto(h2systRel);
    printHisto(h2syst);

    TH2D *h2systForRnd=cloneHisto(h2stat,"h2systForRnd","h2systForRnd");
    setError(h2systForRnd,h2syst);


    TH2D *h2rndA= (TH2D*)h2stat->Clone("h2rndA");
    TH2D *h2rndB= (TH2D*)h2stat->Clone("h2rndB");

    TH1D *h1A1= new TH1D("h1A1","h1A1",1000,-5.,5.);
    TH1D *h1A2= new TH1D("h1A2","h1A2",1000,-5.,5.);
    TH1D *h1B1= new TH1D("h1B1","h1B1",1000,-5.,5.);
    TH1D *h1B2= new TH1D("h1B2","h1B2",1000,-5.,5.);

    for (int itoy=0; itoy<nToys; itoy++) {
      randomizeWithinErr(h2systForRnd,h2rndA, 0,0,1);
      randomizeWithinRelErrSetCentral(h2stat,h2systRel, h2rndB,0,1);
      h1A1->Fill( h2rndA->GetBinContent(1,1) );
      h1A2->Fill( h2rndA->GetBinContent(1,2) );
      h1B1->Fill( h2rndB->GetBinContent(1,1) );
      h1B2->Fill( h2rndB->GetBinContent(1,2) );
    }

    h1B1->SetLineColor(kRed);
    h1B2->SetLineColor(kRed);
    plotHisto(h1A1,"c1",0,1,"from abs unc");
    plotHistoSame(h1B1,"c1","from rel unc");
    plotHisto(h1A2,"c2",0,1,"from abs unc");
    plotHistoSame(h1B2,"c2","from rel unc");
  }
}
