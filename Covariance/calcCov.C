#include "crossSection.h"

// ---------------------------------------------------------
// ---------------------------------------------------------
// ---------------------------------------------------------

void calcCov(int nRndSamples=100)
{
  CrossSection_t cs("xs","8TeV_1D");
  if (!cs.load("cs_DYee_8TeV.root")) return;

  cs.plotCrossSection();
  TH1D *h1cs_orig= (TH1D*)cs.h1PreFsrCS()->Clone("h1cs_orig");
  //plotHistoSame(h1cs_orig,"cs","LPE");
  plotHisto(h1cs_orig,"csx",1);

  //return ;

  TVaried_t var=_varYield;
  std::vector<TH1D*> rndCS;
  if (0) {
    if (!cs.sampleRndVec(var,nRndSamples, rndCS)) {
      std::cout << "could not calculate sampleRndVec\n";
      return;
    }
  }
  else {
    std::vector<TH1D*> h1rndV;
    h1rndV.reserve(nRndSamples);
    const TH1D *h1Src= cs.getVariedHisto(var);
    if (!h1Src) {
      std::cout << "failed to get h1Src\n";
      return;
    }
    TH1D *h1Rnd=NULL;
    for (int i=0; i<nRndSamples; i++) {
      const int nonNegative=0;
      h1Rnd=(TH1D*)h1Src->Clone(Form("h1tmp_%d",i));
      randomizeWithinErr(h1Src, h1Rnd, nonNegative);
      h1rndV.push_back(h1Rnd);
    }
    if (!cs.sampleRndVec(var,h1rndV, rndCS)) {
      std::cout << "could not calculate sampleRndVec(vec)\n";
      return;
    }
  }


  TH1D *h1avgCS=NULL;
  TH2D *h2cov=NULL;
  if (!cs.deriveCov(rndCS, &h1avgCS, &h2cov)) {
    std::cout << "could not calculate the covariance\n";
    return;
  }

  if (1) {
    int itest=1;
    double avg=0;
    double sum2=0;
    TH1D *hTest=new TH1D("hTest","hTest",200,150.,350.);
    for (unsigned int i=0; i<rndCS.size(); i++) {
      const TH1D *h1=rndCS.at(i);
      avg+= h1->GetBinContent(itest);
      sum2+= h1->GetBinContent(itest)*h1->GetBinContent(itest);
      hTest->Fill(h1->GetBinContent(itest));
      std::cout << "i=" << i << ", " << itest << "st val=" << h1->GetBinContent(itest) << "\n";
    }
    avg *= 1/double(rndCS.size());
    sum2*= 1/double(rndCS.size());
    double err2= sum2-avg*avg;
    std::cout << "avg=" << avg << ", err=" << sqrt(err2) << ", err2=" << err2 << "\n";
    std::cout << "hTest: " << hTest->GetMean() << " +- " << hTest->GetRMS() << "\n";

    double sum2test=0;
    for (unsigned int i=0; i<rndCS.size(); i++) {
      const TH1D *h1=rndCS.at(i);
      sum2test += pow(h1->GetBinContent(itest)-avg, 2);
    }
    sum2test*= 1/double(rndCS.size());
    std::cout << " sum2test=" << sum2test << ", errtest=" << sqrt(sum2test) << "\n";

    std::cout << "value from a routine avg: " << h1avgCS->GetBinContent(itest) << " +- " << h1avgCS->GetBinError(itest) << "\n";
    std::cout << " cov : " << h2cov->GetBinContent(itest,itest) << "\n";
  }

  for (unsigned int i=0; i<rndCS.size(); i++) {
    plotHistoSame(rndCS[i],"csx","hist");
  }

  h1cs_orig->SetLineColor(kBlue);
  plotHistoSame(h1cs_orig,"csx","LPE");
  h1avgCS->SetLineColor(kRed);
  plotHistoSame(h1avgCS,"csx","LPE");

  TH2D *h2corr= cov2corr(h2cov);
  TCanvas *ccov= new TCanvas("ccov","ccov",1200,600);
  ccov->Divide(2,1);
  ccov->cd(1);
  gPad->SetRightMargin(0.18);
  gPad->SetLogx(); gPad->SetLogy();
  h2cov->SetStats(0);
  h2cov->Draw("COLZ");
  ccov->cd(2);
  gPad->SetRightMargin(0.18);
  gPad->SetLogx(); gPad->SetLogy();
  h2corr->Draw("COLZ");
  ccov->Update();
}

// ---------------------------------------------------------
// ---------------------------------------------------------
