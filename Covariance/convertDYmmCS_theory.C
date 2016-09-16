#include <TH1D.h>
#include <TCanvas.h>
#include <TFile.h>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

const Int_t nMassBins43 = 43;
const Double_t massBinEdges43[nMassBins43+1] =
  {15, 20, 25, 30, 35, 40, 45, 50, 55, 60,
   64, 68, 72, 76, 81, 86, 91, 96, 101, 106,
   110, 115, 120, 126, 133, 141, 150, 160, 171, 185,
   200, 220, 243, 273, 320, 380, 440, 510, 600, 700,
   830, 1000, 1500, 3000};


void setXaxis(TH1D *h1)
{
  h1->GetXaxis()->SetNoExponent();
  h1->GetXaxis()->SetMoreLogLabels();
}

void convertDYmmCS_theory()
{
  TString fname="dy_mm_13tev_NNLO_WNLO_4pi_twiki.txt";
  std::ifstream fin(fname);
  if (!fin.is_open()) {
    std::cout << "failed to open the file\n";
    return;
  }

  TH1D *h1cs= new TH1D("h1cs_theory","dy_mm_13tev_NNLO_WNLO_4pi_twiki;M [GeV];#sigma_{4#pi}", nMassBins43,massBinEdges43);
  h1cs->SetDirectory(0);
  TH1D *h1csCT14= new TH1D("h1csCT14_theory","dy_mm_13tev_NNLO_WNLO_4pi_twiki.txt CT14;M [GeV];#sigma_{4#pi}", nMassBins43,massBinEdges43);
  h1csCT14->SetDirectory(0);
  TH1D *h1csMMHT= new TH1D("h1csMMHT_theory","dy_mm_13tev_NNLO_WNLO_4pi_twiki.txt MMHT2014;M [GeV];#sigma_{4#pi}", nMassBins43,massBinEdges43);
  h1csMMHT->SetDirectory(0);

  std::string s;
  int skip=6;
  int im=0;
  double keepCS[3],keepCSerr[3];
  for (int i=0; i<3; i++) keepCS[i]=0;
  for (int i=0; i<3; i++) keepCSerr[i]=0;

  while (!fin.eof() && std::getline(fin,s)) {
    if (s.size()==0) continue;
    if (skip) {
      skip--;
      std::cout << "header: " << s << "\n";
      continue;
    }

    im++;
    //std::cout << ":: im=" << im << " :: " << s << "\n";
    std::stringstream ss(s);
    double mlo,mhi, cs,cs_err,csCT14,csCT14_err,csMMHT,csMMHT_err;
    ss >> mlo >> mhi >> cs >> cs_err >> csCT14 >> csCT14_err >> csMMHT >> csMMHT_err;
    int err=0;
    if (mlo==3000) break;
    if ((mlo!=massBinEdges43[im-1]) || (mhi!=massBinEdges43[im])) err=1;
    if (((mlo==1000) && (mhi==1200)) ||
	((mlo==1500) && (mhi==2000))) {
      std::cout << "chk im=" << im << ", cs=" << cs << " +- " << cs_err << "\n";
      err=2;
      keepCS[0]=cs; keepCS[1]=csCT14; keepCS[2]=csMMHT;
      keepCSerr[0]=cs_err; keepCSerr[1]=csCT14_err; keepCSerr[2]=csMMHT_err;
    }
    else if (((mlo==1200) && (mhi==1500)) ||
	     ((mlo==2000) && (mhi==3000))) {
      std::cout << "chk im=" << im << ", cs=" << cs << " +- " << cs_err << "\n";
      err=2; im--;
      keepCS[0]+=cs; keepCS[1]+=csCT14; keepCS[2]+=csMMHT;
      keepCSerr[0]=sqrt(pow(keepCSerr[0],2) + cs_err*cs_err);
      keepCSerr[1]=sqrt(pow(keepCSerr[1],2) + csCT14_err*csCT14_err);
      keepCSerr[2]=sqrt(pow(keepCSerr[2],2) + csMMHT_err*csMMHT_err);
      cs= keepCS[0]; csCT14=keepCS[1]; csMMHT=keepCS[2];
      cs_err= keepCSerr[0]; csCT14_err=keepCSerr[1]; csMMHT_err=keepCSerr[2];
    }
    if (1 || err) {
      if (err) std::cout << "error: ";
      std::cout << "im=" << im << ", mlo=" << mlo << ", mhi=" << mhi
		<< ", " << massBinEdges43[im-1] << " -- " << massBinEdges43[im] << "\n";
      if (err==1) return;
    }
    std::cout << "im=" << im << ", cs=" << cs << " +- " << cs_err << "\n";
    h1cs->SetBinContent(im, cs);
    h1cs->SetBinError(im, cs_err);
    h1csCT14->SetBinContent(im,csCT14);
    h1csCT14->SetBinError(im,csCT14_err);
    h1csMMHT->SetBinContent(im,csMMHT);
    h1csMMHT->SetBinError(im,csMMHT_err);
  }

  fin.close();

  setXaxis(h1cs);
  setXaxis(h1csCT14);
  setXaxis(h1csMMHT);

  TCanvas *cx= new TCanvas("cx","cx",600,600);
  cx->SetLogx();
  cx->SetLogy();
  h1cs->Draw("LPE1");
  h1csCT14->SetLineColor(kRed);
  h1csCT14->SetMarkerColor(kRed);
  h1csCT14->Draw("LPE1 same");
  h1csMMHT->SetLineColor(kGreen+1);
  h1csMMHT->SetMarkerColor(kGreen+1);
  h1csMMHT->Draw("LPE1 same");
  cx->Update();

  TFile fout("theory13TeVmm.root","RECREATE");
  h1cs->Write();
  //h1csCT14->Write();
  h1csMMHT->Write();
  cx->Write();
  fout.Close();
}
