#ifndef DYbinning
#define DYbinning

#include <TLorentzVector.h>
#include <TH2D.h>

namespace DYtools {

const Int_t nMassBins = 45;
const Double_t massBinEdges[nMassBins+1] =
  {15, 20, 25, 30, 35, 40, 45, 50, 55, 60,
   64, 68, 72, 76, 81, 86, 91, 96, 101, 106,
   110, 115, 120, 126, 133, 141, 150, 160, 171, 185,
   200, 220, 243, 273, 320, 380, 440, 510, 600, 700,
   830, 1000, 1200, 1500, 2000, 3000};

const Double_t minMass= massBinEdges[0];
const Double_t maxMass= massBinEdges[nMassBins];

const int EtaPtFIMax=25; // assummed maximum number of (eta,pt) bins

// -------------------------------------------------------------

inline
int InsideMassRange(Double_t m)
{ return ((m>=minMass) && (m<=maxMass)) ? 1:0; }

// -------------------------------------------------------------

inline
int InAcceptance_mm(const TLorentzVector *v1, const TLorentzVector *v2) {
  if ((fabs(v1->Eta())<2.4) && (fabs(v2->Eta())<2.4)) {
    if (((v1->Pt() > 22.) && (v2->Pt() > 10.)) ||
	((v1->Pt() > 10.) && (v2->Pt() > 22.)))
      return 1;
  }
  return 0;
}

// -------------------------------------------------------------

inline
int FlatIndex(const TH2D* h2BinDefs, double eta, double pt) {
  int ibin=h2BinDefs->GetXaxis()->FindBin(eta);
  int jbin=h2BinDefs->GetYaxis()->FindBin(pt);
  int fi= (jbin-1) + (ibin-1)*h2BinDefs->GetNbinsY();
  if ((ibin==0) || (jbin==0)) fi=-1;
  else if ((ibin==h2BinDefs->GetNbinsX()+1) ||
	   (jbin==h2BinDefs->GetNbinsY()+1)) fi=-1;
  //std::cout << "eta=" << eta << ", pt=" << pt << ", ibin=" << ibin
  //	    << ", jbin=" << jbin << ", fi=" << fi << "\n";
  return fi;
}

// -------------------------------------------------------------

inline
int FlatIndex(const TH2D* h2BinDefs, const TLorentzVector *v) {
  return FlatIndex(h2BinDefs,v->Eta(),v->Pt());
}

// -------------------------------------------------------------

inline
int GetValue(const TH2D* h2, double eta, double pt, double &v, double &err,
	     double defaultVal=0., double defaultErr=0.)
{
  int ibin=h2->GetXaxis()->FindBin(eta);
  int jbin=h2->GetYaxis()->FindBin(pt);
  int ok=1;
  if ((ibin==0) || (jbin==0)) ok=0;
  else if ((ibin==h2->GetNbinsX()+1) ||
	   (jbin==h2->GetNbinsY()+1)) ok=0;
  if (ok) {
    v=  h2->GetBinContent(ibin,jbin);
    err=h2->GetBinError  (ibin,jbin);
  }
  else {
    v=defaultVal;
    err=defaultErr;
  }
  //std::cout << "eta=" << eta << ", pt=" << pt << ", ibin=" << ibin
  //	    << ", jbin=" << jbin << ", ok=" << ok
  //	    << ": " << v << " +- " << err
  //	    << "\n";
  return ok;
}

// -------------------------------------------------------------

};


#endif
