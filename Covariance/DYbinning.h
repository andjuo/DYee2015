#ifndef DYbinning
#define DYbinning

#include <TLorentzVector.h>
#include <TH2D.h>
#include <stdarg.h>
#include <TMath.h>
#include <iostream>
#include <cmath>

namespace DYtools {

const Int_t nMassBins45 = 45;
const Double_t massBinEdges45[nMassBins45+1] =
  {15, 20, 25, 30, 35, 40, 45, 50, 55, 60,
   64, 68, 72, 76, 81, 86, 91, 96, 101, 106,
   110, 115, 120, 126, 133, 141, 150, 160, 171, 185,
   200, 220, 243, 273, 320, 380, 440, 510, 600, 700,
   830, 1000, 1200, 1500, 2000, 3000};

const Int_t nMassBins43 = 43;
const Double_t massBinEdges43[nMassBins43+1] =
  {15, 20, 25, 30, 35, 40, 45, 50, 55, 60,
   64, 68, 72, 76, 81, 86, 91, 96, 101, 106,
   110, 115, 120, 126, 133, 141, 150, 160, 171, 185,
   200, 220, 243, 273, 320, 380, 440, 510, 600, 700,
   830, 1000, 1500, 3000};


#ifdef __CXX__
 extern const Int_t nMassBins;
 extern const Double_t *massBinEdges;
#else
 const Int_t nMassBins= nMassBins43;
 const Double_t *massBinEdges= massBinEdges43;
#endif

const Double_t minMass= massBinEdges[0];
const Double_t maxMass= massBinEdges[nMassBins];


// electron data samples
const int nEEDatasetMassBins=11;
const int eeDatasetMass[nEEDatasetMassBins+1] = {10,50,100,200,400,500,700,800,
						 1000,1500,2000,3000};

// -------------------------------------------------------------

inline
  TString massStr(int iMass, int changeToUnderscore=0)
{
  TString s;
  if ((iMass>=0) && (iMass<nMassBins)) {
    char c=(changeToUnderscore) ? '_' : '-';
    s=Form("%1.0lf%c%1.0lf",massBinEdges[iMass],c,massBinEdges[iMass+1]);
  }
  else if (iMass<0) s=Form("lt%1.0lf",massBinEdges[0]);
  else  s=Form("gt%1.0lf",massBinEdges[nMassBins]);
  return s;
}

// -------------------------------------------------------------

inline
int massIdx(Double_t mass) {
  int idx=-1;
  if ((mass>= minMass) && (mass<maxMass)) {
    for (int i=0; i<nMassBins; i++) {
      if ((mass>=massBinEdges[i]) && (mass<massBinEdges[i+1])) {
	idx=i;
	break;
      }
    }
  }
  return idx;
}

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
int InECALSCGap(double eta)
{ return ((fabs(eta)>1.4442) && (fabs(eta)<1.566)) ? 1:0; }

// -------------------------------------------------------------

inline
int InAcceptance_ee(const TLorentzVector *v1, const TLorentzVector *v2,
		    int testSCGap=0) {
  if ((fabs(v1->Eta())<2.5) && (fabs(v2->Eta())<2.5) &&
      (!testSCGap || (!InECALSCGap(v1->Eta()) && !InECALSCGap(v2->Eta())))) {
    if (((v1->Pt() > 30.) && (v2->Pt() > 10.)) ||
	((v1->Pt() > 10.) && (v2->Pt() > 30.)))
      return 1;
  }
  return 0;
}

// -------------------------------------------------------------

inline
int compareBinning(int count, const TH2D* h2a, ...) {
  int ok=1;
  va_list vl;
  va_start(vl,h2a);
  std::cout << "compareBinning of :\n";
  for (int i=0; (i<count-1) && ok; ++i) {
    typedef const TH2D* constTH2D_ptr;
    const TH2D* h2b = (const TH2D*)(va_arg(vl,const TH2D*));
    std::cout << " " << h2a->GetName() << " to " << h2b->GetName() << "\n";
    if (h2a->GetNbinsX() != h2b->GetNbinsX()) {
      std::cout << " x binning is different\n";
      ok=0;
    }
    if (h2a->GetNbinsY() != h2b->GetNbinsY()) {
      std::cout << " y binning is different\n";
      ok=0;
    }
    if (ok) {
      for (int ibin=1; ibin<=h2a->GetNbinsX()+1; ibin++) {
	const double xa=h2a->GetXaxis()->GetBinLowEdge(ibin);
	const double xb=h2b->GetXaxis()->GetBinLowEdge(ibin);
	if (xa!=xb) {
	  std::cout << " x bin at " << ibin << " is different "
		    << xa << " vs " << xb << "\n";
	  ok=0;
	}
      }

      for (int jbin=1; jbin<=h2a->GetNbinsY()+1; jbin++) {
	const double ya=h2a->GetYaxis()->GetBinLowEdge(jbin);
	const double yb=h2b->GetYaxis()->GetBinLowEdge(jbin);
	if (ya!=yb) {
	  std::cout << " y bin at " << jbin << " is different "
		    << ya << " vs " << yb << "\n";
	}
      }
    }
    if (ok) std::cout << " -- check ok\n";
  }
  va_end(vl);
  return ok;
}

// -------------------------------------------------------------

inline
double deltaPhi(double phi1, double phi2) {
  double dphi = fabs(phi1 - phi2);
  if(dphi >= TMath::Pi()) dphi = TMath::TwoPi() - dphi;
  return dphi;
}

// -------------------------------------------------------------

inline
double deltaR(double eta1, double phi1, double eta2, double phi2)
{
  double dEta = fabs(eta1-eta2);
  double dPhi = deltaPhi(phi1, phi2);
  double dr = sqrt(dEta*dEta + dPhi*dPhi);
  return dr;
}

// -------------------------------------------------------------

};


#endif
