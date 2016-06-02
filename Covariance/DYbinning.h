#ifndef DYbinning
#define DYbinning

#include <TLorentzVector.h>
#include <TH2D.h>

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


const Int_t nMassBins= nMassBins43;
const Double_t *massBinEdges= massBinEdges43;

const Double_t minMass= massBinEdges[0];
const Double_t maxMass= massBinEdges[nMassBins];



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

};


#endif
