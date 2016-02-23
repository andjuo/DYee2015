#ifndef Inputs_H
#define Inputs_H

#include <TROOT.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <TObjString.h>
#include <TFile.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TArrayD.h>
#include <TSystem.h>
#include <TRandom3.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>

// -----------------------------------------------------------

TString DayAndTimeTag(int eliminateSigns=1);

void plotHisto(TH1D* h1, TString cName, int logX=0, int logY=0, TString drawOpt="LPE");
void plotHisto(TH2D* h2, TString cName, int logX=0, int logY=0);
void plotHistoSame(TH1D *h1, TString canvName, TString drawOpt);

void printHisto(const TH1D* h1);
void printHisto(const TH2D* h1);
void printRatio(const TH1D* h1a, const TH1D* h1b);
void printField(TString keyName);

inline void plotHisto(TH1* h1, TString cName, int logX=0, int logY=0, TString drawOpt="hist")
{  plotHisto((TH1D*)h1,cName,logX,logY,drawOpt); }
inline void plotHisto(TH2* h2, TString cName, int logX=0, int logY=0)
{  plotHisto((TH2D*)h2,cName,logX,logY); }

inline void printHisto(TH1* h1)
{  printHisto((TH1D*)h1); }
inline void printHisto(TH2* h2)
{  printHisto((TH2D*)h2); }


// -----------------------------------------------------------

TH1D* perMassBinWidth(const TH1D* h1, int prnBinW=0);
TH1D *flattenHisto(const TH2D *h2, TString setName);

// -----------------------------------------------------------

void printObjStringField(TFile &f, TString keyName);

// -----------------------------------------------------------
// -----------------------------------------------------------

template<class histo_t>
histo_t* loadHisto(TFile &fin, TString histoNameOnFile, TString histoName,
		   int absenseIsError, histo_t *dummy)
{
  if (!dummy) {
    // The type is defined
  }

  if (!fin.IsOpen()) {
    std::cout << "file is not open <" << fin.GetName() << ">\n";
    return NULL;
  }

  histo_t *h= (histo_t*) fin.Get(histoNameOnFile);
  if (!h) {
    if (absenseIsError) {
      std::cout << "failed to find histo <" << histoNameOnFile << "> in file <"
		<< fin.GetName() << ">\n";
    }
    return NULL;
  }
  h->SetName(histoName);
  h->SetDirectory(0);

  return h;
}

// -----------------------------------------------------------

template<class histo_t>
histo_t* loadHisto(TString fname, TString histoNameOnFile, TString histoName,
		   histo_t *dummy)
{
  if (!dummy) {
    // The type is defined
  }

  TFile fin(fname);
  if (!fin.IsOpen()) {
    std::cout << "failed to open the file <" << fname << ">\n";
    return NULL;
  }

  histo_t *h= (histo_t*) fin.Get(histoNameOnFile);
  if (!h) {
    std::cout << "failed to find histo <" << histoNameOnFile << "> in file <"
	      << fname << ">\n";
    fin.Close();
    return NULL;
  }
  h->SetName(histoName);
  h->SetDirectory(0);
  fin.Close();
  return h;
}

// -----------------------------------------------------------

template<class histo_t>
inline
histo_t* cloneHisto(const histo_t *hSrc, TString histoName, TString histoTitle)
{
  histo_t* h=(histo_t*) hSrc->Clone(histoName);
  h->SetDirectory(0);
  h->SetName(histoName);
  h->SetTitle(histoTitle);
  return h;
}

// -----------------------------------------------------------

inline
int copyContents(TH1D *h1Dest, const TH1D *h1Src)
{
  if (h1Dest->GetNbinsX() != h1Src->GetNbinsX()) {
    std::cout << "copyContents(TH1D): number of bins is different :"
	      << " h1Dest[" << h1Dest->GetNbinsX() << "],"
	      << "h1Src[" << h1Src->GetNbinsX() << "]\n";
    return 0;

  }
  for (int ibin=1; ibin<=h1Src->GetNbinsX(); ibin++) {
    h1Dest->SetBinContent(ibin, h1Src->GetBinContent(ibin));
    h1Dest->SetBinError  (ibin, h1Src->GetBinError  (ibin));
  }
  return 1;
}

// -----------------------------------------------------------
// -----------------------------------------------------------

TH1D* loadVectorD(TString fname, TString valueField, TString errorField,
		  TString setName, TString setTitle, const TH1D *h1protoType);

// -----------------------------------------------------------
// -----------------------------------------------------------

inline
void randomizeWithinErr(const TH1D *hSrc, TH1D *hDest, int nonNegative)
{
  if (!hSrc) { std::cout << "randomizeWithinErr: hSrc is null\n"; return; }
  if (!hDest) { std::cout << "randomizeWithinErr: hDest is null\n"; return; }
  for (int ibin=1; ibin<=hSrc->GetNbinsX(); ibin++) {
    double val=	gRandom->Gaus(hSrc->GetBinContent(ibin),
			      hSrc->GetBinError(ibin));
    if (nonNegative && (val<0)) val=0;
    hDest->SetBinContent(ibin,val);
    hDest->SetBinError  (ibin,0); //hSrc->GetBinError(ibin));
  }
}


// -----------------------------------------------------------
// -----------------------------------------------------------

int deriveCovariance(const std::vector<TH1D*> &rndCS,
		     TString histoNameTag, TString histoTitle,
		     TH1D **h1avgCS_out, TH2D **h2cov_out);

TH2D* cov2corr(const TH2D* h2cov);

TCanvas *plotCovCorr(TH2D* h2cov, TString canvName,
		     TH2D** h2corr_out=NULL,
		     int autoZRangeCorr=1);

TCanvas *findCanvas(TString canvName);
int findCanvases(TString canvNames, std::vector<TCanvas*> &cV);

void SaveCanvas(TCanvas* canv, const TString &canvName, TString destDir);

// -----------------------------------------------------------
// -----------------------------------------------------------


inline
void HERE(const char *msg)
{ std::cout << msg << std::endl; }

// -----------------------------------------------------------
// -----------------------------------------------------------


#endif
