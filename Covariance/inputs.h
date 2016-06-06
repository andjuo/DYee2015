#ifndef Inputs_H
#define Inputs_H

#include <TROOT.h>
#include <TClass.h>
#include <TCanvas.h>
#include <TLegend.h>
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
#include <TGraphAsymmErrors.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>

// -----------------------------------------------------------

typedef enum { _verUndef=0, _verMu1=100, _verMu76X=101, _verEl1=200 } TVersion_t;

TString versionName(TVersion_t);

// -----------------------------------------------------------

TString DayAndTimeTag(int eliminateSigns=1);

TCanvas* plotHisto(TH1D* h1, TString cName, int logX=0, int logY=0, TString drawOpt="LPE", TString explain="");
TCanvas* plotHisto(TH2D* h2, TString cName, int logX=0, int logY=0);
TCanvas* plotHistoSame(TH1D *h1, TString canvName, TString drawOpt, TString explain="");

void printHisto(const TH1D* h1);
void printHisto(const TH2D* h1);
void printRatio(const TH1D* h1a, const TH1D* h1b);
void printField(TString keyName);

inline void plotHisto(TH1* h1, TString cName, int logX=0, int logY=0, TString drawOpt="hist", TString explain="")
{  plotHisto((TH1D*)h1,cName,logX,logY,drawOpt,explain); }
inline void plotHisto(TH2* h2, TString cName, int logX=0, int logY=0)
{  plotHisto((TH2D*)h2,cName,logX,logY); }

inline void printHisto(TH1* h1)
{  printHisto((TH1D*)h1); }
inline void printHisto(TH2* h2)
{  printHisto((TH2D*)h2); }


template <class th1_t>
inline
void removeError(th1_t *h1)
{
  for (int ibin=1; ibin<h1->GetNbinsX(); ibin++)
    h1->SetBinError(ibin,0);
}

TH1D* errorAsCentral(const TH1D* h1, int relative=0);

// -----------------------------------------------------------

TH1D* perMassBinWidth(const TH1D* h1, int prnBinW=0);
TH1D* flattenHisto(const TH2D *h2, TString setName);

TH1D* convert(TGraphAsymmErrors *gr, TString hName, TString hTitle,
	      int plotIt=0);

// -----------------------------------------------------------

void printObjStringField(TFile &f, TString keyName);

// -----------------------------------------------------------
// -----------------------------------------------------------

const TH1D* h1dummy=NULL;
const TH2D* h2dummy=NULL;


template<class histo_t>
histo_t* loadHisto(TFile &fin, TString histoNameOnFile, TString histoName,
		   int absenseIsError, const histo_t *dummy)
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
		   const histo_t *dummy)
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

template<class histo_t, class histoTarget_t>
inline
histoTarget_t* cloneHisto(const histo_t *hSrc,
			  TString histoName, TString histoTitle,
			  const histoTarget_t *hDummy)
{
  if (hDummy) {} // do nothing. Just silence the compiler
  histoTarget_t* h=(histoTarget_t*) hSrc->Clone(histoName);
  h->SetDirectory(0);
  h->SetName(histoName);
  h->SetTitle(histoTitle);
  return h;
}

// -----------------------------------------------------------

inline
void histoStyle(TH1D* h1, int color, int markerStyle, int lineStyle=1)
{
  h1->SetLineColor(color);
  h1->SetLineStyle(lineStyle);
  h1->SetMarkerColor(color);
  h1->SetMarkerStyle(markerStyle);
}

// -----------------------------------------------------------

template<class histo1D_t>
inline
int copyContents(TH1D *h1Dest, const histo1D_t *h1Src)
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
  if (!hSrc) { std::cout << "randomizeWithinErr(1D): hSrc is null\n"; return; }
  if (!hDest) { std::cout << "randomizeWithinErr(1D): hDest is null\n"; return; }
  for (int ibin=1; ibin<=hSrc->GetNbinsX(); ibin++) {
    double val=	gRandom->Gaus(hSrc->GetBinContent(ibin),
			      hSrc->GetBinError(ibin));
    if (nonNegative && (val<0)) val=0;
    hDest->SetBinContent(ibin,val);
    hDest->SetBinError  (ibin,0); //hSrc->GetBinError(ibin));
  }
}

// -----------------------------------------------------------

inline
void randomizeWithinErr(const TH2D *hSrc, TH2D *hDest, int nonNegative,
			int poissonRnd=0)
{
  if (!hSrc) { std::cout << "randomizeWithinErr(2D): hSrc is null\n"; return; }
  if (!hDest) { std::cout << "randomizeWithinErr(2D): hDest is null\n"; return; }
  for (int ibin=1; ibin<=hSrc->GetNbinsX(); ibin++) {
    for (int jbin=1; jbin<=hSrc->GetNbinsY(); jbin++) {
      double val=  (poissonRnd) ?
	gRandom->Poisson(hSrc->GetBinContent(ibin,jbin)) :
	gRandom->Gaus(hSrc->GetBinContent(ibin,jbin),
		      hSrc->GetBinError(ibin,jbin));
      if (nonNegative && (val<0)) val=0;
      hDest->SetBinContent(ibin,jbin,val);
      hDest->SetBinError  (ibin,jbin,0); //hSrc->GetBinError(ibin,jbin));
    }
  }
}


// -----------------------------------------------------------
// -----------------------------------------------------------

TMatrixD submatrix(const TMatrixD &M, int idxMin, int idxMax);
void printRatio(TString label1, const TVectorD &v1, TString label2, const TVectorD &v2);

template<class histo1D_t>
inline
TVectorD convert2vec(const histo1D_t* h1)
{
  TVectorD v(h1->GetNbinsX());
  for (int ibin=1; ibin<=h1->GetNbinsX(); ibin++)
    v(ibin-1)= h1->GetBinContent(ibin);
  return v;
}


template<class histo2D_t>
inline
TMatrixD convert2mat(const histo2D_t* h2)
{
  TMatrixD m(h2->GetNbinsX(), h2->GetNbinsY());
  for (int ibin=1; ibin<=h2->GetNbinsX(); ibin++) {
    for (int jbin=1; jbin<=h2->GetNbinsY(); jbin++) {
      m(ibin-1,jbin-1) = h2->GetBinContent(ibin,jbin);
    }
  }
  return m;
}

// chi^2_estimate= (vec1-vec2)^T Mcov^{-1} (vec1-vec2)
double chi2estimate(const TVectorD &vec1, const TVectorD &vec2,
		    const TMatrixD &Mcov);

int deriveCovariance(const std::vector<TH1D*> &rndCS,
		     TString histoNameTag, TString histoTitle,
		     TH1D **h1avgCS_out, TH2D **h2cov_out);

//TH2D* convert(const TMatrixD &m, const
TH2D* cov2corr(const TH2D* h2cov);

// uncertainty from covariance. If h1centralVal is supplied, the returned
// uncertainty is relative
TH1D* uncFromCov(const TH2D *h2cov, const TH1D *h1centralVal=NULL,
		 int zeroCentralMeansZeroRelError=0);

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
