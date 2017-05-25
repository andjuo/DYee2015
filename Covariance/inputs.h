#ifndef Inputs_H
#define Inputs_H

#include <TROOT.h>
#include <TClass.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <TObjString.h>
#include <TFile.h>
#include <TTree.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TArrayD.h>
#include <TSystem.h>
#include <TRandom3.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TLorentzVector.h>
#include <stdarg.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>

// -----------------------------------------------------------

typedef enum { _verUndef=0, _verMu1=100, _verMu76X=101, _verMuApproved=102,
	       _verMuMay2017=103,
	       _verEl1=200,
	       _verEl2=201, _verEl2skim=202, _verEl2skim2=203,
	       _verEl2skim3=204,
	       _verEl3=301, _verEl3mb41=302, _verEl3mb42=303,
	       _verElMay2017=304
} TVersion_t;

TString versionName(TVersion_t);
int leptonIdx(TVersion_t); // 0 - electron, 1 - muon, 2 - lepton

// -----------------------------------------------------------

struct HistoStyle_t
{
  int color, markerStyle,lineStyle;
  double markerSize, lineWidth;
public:
  HistoStyle_t(int set_color, int set_markerStyle, int set_lineStyle=1,
	       double set_markerSize=1., double set_lineWidth=1.);
  HistoStyle_t(const HistoStyle_t &hs);

  template<class graph_t>
  void SetStyle(graph_t *gr) const {
    gr->SetLineColor(color);
    gr->SetLineStyle(lineStyle);
    gr->SetLineWidth(lineWidth);
    gr->SetMarkerColor(color);
    gr->SetMarkerStyle(markerStyle);
    gr->SetMarkerSize(markerSize);
  }

  template<class graph_t>
  void operator()(graph_t *gr) const { this->SetStyle(gr); }

  void Print() const
  { std::cout << "hs: color=" << color << ", marker=" << markerStyle
	      << " (sz=" << markerSize << "), lineStyle=" << lineStyle
	      << " (width=" << lineWidth << ")\n"; }
};

// -----------------------------------------------------------

struct PlotCovCorrOpt_t
{
  int autoZRangeCorr, gridLines, logScaleX, logScaleY;
  int setTickX, setTickY;
  double yTitleOffset, leftMargin, rightMargin;
public:
  PlotCovCorrOpt_t(int set_autoZRangeCorr=1, int set_gridLines=1,
		   int set_logScale=1, double set_yTitleOffset=1.5,
		   double set_leftMargin=0.15, double set_rightMargin=0.15);
  PlotCovCorrOpt_t(const TString key, int value);
  PlotCovCorrOpt_t(const PlotCovCorrOpt_t &o);

  void setLogScale(int logX=1, int logY=1)
  { logScaleX=logX; logScaleY=logY; }
  void noLogScale() { logScaleX=0; logScaleY=0; }

  void setTicksX() { setTickX=1; }
  void setTicksY() { setTickY=1; }
  void setTicks(int setx, int sety) { setTickX=setx; setTickY=sety; }

  void assign(const PlotCovCorrOpt_t &opt);
  int setValue(TString key, int value, int verbose=1);
  void adjustTicks(TCanvas *c) const;
};

// -----------------------------------------------------------

TString DayAndTimeTag(int eliminateSigns=1);


template<class ptr_t>
inline
void clearVec(std::vector<ptr_t*> &vec)
{
  for (unsigned int i=0; i<vec.size(); i++)
    if (vec[i]) delete vec[i];
  vec.clear();
}

// assign values
void addToVector(std::vector<TString> &vec, TString strings);
void addToVector(std::vector<int> &vec, int count, ...);
void addToVector(std::vector<double> &vec, int count, ...);

TCanvas* plotHisto(TH1D* h1, TString cName, int logX=0, int logY=0, TString drawOpt="LPE", TString explain="", int gridLines=3);
TCanvas* plotHisto(TH2D* h2, TString cName, int logX=0, int logY=0,
		   double ytitleOffset=1.5, int gridLines=1,
		   double centralZRange=0.);
TCanvas* plotHisto(TH2D* h2, TString cName, const PlotCovCorrOpt_t &o);
TCanvas* plotHisto(const TH2D* h2, TString cName, int logX=0, int logY=0,
		   double ytitleOffset=1.5, int gridLines=1,
		   double centralZRange=0.);
TCanvas* plotHistoSame(TH1D *h1, TString canvName, TString drawOpt, TString explain="");
// combination of plotHisto & plotHistoSame
TCanvas *plotHistoAuto(TH1D *h1, TString canvName, int logX=0, int logY=0,
	       TString drawOpt="LPE", TString explain="", int gridLines=3);

TCanvas* plotHistoError(const TH1D* h1, TString cName, int logX=0, int logY=0,
		TString drawOpt="LPE", TString explain="", int gridLines=3);
TCanvas* plotHistoErrorSame(const TH1D *h1, TString canvName, TString drawOpt,
			    TString explain="");

TCanvas* plotGraphSame(TGraphErrors *h1, TString canvName, TString drawOpt, TString explain="");
TCanvas* plotRatio(TH1D* h1, TString cName, int logX=0, int logY=0,
		   TString drawOpt="LPE", TString explain="", int gridLines=3,
		   int lineAtOne=1);

int moveLegend(TCanvas *c, double dxNDC, double dyNDC);
int increaseLegend(TCanvas *c, double dxNDC_right, double dyNDC_top);
int changeLegendEntry(TCanvas *c, const TH1D *h1, TString newLabel);
void setLeftMargin(TCanvas *c, double xMargin);
void setRightMargin(TCanvas *c, double xMargin);
void setLeftRightMargins(TCanvas *c, double lMargin, double rMargin);

typedef enum { _massFrameCS=1, _massFrameYield } TMassFrameKind_t;
TCanvas* createMassFrame(int iFrame, TString canvNameBase, TString titleStr,
			 int yRangeSet=_massFrameCS,
			 TString *canvName_out=NULL, TH2D **h2frame_out=NULL);

void printHisto(const TH1D* h1, int extraRange=0);
void printHisto(const TH2D* h2, int extraRange=0, int nonZero=0,
		const TH2D* h2DenomHisto=NULL, double treshold=1e-3);
void printHistoWLimit(const TH2D* h2, int nLinesX, int nLinesY,
		      int nonZero=0, double treshold=1e-3);
void printHistoRange(const TH2D* h2);
void printRatio(const TH1D* h1a, const TH1D* h1b, int extraRange=0,
		int includeErr=0, double markIfDiff_RelTol=0.);
void printRatio(const TH2D* h2a, const TH2D* h2b, int extraRange=0,
		int includeErr=0, double markIfDiffRelTol=0.);
void printField(TString keyName);

inline void plotHisto(TH1* h1, TString cName, int logX=0, int logY=0, TString drawOpt="hist", TString explain="", int gridLines=3)
{  plotHisto((TH1D*)h1,cName,logX,logY,drawOpt,explain,gridLines); }
inline void plotHisto(TH2* h2, TString cName, int logX=0, int logY=0,
		      int gridLines=1)
{  plotHisto((TH2D*)h2,cName,logX,logY,1.5,gridLines,0.); }

inline void printHistoTH1(TH1* h1)
{  printHisto((TH1D*)h1); }
inline void printHistoTH2(TH2* h2)
{  printHisto((TH2D*)h2); }


template <class th1_t>
inline
void removeError(th1_t *h1)
{
  for (int ibin=1; ibin<=h1->GetNbinsX(); ibin++)
    h1->SetBinError(ibin,0);
}

template <class th2_t>
inline
void removeErrorH2(th2_t *h2)
{
  for (int ibin=1; ibin<=h2->GetNbinsX(); ibin++)
    for (int jbin=1; jbin<=h2->GetNbinsY(); jbin++)
      h2->SetBinError(ibin,jbin,0);
}

void nullifyOverflow(TH1D *h1);

// if h1ref, differences are taken wrt to 1st histo in the array
int deriveRelSyst(const TH1D *h1ref, const std::vector<TH1D*> &h1V,
		  std::vector<TH1D*> &h1diffV, int relative, int absValues);

TH1D* errorAsCentral(const TH1D* h1, int relative=0);
TH2D* errorAsCentral(const TH2D* h2, int relative=0);
TH1D* addUncertainty(const TString newName, const TH1D *h1_dest,
		     const TH1D *h1err, const TH1D *h1central=NULL);
int setError(TH1D *h1dest, const TH1D *h1src);
int setError(TH2D *h2dest, const TH2D *h2src);
int allGE(const TH1D *h1a, const TH1D *h1b); // 1 if h1a[i]>=h1b[i] for all i.
int assignDiffAsUnc(TH2D *h2target, const TH2D *h2nominal, int relative,
		    int listRangesOnError=0);
int addShiftByUnc(TH2D *h2target, const TH2D *h2src, double nSigmas,
		  int listRangesOnError=0);
int addInQuadrature(TH1D *h1target, const TH1D *h1addTerm, int nullifyErr=1);
int hasValueAbove(const TH2D *h2, double limit);
int hasValueBelow(const TH2D *h2, double limit);
void removeNegatives(TH1D* h1);
int compareRanges(const TH1D *h1a, const TH1D *h1b, int verbose=0);
int compareRanges(const TH2D *h2a, const TH2D *h2b, int verbose=0);
int checkRange(const TH1D* h1, double rangeMin, double rangeMax, int silent=0);
int checkRange(const std::vector<TH1D*> &h1V,
	       double &rangeMin, double &rangeMax,
	       const std::vector<std::pair<double,double> > &ranges,
	       int ignoreZeroValues=0, int silent=1);


void scaleBin(TH1D *h1, int ibin, double x);
void printBin(const TH1D *h1, int ibin, int newLine=1);
int hasDoubleZero(const TH2D *h2);
void setToOne(TH2D* h2);

// -----------------------------------------------------------

TH1D* perMassBinWidth(const TH1D* h1, int prnBinW=0);
TH2D* perMassBinWidth(const TH2D* h2, int prnBinW=0);
TH1D* timesMassBinWidth(const TH1D* h1, int prnBinW=0);
TH2D* timesMassBinWidth(const TH2D* h2, int prnBinW=0);
TH1D* flattenHisto(const TH2D *h2, TString setName);

TH1D* convert(TGraphAsymmErrors *gr, TString hName, TString hTitle,
	      int plotIt=0, int posNegErrs=0);

// -----------------------------------------------------------

void printObjStringField(TFile &f, TString keyName);

// -----------------------------------------------------------
// -----------------------------------------------------------

TCanvas *loadCanvas(TFile &fin, TString nameOnFile, TString newName="",
		    int warnIfMissing=1);
TCanvas *loadCanvas(TString fileName, TString nameOnFile, TString newName="",
		    int warnIfMissing=1);

// -----------------------------------------------------------

extern const TH1D* h1dummy;
extern const TH2D* h2dummy;


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
histo_t* loadHisto(TFile &fin, TString histoNameOnFile,
		   TString set_histoName, TString set_histoTitle,
		   int absenceIsError, const histo_t *dummy)
{
  histo_t *h=loadHisto(fin,histoNameOnFile,set_histoName,absenceIsError,dummy);
  if (h && set_histoTitle.Length()) h->SetTitle(set_histoTitle);
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
histo_t* loadHisto(TString fname, TString histoNameOnFile, TString histoName,
		   int absenceIsError, const histo_t *dummy)
{
  histo_t *h1=loadHisto(fname,histoNameOnFile,histoName,dummy);
  if (!h1 && absenceIsError) {
    std::cout << "failed to load histo with name=<" << histoNameOnFile << "> "
	      << "from file <" << fname << ">\n";
  }
  return h1;
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
  if (histoTitle.Index(";")==-1) {
    h->GetXaxis()->SetTitle(hSrc->GetXaxis()->GetTitle());
    h->GetYaxis()->SetTitle(hSrc->GetYaxis()->GetTitle());
  }
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

template<class histo_t>
inline
void histoStyle(histo_t* h1, int color, int markerStyle, int lineStyle=1,
		double markerSize=1.)
{
  h1->SetLineColor(color);
  h1->SetLineStyle(lineStyle);
  h1->SetMarkerColor(color);
  h1->SetMarkerStyle(markerStyle);
  h1->SetMarkerSize(markerSize);
}

// -----------------------------------------------------------

inline
void graphStyle(TGraphErrors *gr, int color, int markerStyle, int lineStyle=1,
		double markerSize=1.)
{ return histoStyle(gr,color,markerStyle,lineStyle,markerSize); }

// -----------------------------------------------------------

inline void histoStyle(TH1D *h1, const HistoStyle_t &hs) { hs.SetStyle(h1); }
inline void graphStyle(TGraphErrors *gr, const HistoStyle_t &hs)
{ hs.SetStyle(gr); }

// -----------------------------------------------------------

template<class graph1_t, class graph2_t>
inline
void copyStyle(graph1_t *grDest, const graph2_t *gr)
{
  grDest->SetLineColor(gr->GetLineColor());
  grDest->SetLineStyle(gr->GetLineStyle());
  grDest->SetLineWidth(gr->GetLineWidth());
  grDest->SetMarkerColor(gr->GetMarkerColor());
  grDest->SetMarkerStyle(gr->GetMarkerStyle());
  grDest->SetMarkerSize(gr->GetMarkerSize());
}

// -----------------------------------------------------------

template<class graph_t>
inline
void logAxis(graph_t *gr, int axis=1+2, TString xlabel="", TString ylabel="",
	     TString setTitle="")
{
  gr->GetXaxis()->SetDecimals(true);
  gr->GetYaxis()->SetDecimals(true);
  //gr->GetZaxis()->SetDecimals(true); // does not work for TGraphErrors
  if ((axis & 1)!=0) {
    gr->GetXaxis()->SetNoExponent();
    gr->GetXaxis()->SetMoreLogLabels();
  }
  if (((axis & 2)!=0) || ((axis & 4)!=0) || ((axis & 8)!=0)) {
    gr->GetYaxis()->SetNoExponent();
    if ((axis&8)!=0) gr->GetYaxis()->SetNoExponent(false);
    if ((axis&4)!=0) gr->GetYaxis()->SetMoreLogLabels();
  }
  if (xlabel.Length()) gr->GetXaxis()->SetTitle(xlabel);
  if (ylabel.Length()) gr->GetYaxis()->SetTitle(ylabel);
  if (setTitle.Length()) gr->SetTitle(setTitle);
}

// -----------------------------------------------------------

template<class graph_t>
inline
void ratioAxis(graph_t *gr, int axis=1+2, TString xlabel="", TString ylabel="")
{
  gr->GetXaxis()->SetDecimals(true);
  gr->GetYaxis()->SetDecimals(true);
  gr->GetYaxis()->SetTitleSize(0.07);
  gr->GetYaxis()->SetTitleOffset(0.69);
  gr->GetYaxis()->SetLabelSize(0.07);
  gr->GetXaxis()->SetTitleSize(0.07);
  gr->GetXaxis()->SetLabelSize(0.07);
  if ((axis & 1)!=0) {
    gr->GetXaxis()->SetNoExponent();
    gr->GetXaxis()->SetMoreLogLabels();
  }
  if ((axis & 2)!=0) {
    gr->GetYaxis()->SetNoExponent();
    //gr->GetYaxis()->SetMoreLogLabels();
  }
  if (xlabel.Length()) gr->GetXaxis()->SetTitle(xlabel);
  if (ylabel.Length()) gr->GetYaxis()->SetTitle(ylabel);
}

// -----------------------------------------------------------

inline
TString niceMassAxisLabel(int iLepton, TString extra, int iLabel=0)
{
  TString lepton;
  if (iLepton==0) lepton="e";
  else if (iLepton==1) lepton="\\mu";
  else lepton="\\ell\\!";
  TString sublabel=lepton+lepton;
  TString label=sublabel;
  if (extra.Length()) {
    sublabel+=",\\,";
    label= sublabel + TString("\\text{") + extra + TString("}");
  }
  TString start="M", units="[GeV]";
  if (iLabel==1) { start="\\sigma"; label.Prepend("\\!"); units="[pb/GeV]"; }
  return start + TString("_{") + label + TString("}\\text{ " + units + "}");
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

template<class histo2D_t>
inline
int copyContents(TH2D *h2Dest, const histo2D_t *h2Src)
{
  if ((h2Dest->GetNbinsX() != h2Src->GetNbinsX()) ||
      (h2Dest->GetNbinsY() != h2Src->GetNbinsY())) {
    std::cout << "copyContents(TH2D): number of bins is different :"
	      << " h2Dest[" << h2Dest->GetNbinsX() << "]["
	      << h2Dest->GetNbinsY() << "], "
	      << "h2Src[" << h2Src->GetNbinsX() << "]["
	      << h2Src->GetNbinsY() << "]\n";
    return 0;
  }
  for (int ibin=1; ibin<=h2Src->GetNbinsX(); ibin++) {
    for (int jbin=1; jbin<=h2Src->GetNbinsX(); jbin++) {
      h2Dest->SetBinContent(ibin,jbin, h2Src->GetBinContent(ibin,jbin));
      h2Dest->SetBinError  (ibin,jbin, h2Src->GetBinError  (ibin,jbin));
    }
  }
  return 1;
}

// -----------------------------------------------------------
// -----------------------------------------------------------

TH1D* loadVectorD(TString fname, TString valueField, TString errorField,
		  TString setName, TString setTitle, const TH1D *h1protoType);

TTree* loadTree(TString fname, TString treeName, TFile **fout_to_initialize);

// -----------------------------------------------------------
// -----------------------------------------------------------

void randomizeWithinErr(const TH1D *hSrc, TH1D *hDest, int nonNegative,
			int poissonRnd=0, int systematic=0);
void randomizeWithinRelErr(const TH1D *hSrc, TH1D *hDest, double relErr,
			   int nonNegative, int systematic=0);
void randomizeWithinErr(const TH2D *hSrc, TH2D *hDest, int nonNegative,
			int poissonRnd=0, int systematic=0);
void randomizeWithinRelErrSetCentral(const TH2D *h2Central,
				     const TH2D *h2RelUnc,
				     TH2D *hDest, int nonNegative,
				     int systematic=0);
int randomizeWithinPosNegErr(const TH2D *h2_errPos,
			     const TH2D *h2_errNeg,
			     TH2D *h2Dest, int bound01=0, int systematic=0);

// -----------------------------------------------------------
// -----------------------------------------------------------

void printMDim(const TString name, const TMatrixD &M);
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


template<class histo1D_t>
inline
TMatrixD convert2mat1D(const histo1D_t* h1)
{
  TMatrixD v(h1->GetNbinsX(),1);
  for (int ibin=1; ibin<=h1->GetNbinsX(); ibin++)
    v(ibin-1,0)= h1->GetBinContent(ibin);
  return v;
}

template<class histo1D_t>
inline
TMatrixD convert2cov1D(const histo1D_t* h1, int takeError=1)
{
  TMatrixD m(h1->GetNbinsX(),h1->GetNbinsX());
  for (int ibin=1; ibin<=h1->GetNbinsX(); ibin++) {
    double v= (takeError) ? h1->GetBinError(ibin) : h1->GetBinContent(ibin);
    m(ibin-1,ibin-1)= v*v;
  }
  return m;
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

TH2D* convert2histo(const TMatrixD &m, const TH1D *h1_for_axes,
		    TString h2name, TString h2title,
		    const TMatrixD *mErr=NULL);
TH1D* convert2histo1D(const TMatrixD &m, const TMatrixD &cov,
		      TString h1name, TString h1title,
		      const TH1D *h1_for_axes, int iCol=0);

// Change uncertainties, keeping correlations
int changeCov(TMatrixD &cov, const TH1D *h1unc);

// chi^2_estimate= (vec1-vec2)^T Mcov^{-1} (vec1-vec2)
double chi2estimate(const TVectorD &vec1, const TVectorD &vec2,
		    const TMatrixD &Mcov);

int deriveCovariance(const std::vector<TH1D*> &rndCS,
		     TString histoNameTag, TString histoTitle,
		     TH1D **h1avgCS_out, TH2D **h2cov_out);

//TH2D* convert(const TMatrixD &m, const
TH2D* cov2corr(const TH2D* h2cov);

// covariance of summed quantities, like sigmaZ=sum_a sigma_a
// idxMin..(idxMax_p1-1) correspond to lower and upper indices of the range
// if h1binWidth_def is provided, the covariance is assumed to be normalized
double sumCov(const TMatrixD &cov, int idxMin, int idxMax_p1,
	      const TH1D *h1binWidth_def=NULL, int printNumbers=0);

// uncertainty from covariance. If h1centralVal is supplied, the returned
// uncertainty is relative
TH1D* uncFromCov(const TH2D *h2cov, const TH1D *h1centralVal=NULL,
		 int zeroCentralMeansZeroRelError=0);
TH1D* uncFromCov(const TMatrixD &covM, TString hName,
		 const TH1D *h1_binning,
		 const TH1D *h1centralVal=NULL,
		 int zeroCentralMeansZeroRelError=0);

double totUnc(const std::vector<double> &errV, int returnSqrValue=0);

/*
TCanvas *plotCovCorr(TH2D* h2cov, TString canvName,
		     TH2D** h2corr_out=NULL,
		     int autoZRangeCorr=1,
		     int gridLines=1, int logScale=1,
		     double yTitleOffset=1.5);
TCanvas *plotCovCorr(const TMatrixD &cov, const TH1D *h1_for_axis_def,
		     TString histoName, TString canvName,
		     TH2D **h2cov_out=NULL, TH2D **h2corr_out=NULL,
		     int autoZRangeCorr=1,
		     int gridLines=1, int logScale=1,
		     double yTitleOffset=1.5);
*/

TCanvas *plotCovCorr(TH2D* h2cov, TString canvName,
		     const PlotCovCorrOpt_t o=PlotCovCorrOpt_t(),
		     TH2D** h2corr_out=NULL);
TCanvas *plotCovCorr(const TMatrixD &cov, const TH1D *h1_for_axis_def,
		     TString histoName, TString canvName,
		     const PlotCovCorrOpt_t o=PlotCovCorrOpt_t(),
		     TH2D **h2cov_out=NULL, TH2D **h2corr_out=NULL);

TCanvas *findCanvas(TString canvName);
int findCanvases(TString canvNames, std::vector<TCanvas*> &cV);
void closeCanvases(int itimes=3);

void SaveCanvas(TCanvas* canv, TString canvName, TString destDir,
		int correctName=0);
void SaveCanvases(std::vector<TCanvas*> &cV, TString destDir, TFile *fout=NULL);
// if listOfCanvNames="ALL", saves all canvases
int  SaveCanvases(TString listOfCanvNames, TString destDir, TFile *fout=NULL);
void writeTimeTag(TFile *fout=NULL);
void writeMsg(TString msg, TFile *fout=NULL);
int getHistosFromCanvas(TCanvas *c, std::vector<TH1D*> *h1V=NULL,
			std::vector<TH2D*> *h2V=NULL);

// -----------------------------------------------------------

inline
TCanvas* plotHisto(const TH1D* h1, TString cName, int logX=0, int logY=0,
		   TString drawOpt="LPE", TString explain="", int gridLines=3) {
  TH1D *h1tmp=cloneHisto(h1,h1->GetName() + TString("_plotTmp"),h1->GetTitle());
  return plotHisto(h1tmp,cName,logX,logY,drawOpt,explain,gridLines);
}

inline
TCanvas* plotHistoSame(const TH1D *h1, TString canvName, TString drawOpt,
		       TString explain="") {
  TH1D *h1tmp=cloneHisto(h1,h1->GetName() + TString("_plotTmp"),h1->GetTitle());
  return plotHistoSame(h1tmp,canvName,drawOpt,explain);
}

// -----------------------------------------------------------
// -----------------------------------------------------------

template<class T>
inline
void printVec(const char *msg, const std::vector<T> &vec)
{
  std::cout << msg << " [" << vec.size() << "]: ";
  for (unsigned int i=0; i<vec.size(); i++) std::cout << " " << vec[i];
  std::cout << "\n";
}

// -----------------------------------

template<class T>
inline
void printVec(const char *msg, const std::vector<T> &vec, const std::vector<int> &idx)
{
  std::cout << msg << " (sorted) [" << vec.size() << "]: ";
  for (unsigned int i=0; i<vec.size(); i++) std::cout << " " << vec[idx[i]];
  std::cout << "\n";
}

// -----------------------------------

std::vector<int> getSortDescendingIdx(const std::vector<double> &vec);
std::vector<int> getSortDescendingIdx(const std::vector<float> &vec);
int orderChanged(const std::vector<int> &idx);

// -----------------------------------------------------------
// -----------------------------------------------------------


template<class T>
inline
int valueEquals(const T val, const std::string targets)
{
  std::stringstream ss(targets);
  T probe;
  int yes=0;
  while (!ss.eof()) {
    ss >> probe;
    if (probe==val) { yes=1; break; }
  }
  return yes;
}

// -----------------------------------

template<class T>
inline
int hasValue(const T val, const std::string targets)
{ return valueEquals(val,targets); }

// -----------------------------------

template<>
inline
int hasValue(const char *val, const std::string targets)
{ return valueEquals<std::string>(std::string(val),targets); }

// -----------------------------------

// parse strings like "key1=value1 key2=value2"

template<class T>
inline
int setValue(T &var, const std::string &option, const std::string &key,
	     T default_value)
{
  var=default_value;
  size_t p=option.find(key);
  if (p!=std::string::npos) {
    p+= key.size();
    var=T(atof(option.c_str()+p));
    return 1;
  }
  return 0;
}

// -----------------------------------

inline
void HERE(const char *msg)
{ std::cout << msg << std::endl; }

template<class type_t>
inline
void HERE(const char *format, const type_t x)
{ std::cout << Form(format,x) << std::endl; }

template<class type1_t, class type2_t>
inline
void HERE(const char *format, const type1_t x, const type2_t y)
{ std::cout << Form(format,x,y) << std::endl; }


void ptrOk(const void *ptr);

inline
int PosOk(std::string s, std::string substr, size_t from_pos=0)
{
  size_t p=s.find(substr,from_pos);
  return (p!=std::string::npos) ? int(p) : -1;
}

// -----------------------------------------------------------
// -----------------------------------------------------------

std::ostream& operator<<(std::ostream &out, const TLorentzVector &v);
std::ostream& operator<<(std::ostream &out, const TLorentzVector *v);

// -----------------------------------------------------------
// -----------------------------------------------------------

const HistoStyle_t hsRed(kRed,  5,1,1.0,2.);
const HistoStyle_t hsBlue(kBlue,24,1,0.8,2.);
const HistoStyle_t hsGreen(kGreen+1,20,1,0.8,2.);
const HistoStyle_t hsBlack(kBlack,3,1,1.0,2.);
const HistoStyle_t hsColor6(6 ,   26,1,0.8,2.);
const HistoStyle_t hsColor46(46,   23,1,0.8,2.);
const HistoStyle_t hsBlue2(kBlue,24,1,0.8,2.);
const HistoStyle_t hsViolet(9,25,1,0.8,2.); // square
extern const HistoStyle_t hsVec[6];

// -----------------------------------------------------------

#endif
