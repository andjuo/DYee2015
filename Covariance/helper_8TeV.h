#ifndef helper_8TeV_H
#define helper_8TeV_H

#include <TROOT.h>
#include <TH2D.h>
#include <TString.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

// -------------------------------------------------------------

const int _nMassBins2012 = 41;
const double _massBinLimits2012[_nMassBins2012+1] =
  {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76,
   81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141,
   150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440,
   510, 600, 1000, 1500, 2000 }; // 41 bin

const int _nMassBins2D = 6;
const double _massBinLimits2D[_nMassBins2D+1] =
  {20, 30, 45, 60, 120, 200, 1500}; // 6 mass bins
const int _nFlatBins2D= (_nMassBins2D-1)*24 + 12;

const double _flatBins2D[_nFlatBins2D+1] = {
  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
  10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
  20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
  30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
  40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
  50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
  60, 61, 62, 63, 64, 65, 66, 67, 68, 69,
  70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
  80, 81, 82, 83, 84, 85, 86, 87, 88, 89,
  90, 91, 92, 93, 94, 95, 96, 97, 98, 99,
  100, 101, 102, 103, 104, 105, 106, 107, 108, 109,
  110, 111, 112, 113, 114, 115, 116, 117, 118, 119,
  120, 121, 122, 123, 124, 125, 126, 127, 128, 129,
  130, 131, 132 };

// -------------------------------------------------------------

typedef enum { _20to30=0, _30to45, _45to60, _60to120, _120to200, _200to1500,
    _last }
  TRange_t;

typedef enum { _noDiffUnc, _diffAbsUnc, _diffRelUnc,
	       _ratioUnc, } TUncertainty_t;

// -------------------------------------------------------------

const int nSystUnc=7;
const TString systUncNames[nSystUnc]= { "EScale", "EffRho", "DetRes", "BkgEst",
					"FSR", "CollCS", "Total" };

// -------------------------------------------------------------

inline
TString GetRangeName(TRange_t r) {
  TString name;
  switch(r) {
  case _20to30: name="20to30"; break;
  case _30to45: name="30to45"; break;
  case _45to60: name="45to60"; break;
  case _60to120: name="60to120"; break;
  case _120to200: name="120to200"; break;
  case _200to1500: name="200to1500"; break;
  case _last: name="last"; break;
  default: name="unknownRange";
  }
  return name;
}

// -------------------------------------------------------------

inline
TRange_t next(TRange_t &r) {
  if (r!=_last) r=TRange_t(int(r)+1);
  return r;
}

// -------------------------------------------------------------

inline
TString GetRangeAsString(TRange_t r, int caps=0) {
  TString name;
  switch(r) {
  case _20to30: name="20-30 GeV"; break;
  case _30to45: name="30-45 GeV"; break;
  case _45to60: name="45-60 GeV"; break;
  case _60to120: name="60-120 GeV"; break;
  case _120to200: name="120-200 GeV"; break;
  case _200to1500: name="200-1500 GeV"; break;
  case _last: name="last"; break;
  default: name="unknownRange";
  }
  if (caps) {
    name.ReplaceAll(" GeV",".0");
    name.ReplaceAll("-",".0 TO ");
  }
  return name;
}

// -------------------------------------------------------------

inline
void HERE(const TString &msg) { std::cout << msg << std::endl; }

// -------------------------------------------------------------

inline
TH2D* loadAlexCov1D(TString fName1, TString name, TString title,
		    int nBins=_nMassBins2012,
		    const double *binLimits=_massBinLimits2012)
{
  std::ifstream fin1(fName1.Data());
  if (!fin1.is_open()) {
    std::cout << "failed to open the file <" << fName1 << ">\n";
    return NULL;
  }

  TH2D *h2AlexCov=new TH2D(name,title,
			   nBins,binLimits,
			   nBins,binLimits);
  h2AlexCov->SetDirectory(0);
  h2AlexCov->SetStats(0);

  std::string s;
  for (int i=0; i<10; i++) {
    std::getline(fin1,s);
    //std::cout << "skip <" << s << ">\n";
  }

  double x;
  for (int i=0; i<_nMassBins2012; i++) {
    int iBin = _nMassBins2012 - i;
    for (int j=0; j<_nMassBins2012; j++) {
      int jBin = j+1;
      fin1 >> x;
      h2AlexCov->SetBinContent(iBin,jBin, x);
      h2AlexCov->SetBinError(iBin,jBin, 0.);
    }
  }
  fin1.close();
  return h2AlexCov;
}

// -------------------------------------------------------------

inline
TH2D* loadAlexCov2D(TString fName1, TString name, TString title,
		    int nBins=_nFlatBins2D,
		    int skipLines=1,
		    int formatIJVal=0)
{
  std::ifstream fin1(fName1.Data());
  if (!fin1.is_open()) {
    std::cout << "failed to open the file <" << fName1 << ">\n";
    return NULL;
  }

  TH2D *h2AlexCov=new TH2D(name,title,
			   nBins,0,nBins,
			   nBins,0,nBins);
  h2AlexCov->SetDirectory(0);
  h2AlexCov->SetStats(0);

  std::string s;
  if (skipLines) {
    for (int i=0; i<17; i++) {
      std::getline(fin1,s);
      std::cout << "skip <" << s << ">\n";
    }
  }

  double x;
  int idxI,idxJ;
  for (int i=0; i<nBins; i++) {
    int iBin = nBins - i;
    for (int j=0; j<nBins; j++) {
      int jBin = j+1;
      if (formatIJVal==0) {
	fin1 >> x;
	h2AlexCov->SetBinContent(iBin,jBin, x);
	h2AlexCov->SetBinError  (iBin,jBin, 0.);
      }
      else {
	fin1 >> idxI >> idxJ >> x;
	if (idxI==idxJ) std::cout << "idxI==idxJ=" << idxI << ", cov=" << x << "\n";
	h2AlexCov->SetBinContent(idxI,idxJ, x);
	h2AlexCov->SetBinError  (idxI,idxJ, 0.);
      }
    }
  }
  fin1.close();

  return h2AlexCov;
}

// -------------------------------------------------------------

inline
TH2D* loadBLUECov(TString fName1, TString name, TString title,
		  int nBins=_nMassBins2012,
		  const double *binLimits=_massBinLimits2012,
		  int skipLines=0)
{
  std::ifstream fin1(fName1.Data());
  if (!fin1.is_open()) {
    std::cout << "failed to open the file <" << fName1 << ">\n";
    return NULL;
  }

  TH2D *h2Cov=NULL;
  if (binLimits!=NULL) {
    h2Cov= new TH2D(name,title, nBins,binLimits, nBins,binLimits);
  }
  else h2Cov= new TH2D(name,title, nBins,0,nBins, nBins,0,nBins);
  h2Cov->SetDirectory(0);
  h2Cov->SetStats(0);

  // skip lines
  std::string line;
  for (int i=0; i<skipLines; i++) {
    getline(fin1,line);
    std::cout << "skipping line <" << line << ">\n";
  }

  int i,j;
  double x;
  for (int ii=0; ii<nBins; ii++) {
    for (int jj=0; jj<nBins; jj++) {
      fin1 >> i >> j >> x;
      if ((ii!=nBins-1) && (jj!=nBins-1) && fin1.eof()) {
	std::cout << "file <" << fName1 << "> might be too short\n";
	fin1.close();
	return NULL;
      }
      h2Cov->SetBinContent(i,j,x);
      h2Cov->SetBinError  (i,j,0.);
    }
  }

  fin1.close();
  return h2Cov;
}

// -------------------------------------------------------------

inline
int replaceErrorByCovUncertainty(TH1D *h, const TH2D *h2Cov,
				 int divideByBinWidth) {
  if ((h->GetNbinsX() != h2Cov->GetNbinsX()) ||
      (h->GetNbinsX() != h2Cov->GetNbinsY())) {
    std::cout << "ReplaceErrorByCovUncertainty: different sizes\n";
    return 0;
  }

  for (int ibin=1; ibin<=h->GetNbinsX(); ibin++) {
    double err= sqrt(h2Cov->GetBinContent(ibin,ibin));
    if (divideByBinWidth) err/=h->GetBinWidth(ibin+1);
    h->SetBinError(ibin, err);
  }
  return 1;
}

// -------------------------------------------------------------

inline
int replaceErrorByCovUncertainty2D(std::vector<TH1D*> &hV, const TH2D *h2Cov,
				   int divideByBinWidth) {
  int nBins=0;
  for (unsigned int i=0; i<hV.size(); i++) {
    nBins += hV[i]->GetNbinsX();
  }
  if ((nBins != h2Cov->GetNbinsX()) ||
      (nBins != h2Cov->GetNbinsY())) {
    std::cout << "ReplaceErrorByCovUncertainty2D: different sizes\n";
    return 0;
  }

  unsigned int ih=0;
  for (int icbin=1; icbin<=h2Cov->GetNbinsX(); ih++) {
    for (int ibin=1; ibin<=hV[ih]->GetNbinsX(); ++ibin, ++icbin) {
      double err=sqrt(h2Cov->GetBinContent(icbin,icbin));
      if (divideByBinWidth==1) {
	double factor= (icbin<=120) ? 0.1 : 0.2;
	err/=factor;
      }
      hV[ih]->SetBinError(ibin, err );
    }
  }
  return 1;
}

// -------------------------------------------------------------

TH1D* RemoveError(const TH1D* h) {
  TH1D *hOut = (TH1D*)h->Clone(h->GetName() + TString("noErr"));
  hOut->SetDirectory(0);
  hOut->SetStats(0);
  for (int ibin=1; ibin<=h->GetNbinsX(); ibin++) {
    hOut->SetBinError(ibin, 0.);
  }
  return hOut;
}

// -------------------------------------------------------------

inline
TH1D* GetUncertainty(const TH1D *h) {
  TH1D* hUnc=(TH1D*) h->Clone(h->GetName() + TString("_Unc"));
  hUnc->SetTitle(h->GetTitle() + TString("_Unc"));
  hUnc->SetDirectory(0);
  hUnc->SetStats(0);
  for (int ibin=1; ibin<=h->GetNbinsX(); ++ibin) {
    hUnc->SetBinContent(ibin, h->GetBinError(ibin));
    hUnc->SetBinError  (ibin, 0.);
  }
  return hUnc;
}
// -------------------------------------------------------------

inline
TH2D* cov2corr(const TH2D *h2cov)
{
  TString corrName=h2cov->GetName();
  corrName.ReplaceAll("Cov","Corr");
  corrName.ReplaceAll("cov","corr");
  TH2D* h2corr=(TH2D*)h2cov->Clone(corrName);
  h2corr->SetTitle(corrName);
  h2corr->SetDirectory(0);
  h2corr->SetStats(0);

  for (int ibin=1; ibin<=h2cov->GetNbinsX(); ++ibin) {
    for (int jbin=1; jbin<=h2cov->GetNbinsY(); ++jbin) {
	h2corr->SetBinContent(ibin,jbin,
			      h2cov->GetBinContent(ibin,jbin) /
			      (sqrt(h2cov->GetBinContent(ibin,ibin)) *
			       sqrt(h2cov->GetBinContent(jbin,jbin))) );
	if (fabs(h2corr->GetBinContent(ibin,jbin))>(1.+double(1e-12))) {
	  double rd=(fabs(h2corr->GetBinContent(ibin,jbin))-1.)*100.;
	  std::cout << "h2corr("  << h2corr->GetName() << ") "
		    << "> 1 at ibin=" << ibin << ", jbin=" << jbin
		    << "; cov=" << h2cov->GetBinContent(ibin,jbin)
		    << "; relDiff=" << rd << "%\n";
	}
	h2corr->SetBinError(ibin,jbin,0.);
    }
  }
  return h2corr;
}

// -------------------------------------------------------------

inline
TH1D* CombineHistos(const std::vector<TH1D*> &h1V,
		    TString hName="hComb", TString hTitle="hComb")
{
  int binCount=0;
  for (unsigned int i=0; i<h1V.size(); ++i) {
    binCount += h1V[i]->GetNbinsX();
  }
  TH1D* hComb= new TH1D(hName,hTitle,binCount,0,binCount);
  hComb->SetDirectory(0);

  unsigned int ih=0;
  for (int ibin=1; ibin<=hComb->GetNbinsX(); ih++) {
    for (int ihbin=1; ihbin<=h1V[ih]->GetNbinsX(); ++ihbin,++ibin) {
      hComb->SetBinContent(ibin, h1V[ih]->GetBinContent(ihbin));
      hComb->SetBinError  (ibin, h1V[ih]->GetBinError  (ihbin));
    }
  }
  return hComb;
}

// -------------------------------------------------------------

inline
int LoadXSHistos2D(TFile &f, std::vector<TH1D*> &h1V, TString histoTag="") {
  for (TRange_t r=_20to30; r<=_200to1500; next(r)) {
    TH1D* h= (TH1D*)f.Get(GetRangeName(r));
    if (!h) {
      std::cout << "failed to get histogram for r=" << GetRangeName(r) << "\n";
      return 0;
    }
    h->SetDirectory(0);
    h->SetStats(0);
    if (histoTag.Length()) h->SetName(h->GetName() + histoTag);
    h1V.push_back(h);
  }
  return 1;
}

// -------------------------------------------------------------

inline
TH1D* LoadAndCombineXSHistos2D(TFile &f, std::vector<TH1D*> &h1V,
			       TString hName="hComb", TString hTitle="hComb",
			       TString histoTag="")
{
  if (!LoadXSHistos2D(f,h1V,histoTag)) {
    std::cout << "LoadAndCombineXSHistos2D failed to load" << std::endl;
    return NULL;
  }
  TH1D *hComb= CombineHistos(h1V,hName,hTitle);
  return hComb;
}

// -------------------------------------------------------------

inline
TString formNumber(double val0, double err0, int verbatim=0)
{
  TString str;
  double factor=1.;
  double val=val0, err=err0;
  int dfCount=0;

  while (err>100.) { factor*=0.1; dfCount--; err*=0.1; }
  while (err<10.) { factor*=10.; dfCount++; err*=10.; }


  val*=factor;

  val=round(val);
  err=round(err);
  TString format;
  //std::cout << "dfCount=" << dfCount << "\n";
  if (dfCount>=0) {
    format=Form("%c2.%dlf +- %c2.%dlf",'%',dfCount,'%',dfCount);
    str=Form(format.Data(),val/factor,err/factor);
  }
  else {
    int correction=0;
    if (correction) {
      val/=10.;
      err/=10.;
    }
    format=Form("%c2.%dlfE+%d +- %c2.%dlfE+%d",
		'%',correction,-dfCount+correction,
		'%',correction,-dfCount+correction);
    //std::cout << "format=<" << format << ">\n";
    str=Form(format.Data(),val,err);
  }
  
  if (verbatim) {
    std::cout << "formNumber: ini=" << val0 << " +- " << err0 << ", "
	      << "corr=" << val << " +- " << err
	      << " (factor=" << factor << ")\n";
    std::cout << " - str= " << str << "\n";
  }

  return str;
}

// -------------------------------------------------------------

inline
void printHisto(const TH1D* h)
{
  std::cout << "histo " << h->GetName() << ", " << h->GetTitle() << "\n";
  for (int ibin=1; ibin<=h->GetNbinsX(); ++ibin) {
    std::cout << "ibin=" << ibin << " value=" << h->GetBinContent(ibin)
	      << " +- " << h->GetBinError(ibin) << "\n";
  }
}

// -------------------------------------------------------------

inline
void compareHistos(const TH1D* h1, const TH1D *h2)
{
  std::cout << "compare: \n";
  std::cout << " histo1 " << h1->GetName() << ", " << h1->GetTitle() << "\n";
  std::cout << " histo2 " << h2->GetName() << ", " << h2->GetTitle() << "\n";
  if (h1->GetNbinsX()!=h2->GetNbinsX()) {
    std::cout << "Different number of bins: " <<
      h1->GetNbinsX() << " vs " << h2->GetNbinsX() << "\n";
    return;
  }
  for (int ibin=1; ibin<=h1->GetNbinsX(); ++ibin) {
    std::cout << "ibin=" << ibin << " value=" << h1->GetBinContent(ibin)
	      << " +- " << h1->GetBinError(ibin) << "  ";
    std::cout << h2->GetBinContent(ibin) << " +- " << h2->GetBinError(ibin);
    double dv= h1->GetBinContent(ibin) - h2->GetBinContent(ibin);
    double de= h1->GetBinError(ibin) - h2->GetBinError(ibin);
    std::cout << "  diff=" << dv << " " << de << "\n";
  }
}

// -------------------------------------------------------------
// -----------------------------------------------------------------

inline
TH1D *loadASfile_asHisto1D(TString fname, TString label, int skipFirstLine=1) {
  std::ifstream fin(fname.Data());
  if (!fin.is_open()) {
    std::cout << "failed to open <" << fname
	      << "> for label=<" << label <<">\n";
    return NULL;
  }

  while (skipFirstLine>0) {
    std::string s;
    getline(fin,s);
    skipFirstLine--;
  }

  int ibin;
  double massMin, massMax;
  double val,err;

  TH1D *h1=new TH1D(label,label, _nMassBins2012,_massBinLimits2012);
  h1->SetDirectory(0);

  while (!fin.eof()) {
    fin >> ibin >> massMin >> massMax >> val >> err;
    if (fin.eof()) break;
    h1->SetBinContent(ibin,val);
    h1->SetBinError(ibin,err);
  }
  fin.close();

  return h1;
}

// -----------------------------------------------------------------
// -----------------------------------------------------------------

inline
TH1D *loadASfile_asOneHisto2D(TString fname,
			      TString label,
			      int skipFirstLine=1,
			      int containsMassBinVals=0
			      ) {
  std::ifstream fin(fname.Data());
  if (!fin.is_open()) {
    std::cout << "failed to open <" << fname
	      << "> for label=<" << label <<">\n";
    return NULL;
  }

  while (skipFirstLine>0) {
    std::string s;
    getline(fin,s);
    skipFirstLine--;
  }

  int ibin;
  double x,y,val,err;

  TH1D *h1=new TH1D(label,label, _nFlatBins2D, 0, _nFlatBins2D);
  h1->SetDirectory(0);

  ibin=0;
  while (!fin.eof()) {
    if (containsMassBinVals) {
      double im,iy,mMin,mMax,yMin,yMax;
      fin >> im >> iy >> mMin >> mMax >> yMin >> yMax >> val >> err;
      ibin++;
      val /= (yMax-yMin);
      err /= (yMax-yMin);
    }
    else {
      fin >> ibin >> x >> y >> val >> err;
      std::cout << "ibin=" << ibin << ", val=" << val << " , err=" << err << "\n";
    }
    if (fin.eof()) break;
    h1->SetBinContent(ibin,val);
    h1->SetBinError(ibin,err);
  }
  fin.close();

  return h1;
}

// -----------------------------------------------------------------

inline
int loadASfile2D_as1DHistoVec(TString fname,
			      TString label,
			      std::vector<TH1D*> &hV,
			      int skipFirstLine=0
				) {
  std::ifstream fin(fname.Data());
  if (!fin.is_open()) {
    std::cout << "failed to open <" << fname
	      << "> for label=<" << label <<">\n";
    return 0;
  }

  while (skipFirstLine>0) {
    std::string s;
    getline(fin,s);
    skipFirstLine--;
  }

  int ibin;
  double val,err;
  double itmp0,itmp1,itmp2;

  TH1D *h1=NULL;

  ibin=0;
  while (!fin.eof()) {
    if (ibin%24==0) {
      TString hname=Form("h1_iMass%d",ibin/24+1);
      if (label.Length()) { hname.Append("_"); hname.Append(label); }
      int cnt= (ibin/24==5) ? 12 : 24;
      h1=new TH1D(hname,label, cnt, 0., 2.4);
      h1->SetDirectory(0);
      h1->Reset();
      hV.push_back(h1);
    }

    fin >> itmp0;
    if (fin.eof()) break;
    fin >> itmp1 >> itmp2 >> val >> err;
    //std::cout << "itmp0=" << itmp0 << ", itmp1=" << itmp1 << ", itmp2="
    //	      << itmp2 << ", val=" << val << ", err=" << err << "\n";
    h1->SetBinContent(ibin%24 + 1, val);
    h1->SetBinError  (ibin%24 + 1, err);
    ibin++;
    if (itmp0!=ibin) {
      std::cout << "detected file expectation mismatch: "
		<< "itmp0=" << itmp0 << ", ibin=" << ibin << "\n";
      std::cout << " file name <" << fname << ">\n";
      return 0;
    }
  }
  fin.close();

  return 1;
}

// -----------------------------------------------------------------

//inline
//int PosOk(std::string s, std::string subs, size_t pos=0)
//{ return (s.find(subs,pos)==std::string::npos) ? 0 : 1; }

// -----------------------------------------------------------------


#endif
