#include "inputs.h"
#include <sstream>
#include <algorithm>

// --------------------------------------------------------------

const TH1D *h1dummy=NULL;
const TH2D *h2dummy=NULL;

// --------------------------------------------------------------

TString versionName(TVersion_t ver)
{
  TString name="versionNameUNKNOWN";
  switch(ver) {
  case _verUndef: name="UndefVer"; break;
  case _verMu1: name="Mu1"; break;
  case _verMu76X: name="Mu76X"; break;
  case _verMuApproved: name="MuApproved"; break;
  case _verEl1: name="El1"; break;
  case _verEl2: name="El2"; break;
  case _verEl2skim: name="El2skim"; break;
  case _verEl2skim2: name="El2skim2"; break;
  case _verEl2skim3: name="El2skim3"; break;
  case _verEl3: name="El3"; break;
  default:
    std::cout << "versionName is not ready for this version type\n";
  }
  return name;
}

// -----------------------------------------------------------
// -----------------------------------------------------------

TString DayAndTimeTag(int eliminateSigns)
{
   time_t ltime;
   ltime=time(NULL);
   TString str = TString(asctime(localtime(&ltime)));
   if (str[str.Length()-1]=='\n') str.Remove(str.Length()-1,1);
   if (eliminateSigns) {
     str.ReplaceAll(" ","_");
     str.ReplaceAll(":","");
   }
   return str;
}

// -----------------------------------------------------------

void addToVector(std::vector<TString> &vec, TString strings)
{
  std::stringstream ss(strings.Data());
  TString s;
  while (!ss.eof()) {
    ss >> s;
    if (s.Length()) vec.push_back(s);
  }
}

// -----------------------------------------------------------

void addToVector(std::vector<int> &vec, int count, ...)
{
  va_list vl;
  va_start(vl,count);
  for (int i=0; i<count; i++) {
    vec.push_back(va_arg(vl,int));
  }
  va_end(vl);
}

// -----------------------------------------------------------
// -----------------------------------------------------------

TCanvas* plotHisto(TH1D* h1, TString cName, int logX, int logY, TString drawOpt,
		   TString explain)
{
  TCanvas *cx= new TCanvas(cName,cName,600,600);
  if (logX) cx->SetLogx();
  if (logY) cx->SetLogy();
  h1->Draw(drawOpt);
  if (explain.Length()) {
    const int smaller=1;
    TLegend *leg= new TLegend(0.15,0.15,0.40,0.22-0.02*smaller);
    leg->SetName("myLegend");
    leg->AddEntry(h1,explain);
    leg->Draw();
  }
  cx->Update();
  return cx;
}

// ---------------------------------------------------------

TCanvas* plotHisto(TH2D* h2, TString cName, int logx, int logy,
		   double ytitleOffset, double centralZrange)
{
  TCanvas *cx= new TCanvas(cName,cName,600,600);
  cx->SetLeftMargin(0.15);
  cx->SetRightMargin(0.18);
  cx->SetLogx(logx);
  cx->SetLogy(logy);
  h2->GetYaxis()->SetTitleOffset(ytitleOffset);
  if (centralZrange!=0) {
    h2->GetZaxis()->SetRangeUser(-centralZrange,centralZrange);
    h2->GetZaxis()->SetDecimals(true);
  }
  h2->Draw("COLZ");
  cx->Update();
  return cx;
}

// ---------------------------------------------------------

TCanvas* plotHisto(const TH2D* h2_orig, TString cName, int logx, int logy,
		   double ytitleOffset, double centralZrange)
{
  TH2D *h2=cloneHisto(h2_orig,h2_orig->GetName() + TString("__clone_plotHisto"),
		      h2_orig->GetTitle());
  return plotHisto(h2,cName,logx,logy,ytitleOffset,centralZrange);
}

// ---------------------------------------------------------

TCanvas* plotHistoSame(TH1D *h1, TString canvName, TString drawOpt, TString explain)
{
  TSeqCollection *seq=gROOT->GetListOfCanvases();
  TCanvas *cOut=NULL;
  for (int i=0; i<=seq->LastIndex(); i++) {
    TCanvas *c= (TCanvas*) seq->At(i);
    //std::cout << "i=" << i << " " << c->GetName() << "\n";
    if (c->GetName() == canvName) {
      c->cd();
      h1->Draw(drawOpt + TString(" same"));

      if (explain.Length()) {
	TLegend *leg= (TLegend*)gPad->FindObject("myLegend");
	if (leg) {
	  const int smaller=1;
	  leg->SetY2NDC(leg->GetY2NDC() + 0.07-0.02*smaller);
	  leg->AddEntry(h1,explain);
	  leg->Draw();
	}
      }

      cOut=c;
    }
  }
  return cOut;
}

// ---------------------------------------------------------

int moveLegend(TCanvas *c, double dx, double dy)
{
  c->cd();
  TLegend *leg= (TLegend*)gPad->FindObject("myLegend");
  if (!leg) return 0;
  if (dx!=0.) {
    leg->SetX1NDC(leg->GetX1NDC() + dx);
    leg->SetX2NDC(leg->GetX2NDC() + dx);
  }
  if (dy!=0.) {
    leg->SetY1NDC(leg->GetY1NDC() + dy);
    leg->SetY2NDC(leg->GetY2NDC() + dy);
  }
  return 1;
}

// ---------------------------------------------------------

void printHisto(const TH1D *h1, int extraRange) {
  std::cout << "\nhisto " << h1->GetName() << " " << h1->GetTitle() << "\n";
  int d=(extraRange) ? 1 : 0;
  for (int ibin=1-d; ibin<=h1->GetNbinsX()+d; ibin++) {
    std::cout << "ibin=" << ibin << " " << h1->GetBinLowEdge(ibin)
	      << " " << (h1->GetBinLowEdge(ibin)+h1->GetBinWidth(ibin))
	      << "  " << h1->GetBinContent(ibin) << " +- "
	      << h1->GetBinError(ibin) << "\n";
  }
}

// ---------------------------------------------------------

void printHisto(const TH2D *h2, int extraRange, int nonZero,
		const TH2D *h2DenomHisto, double threshold) {
  int d=(extraRange) ? 1 : 0;
  std::cout << "\nhisto " << h2->GetName() << " " << h2->GetTitle() << "\n";
  if (nonZero) std::cout << "  (only non-zero values printed)\n";
  for (int ibin=1-d; ibin<=h2->GetNbinsX()+d; ibin++) {
    for (int jbin=1-d; jbin<=h2->GetNbinsY()+d; jbin++) {
      if (nonZero) {
	double currVal=h2->GetBinContent(ibin,jbin);
	int check=0;
	if (!h2DenomHisto) check=1;
	else {
	  double denom=h2DenomHisto->GetBinContent(ibin,jbin);
	  if (denom==0) check=1;
	  else { if (fabs(currVal/denom)<threshold) continue; }
	}
	if (check && (fabs(h2->GetBinContent(ibin,jbin))<1e-8)) continue;
      }
      std::cout << "ibin=" << ibin << ",jbin=" << jbin
		<< "  " << h2->GetXaxis()->GetBinLowEdge(ibin)
		<< "-" << (h2->GetXaxis()->GetBinLowEdge(ibin)+
			   h2->GetXaxis()->GetBinWidth(ibin))
		<< "  " << h2->GetYaxis()->GetBinLowEdge(jbin)
		<< "-" << (h2->GetYaxis()->GetBinLowEdge(jbin)+
			   h2->GetYaxis()->GetBinWidth(jbin))
		<< "  " << h2->GetBinContent(ibin,jbin) << " +- "
		<< h2->GetBinError(ibin,jbin) << "\n";
    }
  }
}

// ---------------------------------------------------------

void printRatio(const TH1D *h1a, const TH1D *h1b, int extraRange,
		int includeErr) {
  int d=(extraRange) ? 1:0;
  std::cout << "\nhisto A " << h1a->GetName() << " " << h1a->GetTitle()
	    << " [" << h1a->GetNbinsX() << "]\n";
  std::cout << "histo B " << h1b->GetName() << " " << h1b->GetTitle()
	    << " [" << h1b->GetNbinsX() << "]\n";
  std::cout << std::string(5,' ') << h1a->GetName() << std::string(10,' ')
	    << h1b->GetName() << std::string(10,' ') << "ratio\n";
  for (int ibin=1-d; ibin<=h1a->GetNbinsX()+d; ibin++) {
    if ( (h1a->GetBinLowEdge(ibin) != h1b->GetBinLowEdge(ibin)) ||
	 (h1a->GetBinWidth(ibin) != h1b->GetBinWidth(ibin)) ) {
      std::cout << "bining mismatch at ibin=" << ibin << "\n";
      return;
    }
    std::cout << "ibin=" << ibin << " " << h1a->GetBinLowEdge(ibin)
	      << " " << (h1a->GetBinLowEdge(ibin)+h1a->GetBinWidth(ibin))
	      << "  " << h1a->GetBinContent(ibin) << " +- "
	      << h1a->GetBinError(ibin)
	      << "  " << h1b->GetBinContent(ibin) << " +- "
	      << h1b->GetBinError(ibin)
	      << "  " << h1a->GetBinContent(ibin)/h1b->GetBinContent(ibin);
    if (includeErr)
      std::cout << " " << h1a->GetBinError(ibin)/h1b->GetBinError(ibin);
    std::cout << "\n";
  }
}


// ---------------------------------------------------------

void printRatio(const TH2D *h2a, const TH2D *h2b, int extraRange,
		int includeErr) {
  int d=(extraRange) ? 1:0;
  std::cout << "\nhisto A " << h2a->GetName() << " " << h2a->GetTitle()
	    << " [" << h2a->GetNbinsX() << "]\n";
  std::cout << "histo B " << h2b->GetName() << " " << h2b->GetTitle()
	    << " [" << h2b->GetNbinsX() << "]\n";
  std::cout << std::string(5,' ') << h2a->GetName() << std::string(10,' ')
	    << h2b->GetName() << std::string(10,' ') << "ratio\n";
  for (int ibin=1-d; ibin<=h2a->GetNbinsX()+d; ibin++) {
    if ( (h2a->GetXaxis()->GetBinLowEdge(ibin) !=
	  h2b->GetXaxis()->GetBinLowEdge(ibin)) ||
	 (h2a->GetXaxis()->GetBinWidth(ibin) !=
	  h2b->GetXaxis()->GetBinWidth(ibin)) ) {
      std::cout << "x-bining mismatch at ibin=" << ibin << "\n";
      return;
    }
    for (int jbin=1-d; jbin<=h2a->GetNbinsY()+d; jbin++) {
      if ( (h2a->GetYaxis()->GetBinLowEdge(jbin) !=
	    h2b->GetYaxis()->GetBinLowEdge(jbin)) ||
	   (h2a->GetYaxis()->GetBinWidth(jbin) !=
	    h2b->GetYaxis()->GetBinWidth(jbin)) ) {
	std::cout << "y-bining mismatch at jbin=" << jbin << "\n";
	return;
      }

      std::cout << "xybin=(" << ibin << "," << jbin << ")"
		<< " (" << h2a->GetXaxis()->GetBinLowEdge(ibin)
		<< " " << (h2a->GetXaxis()->GetBinLowEdge(ibin)+
			   h2a->GetXaxis()->GetBinWidth(ibin))
		<< " ; " << h2a->GetYaxis()->GetBinLowEdge(jbin)
		<< " " << (h2a->GetYaxis()->GetBinLowEdge(jbin)+
			   h2a->GetYaxis()->GetBinWidth(jbin))
		<< ")  " << h2a->GetBinContent(ibin,jbin) << " +- "
		<< h2a->GetBinError(ibin,jbin)
		<< "  " << h2b->GetBinContent(ibin,jbin) << " +- "
		<< h2b->GetBinError(ibin,jbin)
		<< "  " << h2a->GetBinContent(ibin,jbin)/h2b->GetBinContent(ibin,jbin);
      if (includeErr) {
	std::cout << " " << h2a->GetBinError(ibin,jbin)/h2b->GetBinError(ibin,jbin);
      }
      std::cout << "\n";
    }
  }
}


// ---------------------------------------------------------

void printField(TString keyName)
{
  TObject *obj= gDirectory->Get(keyName);
  if (!obj) { std::cout << "printField: obj <" << keyName << "> not found\n"; }
  else {
    TObjString *str=(TObjString*)obj;
    if (str) std::cout << "value: " << str->String() << "\n";
  }
}

// ---------------------------------------------------------

TH1D* errorAsCentral(const TH1D* h1, int relative) {
  TH1D* h1err= cloneHisto(h1,
			  h1->GetName() + TString("_err"),
			  h1->GetTitle() + TString(" err"));
  for (int ibin=1; ibin<=h1->GetNbinsX(); ibin++) {
    double err= h1->GetBinError(ibin);
    if (relative) {
      double val=h1->GetBinContent(ibin);
      if (val==double(0)) {
	if (err==double(0)) err=0;
	else if (err<double(0)) err=-1e4;
	else err=1e4;
      }
      else err/=val;
    }
    h1err->SetBinContent(ibin, err);
    h1err->SetBinError (ibin, 0.);
  }
  return h1err;
}

// ---------------------------------------------------------

void removeNegatives(TH1D* h1)
{
  for (int ibin=1; ibin<=h1->GetNbinsX(); ibin++) {
    if (h1->GetBinContent(ibin)<0) h1->SetBinContent(ibin,0.);
  }
}

// ---------------------------------------------------------

void scaleBin(TH1D* h1, int ibin, double x)
{
  h1->SetBinContent(ibin, x * h1->GetBinContent(ibin));
  h1->SetBinError  (ibin, x * h1->GetBinError(ibin));
}

// ---------------------------------------------------------

void printBin(TH1D *h1, int ibin, int newLine)
{
  std::cout << h1->GetName() << "[" << ibin << "]="
	    << h1->GetBinContent(ibin) << " +- "
	    << h1->GetBinError(ibin);
  if (newLine) std::cout << "\n";
}

// ---------------------------------------------------------
// ---------------------------------------------------------

TH1D* perMassBinWidth(const TH1D* h1_orig, int prnBinW)
{
  TH1D* h1=(TH1D*)h1_orig->Clone(h1_orig->GetName() + TString("_perMassBinW"));
  for (int ibin=1; ibin<=h1->GetNbinsX(); ibin++) {
    double w= h1->GetBinWidth(ibin);
    if (prnBinW) std::cout << "ibin=" << ibin << ", w=" << w << "\n";
    h1->SetBinContent( ibin, h1_orig->GetBinContent(ibin)/w );
    h1->SetBinError  ( ibin, h1_orig->GetBinError(ibin)/w   );
  }
  return h1;
}

// ---------------------------------------------------------

TH1D* flattenHisto(const TH2D *h2, TString setName)
{
  if (h2->GetNbinsY()==1) {
    const TArrayD* xb=h2->GetXaxis()->GetXbins();
    const Double_t *x= xb->GetArray();
    TH1D *h1= new TH1D(setName,setName,xb->GetSize()-1,x);
    h1->SetDirectory(0);
    h1->Sumw2();
    //delete x;

    for (int ibin=1; ibin<=h2->GetNbinsX(); ibin++) {
      h1->SetBinContent(ibin, h2->GetBinContent(ibin,1));
      h1->SetBinError  (ibin, h2->GetBinError(ibin,1));
    }

    return h1;
  }
  return NULL;
}

// ---------------------------------------------------------

TH1D* convert(TGraphAsymmErrors *gr, TString hName, TString hTitle,
	      int plotIt, int posNegErrs)
{
  Double_t *xLow= gr->GetEXlow();
  Double_t *xHigh= gr->GetEXhigh();
  Double_t *yLow= gr->GetEYlow();
  Double_t *yHigh= gr->GetEYhigh();
  int n=gr->GetN();
  Double_t *x= gr->GetX();
  Double_t *y= gr->GetY();

  Double_t *xb= new Double_t[n+1];

  //std::cout << "n=" << n << "\n";
  double factor=1.;
  //factor=1.05; // to see the difference
  for (int i=0; i<n; i++) {
    xb[i]=x[i]-xLow[i]*factor;
    //std::cout << "i=" << i << ", xLow=" << xLow[i] << ", xHigh=" << xHigh[i] << ", x=" << x[i] << "\n";
  }
  xb[n]= x[n-1]+xHigh[n-1];
  //std::cout << "last point " << xb[n] << "\n";


  TH1D *h1= new TH1D(hName, hTitle, n, xb);
  h1->SetDirectory(0);
  for (int ibin=1; ibin<=h1->GetNbinsX(); ibin++) {
    h1->SetBinContent(ibin, y[ibin-1]);
    double w=0.5;
    if (posNegErrs==-1) w=1.;
    else if (posNegErrs==1) w=0;
    h1->SetBinError  (ibin, w*yLow[ibin-1] + (1.-w)*yHigh[ibin-1]);
  }

  if (plotIt) {
    TString cName= "c" + hName;
    TCanvas *c= new TCanvas(cName,cName,600,600);
    c->SetLogx();
    h1->SetLineColor(kRed);
    h1->Draw("LPE1");
    gr->Draw("LPE same");
    c->Update();
  }

  return h1;
}

// ---------------------------------------------------------
// ---------------------------------------------------------

void printObjStringField(TFile &f, TString keyName) {
  TObjString *s=(TObjString*)f.Get(keyName);
  if (!s) {
    std::cout << "keyName=<" << keyName << ">, was not found in <"
	      << f.GetName() << ">\n";
    return;
  }
  std::cout << "keyName=" << keyName << ", value=<" << s->String() << ">\n";
  return;
}


// ---------------------------------------------------------
// ---------------------------------------------------------

TH1D* loadVectorD(TString fname, TString valueField, TString errorField,
		  TString setName, TString setTitle, const TH1D *h1protoType)
{
  TFile fin(fname);
  if (!fin.IsOpen()) {
    std::cout << "failed to open the file <" << fname << ">\n";
    return NULL;
  }
  TVectorD *val=(TVectorD*)fin.Get(valueField);
  TVectorD *err=(TVectorD*)fin.Get(errorField);
  fin.Close();
  if (!val || !err) {
    if (!val) std::cout << "failed to get value from <" << valueField << ">\n";
    if (!err) std::cout << "failed to get error from <" << errorField << ">\n";
    std::cout << " .. in file <" << fname << ">\n";
    return NULL;
  }

  if ((val->GetNoElements() != h1protoType->GetNbinsX()) ||
      (err->GetNoElements() != h1protoType->GetNbinsX())) {
    std::cout << "loadVectorD: number of elements is different: \n";
    std::cout << " val[" << val->GetNoElements() << "], "
	      << " err[" << err->GetNoElements() << "], "
	      << " h1protoType[" << h1protoType->GetNbinsX() << "]\n";
    return NULL;
  }

  TH1D *h1= (TH1D*)h1protoType->Clone(setName);
  h1->SetDirectory(0);
  h1->SetTitle(setTitle);

  for (int i=0; i<val->GetNoElements(); i++) {
    h1->SetBinContent(i+1, (*val)[i]);
    h1->SetBinError  (i+1, (*err)[i]);
  }
  return h1;
}

// ---------------------------------------------------------

TTree* loadTree(TString fname, TString treeName, TFile **fout_ptr)
{
  TFile *fout = new TFile(fname);
  if (fout_ptr) *fout_ptr = fout;
  if (!fout || !fout->IsOpen()) {
    std::cout << "failed to open the file <" << fname << ">\n";
    return NULL;
  }
  TTree *tree = (TTree*)fout->Get(treeName);
  if (!tree) {
    std::cout << "failed to get the tree " << treeName << " from a file <"
	      << fname << ">\n";
  }
  return tree;
}

// ---------------------------------------------------------
// ---------------------------------------------------------

TMatrixD submatrix(const TMatrixD &M, int idxMin, int idxMax)
{
  TMatrixD mNew(idxMax-idxMin,idxMax-idxMin);
  for (int i=0; i<idxMax-idxMin; i++) {
    for (int j=0; j<idxMax-idxMin; j++) {
      mNew(i,j) = M(i+idxMin,j+idxMin);
    }
  }
  return mNew;
}

// ---------------------------------------------------------

void printRatio(TString label1, const TVectorD &v1, TString label2, const TVectorD &v2)
{
  std::cout << "\nprintRatio " << std::string(5,' ')
	    << label1 << Form(" [%d]",v1.GetNoElements())
	    << std::string(5,' ')
	    << label2 << Form(" [%d]",v2.GetNoElements())
	    << std::string(5,' ') << "\n";
  if (v1.GetNoElements() != v2.GetNoElements()) {
    std::cout << "as the number of elements is different, the comparison is meaningless\n";
    return;
  }

  for (int i=0; i<v1.GetNoElements(); i++) {
    double r=v1[i]/v2[i];
    std::cout << Form("%3d",i) << "  " << v1[i] << "  " << v2[i]
	      << "  " << r << "\n";
  }
  return;
}

// ---------------------------------------------------------

TH2D* convert2histo(const TMatrixD &m, const TH1D *h1_for_axes,
		    TString h2name, TString h2title,
		    const TMatrixD *mErr)
{
  if (m.GetNrows()!=m.GetNcols()) {
    std::cout << "convert2histo: matrix is not square\n";
    std::cout << "  dims: " << m.GetNrows() << " x " << m.GetNcols() << "\n";
    return NULL;
  }
  if (m.GetNrows() != h1_for_axes->GetNbinsX()) {
    std::cout << "convert2histo: matrix has different number of elements "
	      << "than histogram provides bins\n";
    std::cout << "m[" << m.GetNrows() << "," << m.GetNcols() << "], histo["
	      << h1_for_axes->GetNbinsX() << "\n";
    return NULL;
  }
  if (mErr) {
    if (mErr->GetNrows()!=mErr->GetNcols()) {
      std::cout << "convert2histo: mErr is provided, and the matrix is not square\n";
      std::cout << "  dims: " << mErr->GetNrows() << " x "
		<< mErr->GetNcols() << "\n";
      return NULL;
    }
    if (mErr->GetNrows()!=m.GetNrows()) {
      std::cout << "convert2histo: mErr is provided, and the matrix sizes are different\n";
      std::cout << " dims: m[" << m.GetNrows() << "," << m.GetNcols() << "], "
		<< "mErr[" << mErr->GetNrows() << "," << mErr->GetNcols() << "]\n";
      return NULL;
    }
  }

  const TArrayD *xb= h1_for_axes->GetXaxis()->GetXbins();
  if (xb->GetSize()==0) {
    std::cout << "convert2histo: histogram has fixed axis. This case"
	      << " is not implemented\n";
    return NULL;
  }
  const Double_t *x= xb->GetArray();
  TH2D *h2= new TH2D(h2name,h2title, xb->GetSize()-1,x, xb->GetSize()-1,x);
  h2->SetDirectory(0);
  h2->SetStats(0);
  h2->Sumw2();
  if (h2title.Index(";")==-1) {
    TString label=h1_for_axes->GetXaxis()->GetTitle();
    h2->GetXaxis()->SetTitle(label);
    h2->GetYaxis()->SetTitle(label);
  }

  for (int ibin=1; ibin<xb->GetSize(); ibin++) {
    for (int jbin=1; jbin<xb->GetSize(); jbin++) {
      h2->SetBinContent( ibin, jbin, m(ibin-1,jbin-1) );
      double err= (mErr) ? (*mErr)(ibin-1,jbin-1) : 0.;
      h2->SetBinError( ibin, jbin, err );
    }
  }
  return h2;
}

// ---------------------------------------------------------

// chi^2_estimate= (vec1-vec2)^T Mcov^{-1} (vec1-vec2)
double chi2estimate(const TVectorD &vec1, const TVectorD &vec2,
		    const TMatrixD &Mcov)
{
  if ((vec1.GetNoElements()!=vec2.GetNoElements()) ||
      (vec1.GetNoElements()!=Mcov.GetNrows()) ||
      (vec1.GetNoElements()!=Mcov.GetNcols())) {
    std::cout << "chi2estimate: size mismatch: "
	      << Form("vec1[%d], vec2[%d], Mcov[%d,%d]",
		      vec1.GetNoElements(),vec2.GetNoElements(),
		      Mcov.GetNrows(),Mcov.GetNcols())
	      << "\n";
    return -9999.9999;
  }
  TVectorD vDiff( vec1 - vec2 );
  TMatrixD mCovInv( TMatrixD::kInverted, Mcov );
  double chi2= vDiff * ( mCovInv * vDiff );
  return chi2;
}



// ---------------------------------------------------------

int deriveCovariance(const std::vector<TH1D*> &rndCS,
		     TString histoNameTag, TString histoTitle,
		     TH1D **h1avgCS_out, TH2D **h2cov_out)
{
  int nSize= int(rndCS.size());
  if (nSize<2) return 0;

  TH1D *h1a= (TH1D*)rndCS[0]->Clone("h1_" + histoNameTag);
  h1a->SetDirectory(0);
  h1a->Sumw2();
  h1a->SetTitle(histoTitle);
  h1a->Reset();

  // accumulate the average value
  for (unsigned int i=0; i<rndCS.size(); i++) {
    h1a->Add(rndCS[i]);
  }
  h1a->Scale(1/double(nSize));

  // accumulate sqr sums
  int dim=h1a->GetNbinsX();
  TMatrixD sum2(dim+1,dim+1);
  sum2.Zero();
  for (unsigned int i=0; i<rndCS.size(); i++) {
    const TH1D *h1= rndCS[i];
    for (int ibin=1; ibin<=h1->GetNbinsX(); ibin++) {
      for (int jbin=1; jbin<=h1->GetNbinsX(); jbin++) { // ! 1D histo
	sum2(ibin,jbin) +=
	  ( h1->GetBinContent(ibin) - h1a->GetBinContent(ibin) ) *
	  ( h1->GetBinContent(jbin) - h1a->GetBinContent(jbin) );
      }
    }
  }

  sum2 *= 1/double(nSize);

  // create the covariance histogram
  const TArrayD *xb= h1a->GetXaxis()->GetXbins();
  const Double_t *x= xb->GetArray();
  TString h2name= "h2cov_" + histoNameTag;
  TH2D *h2= new TH2D(h2name,h2name, xb->GetSize()-1,x, xb->GetSize()-1,x);
  h2->SetDirectory(0);
  h2->Sumw2();
  for (int ibin=1; ibin<=h2->GetNbinsX(); ibin++) {
    for (int jbin=1; jbin<=h2->GetNbinsY(); jbin++) { // ! 2D histo
      h2->SetBinContent( ibin,jbin, sum2(ibin,jbin) );
      h2->SetBinError  ( ibin,jbin, 0.);
    }
  }

  // assign uncertainties, taking into account that the uncertainty
  // cannot be negative
  for (int ibin=1; ibin<=h2->GetNbinsX(); ibin++) {
    double unc2= h2->GetBinContent(ibin,ibin);
    if (unc2<0) {
      std::cout << "negative central covariance in ibin=" << ibin << "\n";
      unc2=0;
    }
    h1a->SetBinError(ibin, sqrt(unc2));
  }

  // assign results
  if (h1avgCS_out) *h1avgCS_out= h1a;
  if (h2cov_out) *h2cov_out= h2;
  return 1;
}

// ---------------------------------------------------------

TH2D* cov2corr(const TH2D* h2cov)
{
  TH2D* h2corr=(TH2D*)h2cov->Clone(h2cov->GetName() + TString("_corr"));
  h2corr->SetTitle(h2cov->GetTitle() + TString("_corr"));
  h2corr->SetDirectory(0);
  h2corr->SetStats(0);
  for (int ibin=1; ibin<=h2cov->GetNbinsX(); ibin++) {
    double is= sqrt(h2cov->GetBinContent(ibin,ibin));
    for (int jbin=1; jbin<=h2cov->GetNbinsY(); jbin++) {
      double js= sqrt(h2cov->GetBinContent(jbin,jbin));
      h2corr->SetBinContent(ibin,jbin, h2cov->GetBinContent(ibin,jbin)/is/js );
    }
  }
  return h2corr;
}

// ---------------------------------------------------------

TH1D *uncFromCov(const TH2D* h2cov, const TH1D *h1centralVal,
		 int zeroCentralMeansZeroRelError)
{
  const TArrayD *xb= h2cov->GetXaxis()->GetXbins();
  const Double_t *x= xb->GetArray();
  std::cout << "h2cov : " << h2cov->GetName() << " " << h2cov->GetTitle() << ", dim=" << xb->GetSize() << "\n";
  //printHisto(h2cov);
  TString h1name= TString("h1unc_from_") + h2cov->GetName();
  if (h1centralVal) h1name.ReplaceAll("h1unc","h1relUnc");
  TH1D *h1= NULL;
  if (xb->GetSize()>0) h1=new TH1D(h1name,h1name, xb->GetSize()-1,x );
  else h1= new TH1D(h1name,h1name, h2cov->GetNbinsX(), 0, h2cov->GetNbinsX()-1);
  h1->SetDirectory(0);
  h1->SetStats(0);
  //h1->Sumw2();
  for (int ibin=1; ibin<=h2cov->GetNbinsX(); ibin++) {
    double valSqr= h2cov->GetBinContent(ibin,ibin);
    double val= (valSqr<0) ? -100 : sqrt(valSqr);
    if (h1centralVal) {
      double c= h1centralVal->GetBinContent(ibin);
      if (c!=double(0)) val/=c;
      else { if (zeroCentralMeansZeroRelError) val=0.; else val=9999.; }
    }
    h1->SetBinContent(ibin, val);
    h1->SetBinError(ibin, 0);
  }
  return h1;
}

// ---------------------------------------------------------

TCanvas *plotCovCorr(TH2D* h2Cov, TString canvName, TH2D** h2corr_out,
		     int autoZRange, double ytitleOffset) {
  h2Cov->GetXaxis()->SetMoreLogLabels();
  h2Cov->GetXaxis()->SetNoExponent();
  h2Cov->GetYaxis()->SetMoreLogLabels();
  h2Cov->GetYaxis()->SetNoExponent();
  TH2D* h2Corr= cov2corr(h2Cov);
  h2Cov->SetStats(0);

  h2Cov->GetYaxis()->SetTitleOffset(ytitleOffset);
  h2Corr->GetYaxis()->SetTitleOffset(ytitleOffset);
  h2Cov->GetZaxis()->SetDecimals(true);
  h2Corr->GetZaxis()->SetDecimals(true);

  TCanvas *cx= new TCanvas(canvName,canvName,1200,600);
  cx->Divide(2,1);
  cx->cd(1);
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.15);
  h2Cov->Draw("COLZ");
  cx->cd(2);
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.15);
  if (autoZRange) h2Corr->GetZaxis()->SetRangeUser(-1,1);
  h2Corr->Draw("COLZ");
  cx->Update();
  if (h2corr_out) (*h2corr_out)=h2Corr;
  return cx;
}

// ---------------------------------------------------------

TCanvas *plotCovCorr(const TMatrixD &cov, const TH1D *h1_for_axes,
		     TString histoName, TString canvName,
		     TH2D **h2cov_out, TH2D **h2corr_out,
		     int autoZRangeCorr, double ytitleOffset)
{
  TH2D *h2cov= convert2histo(cov,h1_for_axes,histoName,histoName);
  if (h2cov_out) *h2cov_out = h2cov;
  return plotCovCorr(h2cov,canvName,h2corr_out,autoZRangeCorr,ytitleOffset);
}

// ---------------------------------------------------------
// ---------------------------------------------------------

TCanvas *findCanvas(TString canvName)
{
  TSeqCollection *seq=gROOT->GetListOfCanvases();
  TCanvas *cOut=NULL;
  for (int i=0; i<=seq->LastIndex(); i++) {
    TCanvas *c= (TCanvas*) seq->At(i);
    //std::cout << "i=" << i << " " << c->GetName() << "\n";
    if (c->GetName() == canvName) {
      cOut=c; break;
    }
  }
  return cOut;
}

// ---------------------------------------------------------

int findCanvases(TString canvNames, std::vector<TCanvas*> &cV)
{
  unsigned int iniSize=cV.size();
  TSeqCollection *seq=gROOT->GetListOfCanvases();
  for (int i=0; i<=seq->LastIndex(); i++) {
    TCanvas *c= (TCanvas*) seq->At(i);
    //std::cout << "i=" << i << " " << c->GetName() << "\n";
    int keep=0;
    if (canvNames.Length()==0) keep=1;
    else {
      std::stringstream ss(canvNames.Data());
      TString cName;
      while (!ss.eof()) {
	ss >> cName;
	if (cName.Length() && (c->GetName() == cName)) { keep=1; break; }
      }
    }
    if (keep) cV.push_back(c);
  }
  return int(cV.size()-iniSize);
}

// ---------------------------------------------------------

void closeCanvases()
{
  TSeqCollection *seq=gROOT->GetListOfCanvases();
  for (int i=0; i<=seq->LastIndex(); i++) {
    TCanvas *c= (TCanvas*) seq->At(i);
    delete c;
  }
  gSystem->ProcessEvents();
}

// ---------------------------------------------------------

void SaveCanvas(TCanvas* canv, const TString &canvName, TString destDir)
{
  gSystem->mkdir(destDir,kTRUE);
  gSystem->mkdir(destDir+TString("/png"),kTRUE);
  gSystem->mkdir(destDir+TString("/pdf"),kTRUE);
  gSystem->mkdir(destDir+TString("/root"),kTRUE);

  TString saveName=destDir+TString("/png/");
  saveName+=canvName;
  saveName+=".png";
  canv->SaveAs(saveName);
  saveName.ReplaceAll("png","pdf");
  canv->SaveAs(saveName);
  saveName.ReplaceAll("pdf","root");
  canv->SaveAs(saveName);
  return;
}

// ---------------------------------------------------------

void SaveCanvases(std::vector<TCanvas*> &cV, TString destDir, TFile *fout)
{
  if (cV.size()==0) {
    std::cout << "SaveCanvases: vector is empty\n";
    return;
  }

  gSystem->mkdir(destDir,kTRUE);
  gSystem->mkdir(destDir+TString("/png"),kTRUE);
  gSystem->mkdir(destDir+TString("/pdf"),kTRUE);
  gSystem->mkdir(destDir+TString("/root"),kTRUE);

  for (unsigned int i=0; i<cV.size(); i++) {
    TCanvas *canv= cV[i];
    TString saveName=destDir+TString("/png/");
    saveName+=canv->GetName();
    saveName+=".png";
    canv->SaveAs(saveName);
    saveName.ReplaceAll("png","pdf");
    canv->SaveAs(saveName);
    saveName.ReplaceAll("pdf","root");
    canv->SaveAs(saveName);
    if (fout) canv->Write();
  }
  return;
}

// ---------------------------------------------------------
// ---------------------------------------------------------

template<class T>
struct CompareIndicesByVectorValues {
  const std::vector<T> *_values;
  CompareIndicesByVectorValues(const std::vector<T> *v) : _values(v) {}
  bool operator()(const int &i, const int &j)
  { return _values->at(i) > _values->at(j); }
};

// ----------------------------

template<class T>
std::vector<int> getSortDescendingIdxT(const std::vector<T> &vec)
{
  std::vector<int> idx;
  for (unsigned int i=0; i<vec.size(); i++) idx.push_back(i);
  CompareIndicesByVectorValues<T> vecCmp(&vec);
  std::sort(idx.begin(),idx.end(),vecCmp);
  return idx;
}

// ---------------------------------------------------------

std::vector<int> getSortDescendingIdx(const std::vector<double> &vec)
{ return getSortDescendingIdxT(vec); }
std::vector<int> getSortDescendingIdx(const std::vector<float> &vec)
{ return getSortDescendingIdxT(vec); }

// ---------------------------------------------------------

int orderChanged(const std::vector<int> &idx) {
  int yes=0;
  for (unsigned int i=1; !yes && (i<idx.size()); i++) {
    if (idx[i-1]+1!=idx[i]) yes=1;
  }
  return yes;
}

// ---------------------------------------------------------
// ---------------------------------------------------------

std::ostream& operator<<(std::ostream &out, const TLorentzVector &v)
{
  out << Form("(Pt,Eta,Phi,En)=(%lf,%lf,%lf,%lf)",v.Pt(),v.Eta(),v.Phi(),v.E());
  return out;
}

// ---------------------------------------------------------

std::ostream& operator<<(std::ostream &out, const TLorentzVector *v)
{
  if (!v) out << "(null ptr to TLorentzVector)";
  else out << (*v);
  return out;
}

// ---------------------------------------------------------
// ---------------------------------------------------------

