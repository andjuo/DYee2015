#include "inputs.h"
#include <sstream>

// --------------------------------------------------------------

TString versionName(TVersion_t ver)
{
  TString name="versionNameUNKNOWN";
  switch(ver) {
  case _verUndef: name="UndefVer"; break;
  case _verMu1: name="Mu1"; break;
  case _verMu76X: name="Mu76X"; break;
  case _verEl1: name="El1"; break;
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

TCanvas* plotHisto(TH1D* h1, TString cName, int logX, int logY, TString drawOpt,
		   TString explain)
{
  TCanvas *cx= new TCanvas(cName,cName,600,600);
  if (logX) cx->SetLogx();
  if (logY) cx->SetLogy();
  h1->Draw(drawOpt);
  if (explain.Length()) {
    TLegend *leg= new TLegend(0.15,0.15,0.40,0.22);
    leg->SetName("myLegend");
    leg->AddEntry(h1,explain);
    leg->Draw();
  }
  cx->Update();
  return cx;
}

// ---------------------------------------------------------

TCanvas* plotHisto(TH2D* h2, TString cName, int logx, int logy)
{
  TCanvas *cx= new TCanvas(cName,cName,600,600);
  cx->SetRightMargin(0.18);
  cx->SetLogx(logx);
  cx->SetLogy(logy);
  h2->Draw("COLZ");
  cx->Update();
  return cx;
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
	  leg->SetY2NDC(leg->GetY2NDC() + 0.07);
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

void printHisto(const TH1D *h1) {
  std::cout << "\nhisto " << h1->GetName() << " " << h1->GetTitle() << "\n";
  for (int ibin=1; ibin<=h1->GetNbinsX(); ibin++) {
    std::cout << "ibin=" << ibin << " " << h1->GetBinLowEdge(ibin)
	      << " " << (h1->GetBinLowEdge(ibin)+h1->GetBinWidth(ibin))
	      << "  " << h1->GetBinContent(ibin) << " +- "
	      << h1->GetBinError(ibin) << "\n";
  }
}

// ---------------------------------------------------------

void printHisto(const TH2D *h2) {
  std::cout << "\nhisto " << h2->GetName() << " " << h2->GetTitle() << "\n";
  for (int ibin=1; ibin<=h2->GetNbinsX(); ibin++) {
    for (int jbin=1; jbin<=h2->GetNbinsY(); jbin++) {
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

void printRatio(const TH1D *h1a, const TH1D *h1b) {
  std::cout << "\nhisto A " << h1a->GetName() << " " << h1a->GetTitle()
	    << " [" << h1a->GetNbinsX() << "]\n";
  std::cout << "histo B " << h1b->GetName() << " " << h1b->GetTitle()
	    << " [" << h1b->GetNbinsX() << "]\n";
  std::cout << std::string(5,' ') << h1a->GetName() << std::string(10,' ')
	    << h1b->GetName() << std::string(10,' ') << "ratio\n";
  for (int ibin=1; ibin<=h1a->GetNbinsX(); ibin++) {
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
	      << "  " << h1a->GetBinContent(ibin)/h1b->GetBinContent(ibin)
	      << "\n";
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
	      int plotIt)
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
    h1->SetBinError  (ibin, 0.5*(yLow[ibin-1] + yHigh[ibin-1]));
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
		     int autoZRange) {
  h2Cov->GetXaxis()->SetMoreLogLabels();
  h2Cov->GetXaxis()->SetNoExponent();
  h2Cov->GetYaxis()->SetMoreLogLabels();
  h2Cov->GetYaxis()->SetNoExponent();
  TH2D* h2Corr= cov2corr(h2Cov);
  h2Cov->SetStats(0);
  TCanvas *cx= new TCanvas(canvName,canvName,1200,600);
  cx->Divide(2,1);
  cx->cd(1);
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetRightMargin(0.18);
  h2Cov->Draw("COLZ");
  cx->cd(2);
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetRightMargin(0.18);
  if (autoZRange) h2Corr->GetZaxis()->SetRangeUser(-1,1);
  h2Corr->Draw("COLZ");
  cx->Update();
  if (h2corr_out) (*h2corr_out)=h2Corr;
  return cx;
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
  std::stringstream ss(canvNames.Data());
  TString cName;
  while (!ss.eof()) {
    ss >> cName;
    TCanvas *c= findCanvas(cName);
    if (c) cV.push_back(c);
  }
  return int(cV.size()-iniSize);
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
// ---------------------------------------------------------

