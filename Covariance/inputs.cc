#include "inputs.h"

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

void plotHisto(TH1D* h1, TString cName, int logX, int logY)
{
  TCanvas *cx= new TCanvas(cName,cName,600,600);
  if (logX) cx->SetLogx();
  if (logY) cx->SetLogy();
  h1->Draw("LPE");
  cx->Update();
}

// ---------------------------------------------------------

void plotHisto(TH2D* h2, TString cName)
{
  TCanvas *cx= new TCanvas(cName,cName,600,600);
  cx->SetRightMargin(0.18);
  h2->Draw("COLZ");
  cx->Update();
}

// ---------------------------------------------------------

void plotHistoSame(TH1D *h1, TString canvName, TString drawOpt)
{
  TSeqCollection *seq=gROOT->GetListOfCanvases();
  for (int i=0; i<=seq->LastIndex(); i++) {
    TCanvas *c= (TCanvas*) seq->At(i);
    //std::cout << "i=" << i << " " << c->GetName() << "\n";
    if (c->GetName() == canvName) {
      c->cd();
      h1->Draw(drawOpt + TString(" same"));
    }
  }
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
// ---------------------------------------------------------

