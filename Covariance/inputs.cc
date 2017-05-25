#include "inputs.h"
#include <sstream>
#include <algorithm>
#include <cstring>


// --------------------------------------------------------------

const TH1D *h1dummy=NULL;
const TH2D *h2dummy=NULL;

const HistoStyle_t hsVec[6] = { hsBlack, hsBlue, hsRed, hsGreen, hsColor6,
				hsColor46 };

// --------------------------------------------------------------

TString versionName(TVersion_t ver)
{
  TString name="versionNameUNKNOWN";
  switch(ver) {
  case _verUndef: name="UndefVer"; break;
  case _verMu1: name="Mu1"; break;
  case _verMu76X: name="Mu76X"; break;
  case _verMuApproved: name="MuApproved"; break;
  case _verMuMay2017: name="MuMay2017"; break;
  case _verEl1: name="El1"; break;
  case _verEl2: name="El2"; break;
  case _verEl2skim: name="El2skim"; break;
  case _verEl2skim2: name="El2skim2"; break;
  case _verEl2skim3: name="El2skim3"; break;
  case _verEl3: name="El3"; break;
  case _verEl3mb41: name="El3mb41"; break;
  case _verEl3mb42: name="El3mb42"; break;
  case _verElMay2017: name="ElMay2017"; break;
  default:
    std::cout << "versionName is not ready for this version type\n";
  }
  return name;
}

// -----------------------------------------------------------

int leptonIdx(TVersion_t v) {
  int i=-1;
  switch(v) {
  case _verUndef: i=-1;
    break;
  case _verMu1: case _verMu76X: case _verMuApproved: case _verMuMay2017:
    i=1;
    break;
  case _verEl1: case _verEl2: case _verEl2skim: case _verEl2skim2:
  case _verEl2skim3: case _verEl3:
  case _verEl3mb41: case _verEl3mb42:
  case _verElMay2017:
    i=0;
    break;
  default:
    std::cout << "leptonIdx is not ready for version ="<< versionName(v)<< "\n";
  }
  return i;
}

// -----------------------------------------------------------
// -----------------------------------------------------------

HistoStyle_t::HistoStyle_t(int set_color, int set_markerStyle,
			   int set_lineStyle, double set_markerSize,
			   double set_lineWidth) :
  color(set_color), markerStyle(set_markerStyle), lineStyle(set_lineStyle),
  markerSize(set_markerSize), lineWidth(set_lineWidth)
{}

// -----------------------------------------------------------

HistoStyle_t::HistoStyle_t(const HistoStyle_t &hs) :
  color(hs.color), markerStyle(hs.markerStyle), lineStyle(hs.lineStyle),
  markerSize(hs.markerSize), lineWidth(hs.lineWidth)
{}

// -----------------------------------------------------------
// -----------------------------------------------------------

PlotCovCorrOpt_t:: PlotCovCorrOpt_t(int set_autoZRangeCorr, int set_gridLines,
		    int set_logScale, double set_yTitleOffset,
		    double set_leftMargin, double set_rightMargin) :
  autoZRangeCorr(set_autoZRangeCorr), gridLines(set_gridLines),
  logScaleX(set_logScale), logScaleY(set_logScale),
  setTickX(0), setTickY(1),
  yTitleOffset(set_yTitleOffset),
  leftMargin(set_leftMargin), rightMargin(set_rightMargin)
{}

// -----------------------------------------------------------

PlotCovCorrOpt_t:: PlotCovCorrOpt_t(const TString key, int value) :
  autoZRangeCorr(1), gridLines(1),
  logScaleX(1), logScaleY(1),
  setTickX(0), setTickY(1),
  yTitleOffset(1.5),
  leftMargin(0.15), rightMargin(0.15)
{
  setValue(key,value,0);
}

// -----------------------------------------------------------

PlotCovCorrOpt_t::PlotCovCorrOpt_t(const PlotCovCorrOpt_t &o) :
  autoZRangeCorr(o.autoZRangeCorr), gridLines(o.gridLines),
  logScaleX(o.logScaleX), logScaleY(o.logScaleY),
  setTickX(o.setTickX), setTickY(setTickY),
  yTitleOffset(o.yTitleOffset),
  leftMargin(o.leftMargin), rightMargin(o.rightMargin)
{}

// -----------------------------------------------------------

void PlotCovCorrOpt_t::assign(const PlotCovCorrOpt_t &opt)
{
  autoZRangeCorr= opt.autoZRangeCorr;
  gridLines= opt.gridLines;
  logScaleX= opt.logScaleX;
  logScaleY= opt.logScaleY;
  setTickX= opt.setTickX;
  setTickY= opt.setTickY;
  yTitleOffset= opt.yTitleOffset;
  leftMargin= opt.leftMargin;
  rightMargin= opt.rightMargin;
}

// -----------------------------------------------------------

int PlotCovCorrOpt_t::setValue(TString key, int value, int verbose)
{
  int res=1;
  if (key==TString("logScale")) {
    logScaleX= value;
    logScaleY= value;
    if (verbose) {
      std::cout << "PlotCovCorrOpt_t::setValue got key=" << key
		<< " changing value of logScaleX,logScaleY\n";
    }
  }
  else {
    res=0;
    if (verbose) {
      std::cout << "PlotCovCorrOpt_t::setValue got unrecognized key=" << key
		<< ". No action taken\n";
    }
  }
  return res;
}

// -----------------------------------------------------------

void PlotCovCorrOpt_t::adjustTicks(TCanvas *c) const
{
  c->SetTickx(setTickX);
  c->SetTicky(setTickY);
  c->Modified();
  c->Update();
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
    int val=va_arg(vl,int);
    //std::cout << " adding int " << val << "\n";
    vec.push_back(val);
  }
  va_end(vl);
}

// -----------------------------------------------------------

void addToVector(std::vector<double> &vec, int count, ...)
{
  va_list vl;
  va_start(vl,count);
  for (int i=0; i<count; i++) {
    double val=va_arg(vl,double);
    std::cout << " adding double " << val << "\n";
    vec.push_back(val);
  }
  va_end(vl);
}

// -----------------------------------------------------------
// -----------------------------------------------------------

TCanvas* plotHisto(TH1D* h1, TString cName, int logX, int logY, TString drawOpt,
		   TString explain, int gridLines)
{
  TCanvas *cx= new TCanvas(cName,cName,600,600);
  if (logX) cx->SetLogx();
  if (logY) cx->SetLogy();
  if (gridLines) cx->SetGrid((gridLines&1)!=0,(gridLines&2)!=0);
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
		   double ytitleOffset, int gridLines, double centralZrange)
{
  TCanvas *cx= new TCanvas(cName,cName,600,600);
  if (gridLines) cx->SetGrid(1,1);
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

TCanvas *plotHisto(TH2D* h2, TString cName, const PlotCovCorrOpt_t &o)
{
  return plotHisto(h2,cName,o.logScaleX,o.logScaleY,o.yTitleOffset,
		   o.gridLines,(o.autoZRangeCorr) ? 1. : 0.);
}

// ---------------------------------------------------------

TCanvas* plotHisto(const TH2D* h2_orig, TString cName, int logx, int logy,
		   double ytitleOffset, int gridLines, double centralZrange)
{
  TH2D *h2=cloneHisto(h2_orig,h2_orig->GetName() + TString("__clone_plotHisto"),
		      h2_orig->GetTitle());
  return plotHisto(h2,cName,logx,logy,ytitleOffset,gridLines,centralZrange);
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
	const int smaller=1;
	TLegend *leg= (TLegend*)gPad->FindObject("myLegend");
	if (!leg) {
	  leg= new TLegend(0.15,0.15,0.40,0.22-0.02*smaller);
	  leg->SetName("myLegend");
	  leg->AddEntry(h1,explain);
	  leg->Draw();
	}
	else {
	  leg->SetY2NDC(leg->GetY2NDC() + 0.07-0.02*smaller);
	  leg->AddEntry(h1,explain);
	  leg->Draw();
	}
      }

      cOut=c;
    }
  }
  cOut->Update();
  return cOut;
}

// ---------------------------------------------------------

TCanvas* plotHistoAuto(TH1D *h1, TString canvName, int logX, int logY,
		       TString drawOpt, TString explain, int gridLines) {
  TCanvas *c= findCanvas(canvName);
  if (!c) c=plotHisto(h1,canvName,logX,logY,drawOpt,explain,gridLines);
  else plotHistoSame(h1,canvName,drawOpt,explain);
  return c;
}

// ---------------------------------------------------------

TCanvas* plotHistoError(const TH1D* h1, TString cName, int logX, int logY,
			TString drawOpt, TString explain, int gridLines)
{
  TH1D *h1err= errorAsCentral(h1);
  return plotHisto(h1err,cName,logX,logY,drawOpt,explain,gridLines);
}

// ---------------------------------------------------------

TCanvas* plotHistoErrorSame(const TH1D *h1, TString canvName, TString drawOpt,
			    TString explain)
{
  TH1D *h1err= errorAsCentral(h1);
  return plotHistoSame(h1err,canvName,drawOpt,explain);
}

// ---------------------------------------------------------

TCanvas* plotGraphSame(TGraphErrors *gr, TString canvName, TString drawOpt, TString explain)
{
  TSeqCollection *seq=gROOT->GetListOfCanvases();
  TCanvas *cOut=NULL;
  for (int i=0; i<=seq->LastIndex(); i++) {
    TCanvas *c= (TCanvas*) seq->At(i);
    //std::cout << "i=" << i << " " << c->GetName() << "\n";
    if (c->GetName() == canvName) {
      c->cd();
      gr->Draw(drawOpt + TString(" same"));

      if (explain.Length()) {
	const int smaller=1;
	TLegend *leg= (TLegend*)gPad->FindObject("myLegend");
	if (!leg) {
	  leg= new TLegend(0.15,0.15,0.40,0.22-0.02*smaller);
	  leg->SetName("myLegend");
	  leg->AddEntry(gr,explain,drawOpt);
	  leg->Draw();
	}
	else {
	  leg->SetY2NDC(leg->GetY2NDC() + 0.07-0.02*smaller);
	  leg->AddEntry(gr,explain,drawOpt);
	  leg->Draw();
	}
      }

      cOut=c;
    }
  }
  cOut->Update();
  return cOut;
}

// -----------------------------------------------------------

TCanvas* plotRatio(TH1D* h1, TString cName, int logX, int logY, TString drawOpt,
		   TString explain, int gridLines, int lineAtOne)
{
  TCanvas *cx= new TCanvas(cName,cName,600,300);
  cx->SetBottomMargin(0.15);
  if (logX) cx->SetLogx();
  if (logY) cx->SetLogy();
  if (gridLines) cx->SetGrid((gridLines&1)!=0,(gridLines&2)!=0);

  ratioAxis(h1,logX+2*logY);
  h1->Draw(drawOpt);

  if (explain.Length()) {
    const int smaller=1;
    TLegend *leg= new TLegend(0.15,0.15,0.40,0.22-0.02*smaller);
    leg->SetName("myLegend");
    leg->AddEntry(h1,explain);
    leg->Draw();
  }
  if (lineAtOne) {
    int nBins= h1->GetNbinsX();
    TLine *line= new TLine(h1->GetBinLowEdge(1),1.,
			   h1->GetBinLowEdge(nBins) + h1->GetBinWidth(nBins),1.);
    line->SetLineStyle(3);
    line->SetLineColor(kBlack);
    line->Draw();
  }
  cx->Update();
  return cx;
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
  c->Modified(); c->Update();
  return 1;
}

// ---------------------------------------------------------

int increaseLegend(TCanvas *c, double dx, double dy)
{
  c->cd();
  TLegend *leg= (TLegend*)gPad->FindObject("myLegend");
  if (!leg) return 0;
  if (dx!=0.) {
    //leg->SetX1NDC(leg->GetX1NDC() + dx);
    leg->SetX2NDC(leg->GetX2NDC() + dx);
  }
  if (dy!=0.) {
    //leg->SetY1NDC(leg->GetY1NDC() + dy);
    leg->SetY2NDC(leg->GetY2NDC() + dy);
  }
  c->Modified(); c->Update();
  return 1;
}

// ---------------------------------------------------------

void setLeftMargin(TCanvas *c, double xMargin)
{
  c->cd();
  c->SetLeftMargin(xMargin);
  gPad->Modified();
  c->Modified();
  c->Update();
}

// ---------------------------------------------------------

void setRightMargin(TCanvas *c, double xMargin)
{
  c->cd();
  c->SetRightMargin(xMargin);
  gPad->Modified();
  c->Modified();
  c->Update();
}

// ---------------------------------------------------------

void setLeftRightMargins(TCanvas *c, double lMargin, double rMargin)
{
  c->cd();
  c->SetLeftMargin(lMargin);
  c->SetRightMargin(rMargin);
  gPad->Modified();
  c->Modified();
  c->Update();
}

// ---------------------------------------------------------

TCanvas *createMassFrame(int iFrame, TString canvNameBase, TString titleStr,
			 int yRangeSet,
			 TString *canvName_out, TH2D **h2frame_out)
{
  const int nDivs=5;
  const double divX[nDivs][2] = { { 15., 60.},
				  { 60., 102.},
				  { 100., 200.},
				  { 200., 600.},
				  { 500., 3000.} };
  const double divY_cs[nDivs][2] = { { 3., 300. }, // 10-60
				     { 3., 300. }, // 60-100
				     { 2e-2, 10.}, // 100-200
				     { 1e-4, 0.1}, // 200-600
				     { 1e-9,1e-3} }; // >600

  const double divY_yield[nDivs][2] = { { 0.1, 10. }, // 10-60
					{ 1., 400. }, // 60-100
					{ 0.1, 30.}, // 100-200
					{ 0.01, 1.}, // 200-600
					{ 1e-5,0.1} }; // >600

  double divY[nDivs][2];
  if (yRangeSet==_massFrameCS) memcpy(divY,divY_cs,nDivs*2*sizeof(double));
  else memcpy(divY,divY_yield,nDivs*2*sizeof(double));
  int iD=iFrame;
  TString rangeStr=Form("%d_%d",int(divX[iD][0]),int(divX[iD][1]));
  if (titleStr.Index(";")!=-1) {
    titleStr.Insert(titleStr.Index(";")," "+rangeStr);
  }
  else titleStr.Append(" " +rangeStr);
  TH2D* h2frame= new TH2D("h2frame" + rangeStr,titleStr,
			  100,divX[iD][0],divX[iD][1],
			  100,divY[iD][0],divY[iD][1]);
  h2frame->SetDirectory(0);
  h2frame->SetStats(0);
  TString cName= canvNameBase + rangeStr;
  TCanvas *cx= new TCanvas(cName,cName,600,600);
  cx->SetGrid(1,1);
  cx->SetLogy();
  logAxis(h2frame);
  if (divY[iD][0]<1e-5) h2frame->GetYaxis()->SetNoExponent(0);
  if (yRangeSet==1) {
    if (0 && (divX[iD][0]>490)) {
      cx->SetLeftMargin(0.15);
      h2frame->GetYaxis()->SetTitleOffset(2.1);
    }
    else {
      cx->SetLeftMargin(0.11);
      h2frame->GetYaxis()->SetTitleOffset(1.4);
    }
  }
  h2frame->Draw();
  cx->Update();
  if (canvName_out) *canvName_out= cName;
  if (h2frame_out) *h2frame_out= h2frame;
  return cx;
}

// ---------------------------------------------------------
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

void printHistoWLimit(const TH2D *h2, int nLinesX, int nLinesY,
		      int nonZero, double threshold)
{
  std::cout << "\nhisto " << h2->GetName() << " " << h2->GetTitle()
	    << " print maxLines " << Form("(%d,%d)",nLinesX,nLinesY) << "\n";
  for (int ibin=1; ibin<=h2->GetNbinsX(); ibin++) {
    if ((nLinesX>0) && (ibin>nLinesX)) {
      std::cout << " (discontinued)\n";
      break;
    }
    for (int jbin=1; jbin<=h2->GetNbinsY(); jbin++) {
      if ((nLinesY>0) && (jbin>nLinesY)) {
	std::cout << " ...\n";
	break;
      }
      if (nonZero) {
	double currVal=h2->GetBinContent(ibin,jbin);
	if (fabs(currVal)<threshold) continue;
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

void printHistoRange(const TH2D *h2)
{
  std::cout << "\nhisto " << h2->GetName() << " " << h2->GetTitle() << "\n";
  std::cout << " x range: ";
  for (int ibin=1; ibin<=h2->GetNbinsX(); ibin++) {
    std::cout << "  " << h2->GetXaxis()->GetBinLowEdge(ibin);
  }
  int nBinsX=h2->GetNbinsX();
  std::cout << "  " << (h2->GetXaxis()->GetBinLowEdge(nBinsX)+
			h2->GetXaxis()->GetBinWidth(nBinsX)) << "\n";
  std::cout << " y range: ";
  for (int ibin=1; ibin<=h2->GetNbinsY(); ibin++) {
    std::cout << "  " << h2->GetYaxis()->GetBinLowEdge(ibin);
  }
  int nBinsY=h2->GetNbinsY();
  std::cout << "  " << (h2->GetYaxis()->GetBinLowEdge(nBinsY)+
			h2->GetYaxis()->GetBinWidth(nBinsY)) << "\n";
}

// ---------------------------------------------------------

void printRatio(const TH1D *h1a, const TH1D *h1b, int extraRange,
		int includeErr, double markIfDiff_relTol) {
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
      std::cout << "printRatio: bining mismatch at ibin=" << ibin << "\n";
      return;
    }
    std::cout << "ibin=" << ibin << " " << h1a->GetBinLowEdge(ibin)
	      << " " << (h1a->GetBinLowEdge(ibin)+h1a->GetBinWidth(ibin))
	      << "  " << h1a->GetBinContent(ibin) << " +- "
	      << h1a->GetBinError(ibin)
	      << "  " << h1b->GetBinContent(ibin) << " +- "
	      << h1b->GetBinError(ibin);
    double r=h1a->GetBinContent(ibin)/h1b->GetBinContent(ibin);
    int isDiff=0;
    // correct for nan
    if (r!=r) {
      if ((h1a->GetBinContent(ibin)==0.) && (h1b->GetBinContent(ibin)==0.)) {
	r=1.;
      }
      else isDiff=1;
    }
    std::cout << "  " << r;
    if (includeErr)
      std::cout << " " << h1a->GetBinError(ibin)/h1b->GetBinError(ibin);

    if ((markIfDiff_relTol!=0.) && !isDiff) {
      double d= h1a->GetBinContent(ibin);
      if (d==0.) d=h1b->GetBinContent(ibin);
      if (fabs((h1a->GetBinContent(ibin)-h1b->GetBinContent(ibin))/d) >
	  markIfDiff_relTol) isDiff=1;
      if (!isDiff && includeErr) {
	double d= h1a->GetBinError(ibin);
	if (d==0.) d=h1b->GetBinError(ibin);
	if ((d!=0.) &&
	    (fabs((h1a->GetBinError(ibin)-h1b->GetBinError(ibin))/d) >
	     markIfDiff_relTol)) isDiff=2;
      }
    }
    if (isDiff) std::cout << std::string(isDiff,'*');
    std::cout << "\n";
  }
}


// ---------------------------------------------------------

void printRatio(const TH2D *h2a, const TH2D *h2b, int extraRange,
		int includeErr, double markIfDiff_relTol) {
  int d=(extraRange) ? 1:0;
  std::cout << "\nhisto A " << h2a->GetName() << " " << h2a->GetTitle()
	    << " [" << h2a->GetNbinsX() << "]\n";
  std::cout << "histo B " << h2b->GetName() << " " << h2b->GetTitle()
	    << " [" << h2b->GetNbinsX() << "]\n";
  std::cout << std::string(5,' ') << h2a->GetName() << std::string(10,' ')
	    << h2b->GetName() << std::string(10,' ') << "ratio\n";
  for (int ibin=1-d; ibin<=h2a->GetNbinsX()+d; ibin++) {
    if (ibin-d>h2b->GetNbinsX()) {
      std::cout << "printRatio: different number of x bins\n";
      continue;
    }
    if ( (h2a->GetXaxis()->GetBinLowEdge(ibin) !=
	  h2b->GetXaxis()->GetBinLowEdge(ibin)) ||
	 (h2a->GetXaxis()->GetBinWidth(ibin) !=
	  h2b->GetXaxis()->GetBinWidth(ibin)) ) {
      std::cout << "printRatio: x-bining mismatch at ibin=" << ibin << "\n";
      //return;
    }
    for (int jbin=1-d; jbin<=h2a->GetNbinsY()+d; jbin++) {
      if (jbin-d>h2b->GetNbinsY()) {
	std::cout << "different number of y bins\n";
	continue;
      }
      if ( (h2a->GetYaxis()->GetBinLowEdge(jbin) !=
	    h2b->GetYaxis()->GetBinLowEdge(jbin)) ||
	   (h2a->GetYaxis()->GetBinWidth(jbin) !=
	    h2b->GetYaxis()->GetBinWidth(jbin)) ) {
	std::cout << "printRatio: y-bining mismatch at jbin=" << jbin << "\n";
	//return;
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
		<< h2b->GetBinError(ibin,jbin);
      double r=h2a->GetBinContent(ibin,jbin)/h2b->GetBinContent(ibin,jbin);
      int isDiff=0;
      if (r!=r) {
	if ((h2a->GetBinContent(ibin,jbin)==0.) &&
	    (h2b->GetBinContent(ibin,jbin)==0.)) {
	  r=1;
	}
	else isDiff=1;
      }
      std::cout << "  " << r;
      if (includeErr) {
	std::cout << " " << h2a->GetBinError(ibin,jbin)/h2b->GetBinError(ibin,jbin);
      }

      if ((markIfDiff_relTol!=0.) && !isDiff) {
	double d= h2a->GetBinContent(ibin,jbin);
	if (d==0.) d=h2b->GetBinContent(ibin,jbin);
	if (fabs((h2a->GetBinContent(ibin,jbin)-h2b->GetBinContent(ibin,jbin))/d)
	    > markIfDiff_relTol) {
	  //std::cout << "\ndebug: h2a=" << h2a->GetBinContent(ibin,jbin) << ", h2b=" << h2b->GetBinContent(ibin,jbin) << ", d=" << d << ", diff=" << (h2a->GetBinContent(ibin,jbin)/h2b->GetBinContent(ibin,jbin)) << ", the ratio=" << fabs((h2a->GetBinContent(ibin,jbin)-h2b->GetBinContent(ibin,jbin))/d) << ", markIfDiff_relTol=" << markIfDiff_relTol << "\n";
	  isDiff=1;
	}
	if (!isDiff && includeErr) {
	  double d= h2a->GetBinError(ibin,jbin);
	  if (d==0.) d= h2b->GetBinError(ibin,jbin);
	  if ((d!=0.) &&
	      (fabs((h2a->GetBinError(ibin,jbin)-h2b->GetBinError(ibin,jbin))/d)
	       > markIfDiff_relTol)) isDiff=2;
	}
      }
      if (isDiff) std::cout << " " << std::string(isDiff,'*');
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

// if h1ref, differences are taken wrt to 1st histo in the array
int deriveRelSyst(const TH1D *h1ref_inp, const std::vector<TH1D*> &h1V,
		  std::vector<TH1D*> &h1diffV, int relative, int absValues)
{
  if ((h1ref_inp==NULL) && (h1V.size()<2)) {
    std::cout << "deriveRelSyst: too few input histos\n";
    return 0;
  }
  unsigned int startIdx= (h1ref_inp==NULL) ? 1 : 0;
  const TH1D *h1ref= (h1ref_inp==NULL) ? h1V[0] : h1ref_inp;
  h1diffV.reserve(h1V.size());
  for (unsigned int i=startIdx; i<h1V.size(); i++) {
    TString op= (relative) ? "_relDiffDiv_" : "_diff_";
    TString hname=h1V[i]->GetName() + op + h1ref->GetName();
    TH1D *h1diff= cloneHisto(h1V[i],hname,hname);
    h1diff->Add(h1ref,-1);
    if (relative) h1diff->Divide(h1ref);
    removeError(h1diff);
    if (absValues) {
      for (int ibin=1; ibin<=h1diff->GetNbinsX(); ibin++) {
	h1diff->SetBinContent(ibin, fabs(h1diff->GetBinContent(ibin)));
      }
    }
    h1diffV.push_back(h1diff);
  }
  return 1;
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

TH2D* errorAsCentral(const TH2D* h2, int relative) {
  TH2D* h2err= cloneHisto(h2,
			  h2->GetName() + TString("_err"),
			  h2->GetTitle() + TString(" err"));
  for (int ibin=1; ibin<=h2->GetNbinsX(); ibin++) {
    for (int jbin=1; jbin<=h2->GetNbinsY(); jbin++) {
      double err= h2->GetBinError(ibin,jbin);
      if (relative) {
	double val=h2->GetBinContent(ibin,jbin);
	if (val==double(0)) {
	  if (err==double(0)) err=0;
	  else if (err<double(0)) err=-1e4;
	  else err=1e4;
	}
	else err/=val;
      }
      h2err->SetBinContent(ibin,jbin, err);
      h2err->SetBinError  (ibin,jbin, 0.);
    }
  }
  return h2err;
}

// ---------------------------------------------------------

TH1D *addUncertainty(const TString newName, const TH1D *h1_dest,
		     const TH1D *h1err, const TH1D *h1central)
{
  TH1D *h1= cloneHisto(h1_dest,newName,newName);
  for (int ibin=1; ibin<=h1->GetNbinsX(); ibin++) {
    double err1= h1_dest->GetBinError(ibin);
    double err2= h1err->GetBinContent(ibin);
    if (h1central) err2*= h1central->GetBinContent(ibin);
    h1->SetBinError(ibin, sqrt(err1*err1 + err2*err2));
  }
  return h1;
}

// ---------------------------------------------------------

int setError(TH1D *h1dest, const TH1D* h1src)
{
  if (!compareRanges(h1dest,h1src,0)) {
    std::cout << "setError(TH1D): different ranges\n";
    return 0;
  }
  for (int ibin=1; ibin<=h1src->GetNbinsX(); ibin++) {
    h1dest->SetBinError(ibin, h1src->GetBinError(ibin));
  }
  return 1;
}

// ---------------------------------------------------------

int setError(TH2D *h2dest, const TH2D* h2src)
{
  if (!compareRanges(h2dest,h2src,0)) {
    std::cout << "setError(TH2D): different ranges\n";
    return 0;
  }
  for (int ibin=1; ibin<=h2src->GetNbinsX(); ibin++) {
    for (int jbin=1; jbin<=h2src->GetNbinsY(); jbin++) {
      h2dest->SetBinError(ibin,jbin, h2src->GetBinError(ibin,jbin));
    }
  }
  return 1;
}

// ---------------------------------------------------------

int allGE(const TH1D *h1a, const TH1D *h1b)
{
  if (h1a->GetNbinsX() != h1b->GetNbinsX()) {
    std::cout << "allGE: different number of bins: "
	      << h1a->GetNbinsX() << " vs " << h1b->GetNbinsX() << "\n";
    return 0;
  }
  int ok=1;
  for (int ibin=1; ok && (ibin<=h1a->GetNbinsX()); ibin++) {
    if (h1a->GetBinContent(ibin) < h1b->GetBinContent(ibin)) ok=0;
  }
  return ok;
}

// ---------------------------------------------------------

int assignDiffAsUnc(TH2D *h2target, const TH2D *h2nominal, int relative,
		    int listRangesOnError)
{
  if (!compareRanges(h2target,h2nominal,listRangesOnError)) {
    std::cout << "assignDiffAsUnc: ranges are different\n";
    std::cout << "nominal:  "; printHisto(h2nominal);
    std::cout << "target :  "; printHisto(h2target);
    return 0;
  }
  for (int ibin=1; ibin<=h2target->GetNbinsX(); ibin++) {
    for (int jbin=1; jbin<=h2target->GetNbinsY(); jbin++) {
      double diff= h2target->GetBinContent(ibin,jbin) -
	h2nominal->GetBinContent(ibin,jbin);
      double vCentral= 0;
      if (relative) {
	vCentral=1;
	if (h2nominal->GetBinContent(ibin,jbin)!=0) {
	  diff/= h2nominal->GetBinContent(ibin,jbin);
	}
	else if (diff!=0) {
	  diff/= h2target->GetBinContent(ibin,jbin);
	}
      }
      h2target->SetBinContent(ibin,jbin, vCentral);
      h2target->SetBinError  (ibin,jbin, diff);
    }
  }
  return 1;
}

// ---------------------------------------------------------

int addShiftByUnc(TH2D *h2target, const TH2D *h2src, double nSigmas,
		  int listRangesOnError)
{
  if (!compareRanges(h2target,h2src,listRangesOnError)) {
    std::cout << "addShiftByUnc: ranges are different\n";
    std::cout << "src:  "; printHisto(h2src);
    std::cout << "target :  "; printHisto(h2target);
    return 0;
  }
  for (int ibin=1; ibin<=h2target->GetNbinsX(); ibin++) {
    for (int jbin=1; jbin<=h2target->GetNbinsY(); jbin++) {
      double v= h2target->GetBinContent(ibin,jbin);
      h2target->SetBinContent(ibin,jbin,
			      v + nSigmas * h2src->GetBinError(ibin,jbin) );
    }
  }
  return 1;
}

// ---------------------------------------------------------

int addInQuadrature(TH1D *h1target, const TH1D *h1addTerm, int nullifyErr)
{
  if (!compareRanges(h1target,h1addTerm,1)) {
    std::cout << "error in addInQuadrature\n";
    return 0;
  }
  for (int ibin=1; ibin<=h1target->GetNbinsX(); ibin++) {
    double a= h1target->GetBinContent(ibin);
    const double b= h1addTerm->GetBinContent(ibin);
    h1target->SetBinContent(ibin, sqrt(a*a + b*b));
    if (nullifyErr) h1target->SetBinError  (ibin, 0.);
  }
  return 1;
}


// ---------------------------------------------------------

int hasValueAbove(const TH2D* h2, double limit)
{
  int ok=1;
  for (int ibin=1; ok && (ibin<=h2->GetNbinsX()); ibin++) {
    for (int jbin=1; ok && (jbin<=h2->GetNbinsY()); jbin++) {
      if (h2->GetBinContent(ibin,jbin)>limit) ok=0;
    }
  }
  return (1-ok);
}

// ---------------------------------------------------------

int hasValueBelow(const TH2D* h2, double limit)
{
  int ok=1;
  for (int ibin=1; ok && (ibin<=h2->GetNbinsX()); ibin++) {
    for (int jbin=1; ok && (jbin<=h2->GetNbinsY()); jbin++) {
      if (h2->GetBinContent(ibin,jbin)<limit) ok=0;
    }
  }
  return (1-ok);
}

// ---------------------------------------------------------

void removeNegatives(TH1D* h1)
{
  for (int ibin=1; ibin<=h1->GetNbinsX(); ibin++) {
    if (h1->GetBinContent(ibin)<0) h1->SetBinContent(ibin,0.);
  }
}

// ---------------------------------------------------------

int compareRanges(const TH1D *h1a, const TH1D *h1b, int verbose)
{
  int d=0;
  if (h1a->GetNbinsX()!=h1b->GetNbinsX()) {
    if (verbose) std::cout << "different number of x bins: "
			   << h1a->GetNbinsX() << " vs "
			   << h1b->GetNbinsX() << "\n";
    return 0;
  }
  for (int ibin=1-d; ibin<=h1a->GetNbinsX()+d; ibin++) {
    if (h1a->GetXaxis()->GetBinLowEdge(ibin) !=
	h1b->GetXaxis()->GetBinLowEdge(ibin)) {
      if (verbose) std::cout << "compareRanges: x-bining mismatch at ibin="
			     << ibin << ": "
			     << h1a->GetXaxis()->GetBinLowEdge(ibin)
			     << " vs "
			     << h1b->GetXaxis()->GetBinLowEdge(ibin)
			     << "\n";
      return 0;
    }
  }
  return 1;
}

// ---------------------------------------------------------

int compareRanges(const TH2D *h2a, const TH2D *h2b, int verbose)
{
  int d=0;
  for (int ibin=1-d; ibin<=h2a->GetNbinsX()+d; ibin++) {
    if ( (h2a->GetXaxis()->GetBinLowEdge(ibin) !=
	  h2b->GetXaxis()->GetBinLowEdge(ibin)) ||
	 (h2a->GetXaxis()->GetBinWidth(ibin) !=
	  h2b->GetXaxis()->GetBinWidth(ibin)) ) {
      if (verbose) {
	std::cout << "compareRanges: x-bining mismatch at ibin=" << ibin << ": "
		  << h2a->GetXaxis()->GetBinLowEdge(ibin) << " w:"
		  << h2a->GetXaxis()->GetBinWidth(ibin) << " vs "
		  << h2b->GetXaxis()->GetBinLowEdge(ibin) << " w:"
		  << h2b->GetXaxis()->GetBinWidth(ibin) << "\n";
      }
      return 0;
    }
    for (int jbin=1-d; jbin<=h2a->GetNbinsY()+d; jbin++) {
      if ( (h2a->GetYaxis()->GetBinLowEdge(jbin) !=
	    h2b->GetYaxis()->GetBinLowEdge(jbin)) ||
	   (h2a->GetYaxis()->GetBinWidth(jbin) !=
	    h2b->GetYaxis()->GetBinWidth(jbin)) ) {
	if (verbose) {
	  std::cout << "compareRanges: y-bining mismatch at jbin="<<jbin << ": "
		    << h2a->GetYaxis()->GetBinLowEdge(jbin) << " w:"
		    << h2a->GetYaxis()->GetBinWidth(jbin) << " vs "
		    << h2b->GetYaxis()->GetBinLowEdge(jbin) << " w:"
		    << h2b->GetYaxis()->GetBinWidth(jbin) << "\n";
	}
	return 0;
      }
    }
  }
  return 1;
}

// ---------------------------------------------------------

int checkRange(const TH1D *h1, double yrangeMin, double yrangeMax, int silent)
{
  int ok=1;
  for (int ibin=1; ibin<=h1->GetNbinsX(); ibin++) {
    double v=h1->GetBinContent(ibin);
    double dv=h1->GetBinError(ibin);
    if ((v-dv<yrangeMin) || (v+dv>yrangeMax)) {
      ok=0;
      if (!silent) {
	std::cout << "Warning histo " << h1->GetName() << " has value v=" << v
		  << " that is not in range "
		  << yrangeMin << " .. " << yrangeMax << "\n";
      }
    }
  }
  return ok;
}

// ---------------------------------------------------------

int checkRange(const std::vector<TH1D*> &h1V,
	       double &rangeMin, double &rangeMax,
	       const std::vector<std::pair<double,double> > &ranges,
	       int ignoreZeroValues, int silent)
{
  // determine the widest y range
  double vMin=1e9, vMax=-1e9;
  for (unsigned int i=0; i<h1V.size(); i++) {
    const TH1D *h1= h1V[i];
    for (int ibin=1; ibin<=h1->GetNbinsX(); ibin++) {
      double v =h1->GetBinContent(ibin);
      double dv=h1->GetBinError(ibin);
      if (ignoreZeroValues && (v==0)) continue;
      if (v-dv < vMin) vMin= v-dv;
      if (v+dv > vMax) vMax= v+dv;
    }
  }

  int ok=0;
  for (std::vector<std::pair<double,double> >::const_iterator it=ranges.begin();
       !ok && (it!=ranges.end()); it++) {
    std::pair<double,double> p=*it;
    if ((p.first<vMin) && (p.second>vMax)) {
      ok=1;
      rangeMin=p.first;
      rangeMax=p.second;
    }
  }

  if (!ok && !silent) {
    std::cout << "checkRange(h1V): failed to find suitable range for "
	      << vMin << " " << vMax << "\n";
  }

  return ok;
}

// ---------------------------------------------------------

void scaleBin(TH1D* h1, int ibin, double x)
{
  h1->SetBinContent(ibin, x * h1->GetBinContent(ibin));
  h1->SetBinError  (ibin, x * h1->GetBinError(ibin));
}

// ---------------------------------------------------------

void printBin(const TH1D *h1, int ibin, int newLine)
{
  std::cout << h1->GetName() << "[" << ibin << "]="
	    << h1->GetBinContent(ibin) << " +- "
	    << h1->GetBinError(ibin);
  if (newLine) std::cout << "\n";
}

// ---------------------------------------------------------

int hasDoubleZero(const TH2D *h2)
{
  int ok=1;
  for (int ibin=1; ok && (ibin<=h2->GetNbinsX()); ibin++) {
    for (int jbin=1; ok && (jbin<=h2->GetNbinsY()); jbin++) {
      if ((h2->GetBinContent(ibin,jbin) == 0) &&
	  (h2->GetBinError(ibin,jbin) == 0)) ok=0;
    }
  }
  return (!ok);
}

// ---------------------------------------------------------

void setToOne(TH2D *h2)
{
  for (int ibin=1; ibin<=h2->GetNbinsX(); ibin++) {
    for (int jbin=1; jbin<=h2->GetNbinsY(); jbin++) {
      h2->SetBinContent(ibin,jbin, 1.);
      h2->SetBinError  (ibin,jbin, 0.);
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

TH2D* perMassBinWidth(const TH2D* h2_orig, int prnBinW)
{
  TH2D* h2=(TH2D*)h2_orig->Clone(h2_orig->GetName() + TString("_perMassBinW"));
  for (int ibin=1; ibin<=h2->GetNbinsX(); ibin++) {
    double dX= h2->GetXaxis()->GetBinWidth(ibin);
    for (int jbin=1; jbin<=h2->GetNbinsY(); jbin++) {
      double dY= h2->GetYaxis()->GetBinWidth(jbin);
      if (prnBinW) {
	std::cout << "ibin=" << ibin << ", jbin=" << jbin
		  << ", dX=" << dX << ", dY=" << dY << "\n";
      }
      h2->SetBinContent( ibin,jbin, h2_orig->GetBinContent(ibin,jbin)/(dX*dY) );
      h2->SetBinError  ( ibin,jbin, h2_orig->GetBinError(ibin,jbin)/(dX*dY)   );
    }
  }
  return h2;
}

// ---------------------------------------------------------

TH1D* timesMassBinWidth(const TH1D* h1_orig, int prnBinW)
{
  TH1D* h1=(TH1D*)h1_orig->Clone(h1_orig->GetName() + TString("_timesMassBinW"));
  for (int ibin=1; ibin<=h1->GetNbinsX(); ibin++) {
    double w= h1->GetBinWidth(ibin);
    if (prnBinW) std::cout << "ibin=" << ibin << ", w=" << w << "\n";
    h1->SetBinContent( ibin, h1_orig->GetBinContent(ibin)*w );
    h1->SetBinError  ( ibin, h1_orig->GetBinError(ibin)*w   );
  }
  return h1;
}

// ---------------------------------------------------------

TH2D* timesMassBinWidth(const TH2D* h2_orig, int prnBinW)
{
  TH2D* h2=(TH2D*)h2_orig->Clone(h2_orig->GetName() +TString("_timesMassBinW"));
  for (int ibin=1; ibin<=h2->GetNbinsX(); ibin++) {
    double dx= h2->GetXaxis()->GetBinWidth(ibin);
    for (int jbin=1; jbin<=h2->GetNbinsY(); jbin++) {
      double dy= h2->GetYaxis()->GetBinWidth(jbin);
      if (prnBinW) std::cout << "ibin=" << ibin << ", jbin=" << jbin
			     << ", dx*dy=" << dx << "*" << dy << "\n";
      h2->SetBinContent( ibin,jbin, h2_orig->GetBinContent(ibin,jbin)*dx*dy );
      h2->SetBinError  ( ibin,jbin, h2_orig->GetBinError(ibin,jbin)*dx*dy   );
    }
  }
  return h2;
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

TCanvas *loadCanvas(TFile &fin, TString nameOnFile, TString newName,
		    int warnIfMissing)
{
  TCanvas *c= (TCanvas*)fin.Get(nameOnFile);
  if (!c) {
    if (warnIfMissing) std::cout << "failed to find canvas " << nameOnFile
				 << " in file <" << fin.GetName() << ">\n";
    return NULL;
  }
  c->SetName(newName);
  return c;
}

// ---------------------------------------------------------

TCanvas *loadCanvas(TString fileName, TString nameOnFile, TString newName,
		    int warnIfMissing)
{
  TFile fin(fileName);
  if (!fin.IsOpen()) {
    std::cout << "loadCanvas: failed to open the file <" << fileName << ">\n";
    return NULL;
  }
  TCanvas *c= loadCanvas(fin,nameOnFile,newName,warnIfMissing);
  fin.Close();
  return c;
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

void randomizeWithinErr(const TH1D *hSrc, TH1D *hDest, int nonNegative,
			int poissonRnd, int systematic)
{
  if (!hSrc) { std::cout << "randomizeWithinErr(1D): hSrc is null\n"; return; }
  if (!hDest) { std::cout << "randomizeWithinErr(1D): hDest is null\n"; return; }
  if (!systematic) {
    for (int ibin=1; ibin<=hSrc->GetNbinsX(); ibin++) {
      double val=	(poissonRnd) ?
	gRandom->Poisson(hSrc->GetBinContent(ibin)) :
	gRandom->Gaus(hSrc->GetBinContent(ibin),
		      hSrc->GetBinError(ibin));
      if (nonNegative && (val<0)) val=0;
      hDest->SetBinContent(ibin,val);
      hDest->SetBinError  (ibin,0); //hSrc->GetBinError(ibin));
    }
  }
  else {
    if (poissonRnd) {
      hDest->Reset();
      std::cout << "randomizeWithinErr(TH1D*): poissonRnd is not ready for "
		<< "systematic uncertainty\n";
      return;
    }
    double rnd=	gRandom->Gaus(0.,1.);
    for (int ibin=1; ibin<=hSrc->GetNbinsX(); ibin++) {
      double val= hSrc->GetBinContent(ibin) + rnd * hSrc->GetBinError(ibin);
      if (nonNegative && (val<0)) val=0;
      hDest->SetBinContent(ibin,val);
      hDest->SetBinError  (ibin,0);
    }
  }
}

// ---------------------------------------------------------

void randomizeWithinRelErr(const TH1D *hSrc, TH1D *hDest, double relErr,
			   int nonNegative, int systematic)
{
  if (!hSrc) { std::cout << "randomizeWithinRelErr(1D): hSrc is null\n"; return; }
  if (!hDest) { std::cout << "randomizeWithinRelErr(1D): hDest is null\n"; return; }
  double rnd0= gRandom->Gaus(0.,1.);
  for (int ibin=1; ibin<=hSrc->GetNbinsX(); ibin++) {
    double rnd= (systematic) ? rnd0 : gRandom->Gaus(0.,1.);
    double val=	hSrc->GetBinContent(ibin) +
      rnd * relErr * hSrc->GetBinContent(ibin);
    if (nonNegative && (val<0)) val=0;
    hDest->SetBinContent(ibin,val);
    hDest->SetBinError  (ibin,0); //hSrc->GetBinError(ibin));
  }
}

// ---------------------------------------------------------

void randomizeWithinErr(const TH2D *hSrc, TH2D *hDest, int nonNegative,
			int poissonRnd, int systematic)
{
  if (!hSrc) { std::cout << "randomizeWithinErr(2D): hSrc is null\n"; return; }
  if (!hDest) { std::cout << "randomizeWithinErr(2D): hDest is null\n"; return; }
  if (!systematic) {
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
  else {
    if (poissonRnd) {
      hDest->Reset();
      std::cout << "randomizeWithinErr(TH2D): Poisson error is not ready with "
		<< "flag systematic\n";
      return;
    }
    double rnd= gRandom->Gaus(0.,1.);
    for (int ibin=1; ibin<=hSrc->GetNbinsX(); ibin++) {
      for (int jbin=1; jbin<=hSrc->GetNbinsY(); jbin++) {
	double val= hSrc->GetBinContent(ibin,jbin) +
	  rnd * hSrc->GetBinError(ibin,jbin);
	if (nonNegative && (val<0)) val=0;
	hDest->SetBinContent(ibin,jbin,val);
	hDest->SetBinError  (ibin,jbin,0); //hSrc->GetBinError(ibin,jbin));
      }
    }
  }
}

// -----------------------------------------------------------

void randomizeWithinRelErrSetCentral(const TH2D *h2Central,
				     const TH2D *h2RelUnc,
				     TH2D *hDest, int nonNegative,
				     int systematic)
{
  if (!h2Central) { std::cout << "randomizeWithinRelErrSetCentral(2D): h2Central is null\n"; return; }
  if (!h2RelUnc) { std::cout << "randomizeWithinRelErrSetCentral(2D): h2RelUnc is null\n"; return; }
  if (!hDest) { std::cout << "randomizeWithinRelErrSetCentral(2D): hDest is null\n"; return; }
  double rnd0= gRandom->Gaus(0,1);
  for (int ibin=1; ibin<=h2RelUnc->GetNbinsX(); ibin++) {
    for (int jbin=1; jbin<=h2RelUnc->GetNbinsY(); jbin++) {
      double rnd= (systematic) ? rnd0 : gRandom->Gaus(0,1);
      double v0 = h2Central->GetBinContent(ibin,jbin);
      double val= v0 + rnd * v0 * h2RelUnc->GetBinError(ibin,jbin);
      if (nonNegative && (val<0)) val=0;
      hDest->SetBinContent(ibin,jbin,val);
      hDest->SetBinError  (ibin,jbin,0); //hSrc->GetBinError(ibin,jbin));
    }
  }
}

// -----------------------------------------------------------

int randomizeWithinPosNegErr(const TH2D *h2_errPos,
			     const TH2D *h2_errNeg,
			     TH2D *h2Dest,
			     int bound01, int systematic)
{
  if (!h2_errPos || !h2_errNeg || !h2Dest) {
    std::cout << "randomizeWithinPosNegErr: Null ptrs detected\n";
    return 0;
  }
  double rnd0= gRandom->Gaus(0.,1.);
  for (int ibin=1; ibin<=h2Dest->GetNbinsX(); ibin++) {
    for (int jbin=1; jbin<=h2Dest->GetNbinsY(); jbin++) {
      double r= (systematic) ? rnd0 : gRandom->Gaus(0.,1.);
      double err=0;
      if (r<0) err= h2_errNeg->GetBinError(ibin,jbin);
      else if (r>0) err= h2_errPos->GetBinError(ibin,jbin);
      double v= h2_errPos->GetBinContent(ibin,jbin) + r * err;
      if (bound01) {
	if (v<0.) v=0.;
	else if (v>1.) v=1.;
      }
      h2Dest->SetBinContent(ibin,jbin, v);
      h2Dest->SetBinError  (ibin,jbin, 0.);
    }
  }
  return 1;
}


// ---------------------------------------------------------
// ---------------------------------------------------------

void printMDim(const TString name, const TMatrixD &m)
{ std::cout << name << "[" << m.GetNrows() << "x" << m.GetNcols() << "]"; }

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
  TH1D *h1axes=NULL;
  if (h1_for_axes) h1axes= cloneHisto(h1_for_axes,"h1axes","h1axes");
  else h1axes= new TH1D("h1axes","h1axes;binIdx;count",
			m.GetNrows(),1,m.GetNrows()+1);
  if (m.GetNrows() != h1axes->GetNbinsX()) {
    std::cout << "convert2histo: matrix has different number of elements "
	      << "than histogram provides bins\n";
    std::cout << "m[" << m.GetNrows() << "," << m.GetNcols() << "], histo["
	      << h1axes->GetNbinsX() << "\n";
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

  const TArrayD *xb= h1axes->GetXaxis()->GetXbins();
  int nBins=xb->GetSize();
  Double_t *myXArr=NULL;
  if (nBins==0) {
    //std::cout << "convert2histo: histogram has fixed axis. This case"
    //	      << " is not implemented\n";
    //return NULL;
    nBins= h1axes->GetNbinsX();
    myXArr= new Double_t[nBins+1];
    for (int ibin=1; ibin<=nBins; ibin++) {
      myXArr[ibin-1]= h1axes->GetBinLowEdge(ibin);
    }
    myXArr[nBins]= h1axes->GetBinLowEdge(nBins) +
      h1axes->GetBinWidth(nBins);
  }
  const Double_t *x= (xb->GetSize()!=0) ? xb->GetArray() : myXArr;
  TH2D *h2= new TH2D(h2name,h2title, nBins-1,x, nBins-1,x);
  h2->SetDirectory(0);
  h2->SetStats(0);
  h2->Sumw2();
  if (h2title.Index(";")==-1) {
    TString label=h1axes->GetXaxis()->GetTitle();
    h2->GetXaxis()->SetTitle(label);
    h2->GetYaxis()->SetTitle(label);
  }

  for (int ibin=1; ibin<nBins; ibin++) {
    for (int jbin=1; jbin<nBins; jbin++) {
      h2->SetBinContent( ibin, jbin, m(ibin-1,jbin-1) );
      double err= (mErr) ? (*mErr)(ibin-1,jbin-1) : 0.;
      h2->SetBinError( ibin, jbin, err );
    }
  }
  if (myXArr) delete myXArr;
  if (h1axes) delete h1axes;
  return h2;
}

// ---------------------------------------------------------

TH1D* convert2histo1D(const TMatrixD &m, const TMatrixD &cov,
		      TString h1name, TString h1title,
		      const TH1D *h1_for_axes, int iCol)
{
  if (m.GetNcols()<=iCol) {
    std::cout << "convert2histo1D: matrix has fewer columns than requested iCol\n";
    std::cout << "  dims: " << m.GetNrows() << " x " << m.GetNcols()
	      << ", iCol=" << iCol << "\n";
    return NULL;
  }
  if (h1_for_axes && (m.GetNrows() != h1_for_axes->GetNbinsX())) {
    std::cout << "convert2histo1D: matrix has different number of elements "
	      << "than histogram provides bins\n";
    std::cout << "m[" << m.GetNrows() << "," << m.GetNcols() << "], histo["
	      << h1_for_axes->GetNbinsX() << "\n";
    return NULL;
  }
  if ((m.GetNrows() != cov.GetNrows()) || (m.GetNrows() != cov.GetNcols())) {
    std::cout << "convert2histo1D: meas "
	      << m.GetNrows() << "x" << m.GetNcols() << ", cov "
	      << cov.GetNrows() << "x" << cov.GetNcols() << "\n";
    return NULL;
  }

  TH1D *h1= NULL;
  if (h1_for_axes) h1=(TH1D*)h1_for_axes->Clone(h1name);
  else h1= new TH1D(h1name,h1name,m.GetNrows(),1.,m.GetNrows()+1.);
  h1->SetStats(0);
  h1->SetDirectory(0);
  h1->SetTitle(h1title);

  for (int ibin=1; ibin<=h1->GetNbinsX(); ibin++) {
    int jbin=(iCol==-1) ? ibin : (iCol+1);
    h1->SetBinContent( ibin, m(ibin-1,jbin-1) );
    h1->SetBinError( ibin, sqrt(cov(ibin-1,ibin-1)) );
  }
  return h1;
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

int changeCov(TMatrixD &cov, const TH1D *h1unc)
{
  TMatrixD newCov(TMatrixD::kZero,cov);
  if (cov.GetNrows()!=h1unc->GetNbinsX()) {
    std::cout << "changeCov detected dim difference: "
	      << cov.GetNrows() << " vs " << h1unc->GetNbinsX() << "\n";
    cov=newCov;
    return 0;
  }

  for (int ir=0; ir<cov.GetNrows(); ir++) {
    for (int ic=0; ic<cov.GetNcols(); ic++) {
      double corr= cov(ir,ic) / sqrt(cov(ir,ir)*cov(ic,ic));
      newCov(ir,ic) = corr *
	h1unc->GetBinContent(ir+1) * h1unc->GetBinContent(ic+1);
    }
  }

  cov=newCov;
  return 1;
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

double sumCov(const TMatrixD &cov, int idxMin, int idxMax,
	      const TH1D *h1binWidth_def, int printNumbers)
{
  if ((idxMax<idxMin) || (idxMin<0)) {
    std::cout << "sumCov: incorrect idx values: "
	      << idxMin << " " << idxMax << "\n";
    return 0;
  }
  double sum=0;
  for (int ir=idxMin; ir<idxMax; ir++) {
    for (int ic=idxMin; ic<idxMax; ic++) {
      double val= cov(ir,ic);
      if (printNumbers) {
	std::cout << Form("(ir,ic)=(%d,%d) ",ir,ic) << " val=" << val;
      }
      if (h1binWidth_def) {
	val*=
	  h1binWidth_def->GetBinWidth(ir+1) * h1binWidth_def->GetBinWidth(ic+1);
	if (printNumbers) {
	  std::cout << " * [bw] " << h1binWidth_def->GetBinWidth(ir+1) << " * "
		    << h1binWidth_def->GetBinWidth(ic+1);
	}
      }
      if (printNumbers) std::cout << "\n";
      sum+= val;
    }
  }
  return sum;
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

TH1D *uncFromCov(const TMatrixD &covM,
		 TString hName,
		 const TH1D *h1_binning,
		 const TH1D *h1centralVal,
		 int zeroCentralMeansZeroRelError)
{
  TH1D* h1= cloneHisto(h1_binning,hName,hName);
  h1->SetDirectory(0);
  h1->SetStats(0);
  //h1->Sumw2();
  for (int ir=0; ir<covM.GetNrows(); ir++) {
    double uncSqr= covM(ir,ir);
    if (uncSqr<0) {
      std::cout << "uncFromCov(covM): negative diagonal entry\n";
      return NULL;
    }
    double val= sqrt(uncSqr);
    if (h1centralVal) {
      double c= h1centralVal->GetBinContent(ir+1);
      if (c!=double(0)) val/=c;
      else { if (zeroCentralMeansZeroRelError) val=0.; else val=9999.; }
    }
    h1->SetBinContent(ir+1, val);
    h1->SetBinError(ir+1, 0);
  }
  return h1;
}

// ---------------------------------------------------------

double totUnc(const std::vector<double> &errV, int returnSqrValue)
{
  double s=0;
  for (unsigned int i=0; i<errV.size(); i++) {
    s+= pow(errV[i],2);
  }
  if (!returnSqrValue) s=sqrt(s);
  return s;
}

// ---------------------------------------------------------

TCanvas *plotCovCorr(TH2D* h2Cov, TString canvName,
		     const PlotCovCorrOpt_t o,
		     TH2D** h2corr_out)
{
  h2Cov->GetXaxis()->SetMoreLogLabels();
  h2Cov->GetXaxis()->SetNoExponent();
  h2Cov->GetYaxis()->SetMoreLogLabels();
  h2Cov->GetYaxis()->SetNoExponent();
  TH2D* h2Corr= cov2corr(h2Cov);
  h2Cov->SetStats(0);

  h2Cov->GetYaxis()->SetTitleOffset(o.yTitleOffset);
  h2Corr->GetYaxis()->SetTitleOffset(o.yTitleOffset);
  h2Cov->GetZaxis()->SetDecimals(true);
  h2Corr->GetZaxis()->SetDecimals(true);

  TCanvas *cx= new TCanvas(canvName,canvName,1200,600);
  cx->Divide(2,1);
  cx->cd(1);
  if (o.gridLines) gPad->SetGrid(1,1);
  if (o.logScaleX) gPad->SetLogx();
  if (o.logScaleY) gPad->SetLogy();
  gPad->SetRightMargin(o.rightMargin);
  gPad->SetLeftMargin(o.leftMargin);
  h2Cov->Draw("COLZ");
  cx->cd(2);
  if (o.gridLines) gPad->SetGrid(1,1);
  if (o.logScaleX) gPad->SetLogx();
  if (o.logScaleY) gPad->SetLogy();
  gPad->SetRightMargin(o.rightMargin);
  gPad->SetLeftMargin(o.leftMargin);
  if (o.autoZRangeCorr) h2Corr->GetZaxis()->SetRangeUser(-1,1);
  h2Corr->Draw("COLZ");
  cx->Update();
  if (h2corr_out) (*h2corr_out)=h2Corr;
  return cx;
}

// ---------------------------------------------------------

TCanvas *plotCovCorr(const TMatrixD &cov, const TH1D *h1_for_axes,
		     TString histoName, TString canvName,
		     const PlotCovCorrOpt_t o,
		     TH2D **h2cov_out, TH2D **h2corr_out)
{
  TH2D *h2cov= convert2histo(cov,h1_for_axes,histoName,histoName);
  if (h2cov_out) *h2cov_out = h2cov;
  return plotCovCorr(h2cov,canvName,o,h2corr_out);
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

void closeCanvases(int itimes)
{
  for (int i=0; i<itimes; i++) {
    TSeqCollection *seq=gROOT->GetListOfCanvases();
    for (int i=0; i<=seq->LastIndex(); i++) {
      TCanvas *c= (TCanvas*) seq->At(i);
      delete c;
    }
    gSystem->ProcessEvents();
  }
}

// ---------------------------------------------------------

void SaveCanvas(TCanvas* canv, TString canvName, TString destDir,
		int correctName)
{
  if (correctName==1) {
    canvName.ReplaceAll(" ","_");
    canvName.ReplaceAll(":","_");
    canvName.ReplaceAll(";","_");
  }

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

int SaveCanvases(TString listOfCanvNames, TString destDir, TFile *fout)
{
  if (listOfCanvNames.Length()==0) return 0;

  if (listOfCanvNames=="ALL") {
    TSeqCollection *seq=gROOT->GetListOfCanvases();
    std::vector<TCanvas*> cV;
    for (int i=0; i<=seq->LastIndex(); i++) {
      TCanvas *c= (TCanvas*) seq->At(i);
      cV.push_back(c);
    }
    SaveCanvases(cV,destDir,fout);
    return 1;
  }

  std::stringstream ss(listOfCanvNames.Data());
  TString cName;
  std::vector<TCanvas*> cV;
  while (!ss.eof()) {
    ss >> cName;
    if (cName.Length()) {
      TCanvas *c= findCanvas(cName);
      if (!c) { std::cout << "\tCanvas " << cName << " was not found\n"; }
      else cV.push_back(c);
    }
  }
  if (cV.size()) SaveCanvases(cV,destDir,fout);
  return int(cV.size());
}

// ---------------------------------------------------------

void writeTimeTag(TFile *fout)
{
  if (fout) fout->cd();
  TObjString timeTag(DayAndTimeTag(0));
  timeTag.Write(timeTag.String());
}

// ---------------------------------------------------------

void writeMsg(TString msg, TFile *fout)
{
  if (fout) fout->cd();
  TObjString str(msg);
  str.Write(str.String());
}

// --------------------------------------------------------------

int getHistosFromCanvas(TCanvas *c, std::vector<TH1D*> *h1V,
			std::vector<TH2D*> *h2V)
{
  if (!h1V && !h2V) {
    std::cout << "getHistosFromCanvas: either h1V or h2V has to be non-null\n";
    return 0;
  }
  TObject *obj=NULL;
  TIter next(c->GetListOfPrimitives());
  int count=0;
  while ((obj=next())) { // assignment!
    //std::cout << "Reading: " << obj->GetName() << "\n";
    if (h1V && obj->InheritsFrom("TH1")) {
      TH1D* h1=(TH1D*)c->GetPrimitive(obj->GetName());
      if (!h1) std::cout << "h1 is null\n";
      else {
	//std::cout << "converted to " << h1->GetName() << "\n";
	h1V->push_back(h1);
	count++;
      }
    }
    if (h2V && obj->InheritsFrom("TH2")) {
      TH2D* h2=(TH2D*)c->GetPrimitive(obj->GetName());
      if (!h2) std::cout << "h2 is null\n";
      else {
	//std::cout << "converted to " << h2->GetName() << "\n";
	h2V->push_back(h2);
	count++;
      }
    }
  }
  return count;
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

void ptrOk(const void *ptr)
{
  if (!ptr) std::cout << "null ptr\n"; else std::cout << "ptr ok\n";
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

