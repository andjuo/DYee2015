#include "DYmm13TeV_eff.h"
#include <TVectorD.h>

// -------------------------------------------------------------
// -------------------------------------------------------------

DYTnPEff_t::DYTnPEff_t() :
	     h2Eff_RecoID_Data(NULL), h2Eff_Iso_Data(NULL),
	     h2Eff_HLT4p2_Data(NULL), h2Eff_HLT4p3_Data(NULL),
	     h2Eff_RecoID_MC(NULL), h2Eff_Iso_MC(NULL),
	     h2Eff_HLT4p2_MC(NULL), h2Eff_HLT4p3_MC(NULL),
	     h2VReco(), h2VIso(), h2VHLT(),
	     //fDefVal(1.), fDefErr(0.),
	     fEffReco(0.), fEffIso(0.), fEffHlt1(0.), fEffHlt2(0.)
{
}

// -------------------------------------------------------------

DYTnPEff_t::DYTnPEff_t(const DYTnPEff_t &e, TString tag) :
	     h2Eff_RecoID_Data(NULL), h2Eff_Iso_Data(NULL),
	     h2Eff_HLT4p2_Data(NULL), h2Eff_HLT4p3_Data(NULL),
	     h2Eff_RecoID_MC(NULL), h2Eff_Iso_MC(NULL),
	     h2Eff_HLT4p2_MC(NULL), h2Eff_HLT4p3_MC(NULL),
	     h2VReco(), h2VIso(), h2VHLT(),
	     //fDefVal(e.fDefVal), fDefErr(e.fDefErr),
	     fEffReco(0.), fEffIso(0.), fEffHlt1(0.), fEffHlt2(0.)
{
  if (!this->assign(e,tag)) {
    std::cout << "error in clone constructor\n";
  }
}

// -------------------------------------------------------------

int DYTnPEff_t::updateVectors()
{
  h2VReco.clear();
  h2VIso.clear();
  h2VHLT.clear();
  h2VReco.push_back(h2Eff_RecoID_Data);
  h2VReco.push_back(h2Eff_RecoID_MC);
  h2VIso.push_back(h2Eff_Iso_Data);
  h2VIso.push_back(h2Eff_Iso_MC);
  h2VHLT.push_back(h2Eff_HLT4p2_Data);
  h2VHLT.push_back(h2Eff_HLT4p2_MC);
  h2VHLT.push_back(h2Eff_HLT4p3_Data);
  h2VHLT.push_back(h2Eff_HLT4p3_MC);
  return ptrsOk();
}

// -------------------------------------------------------------

int DYTnPEff_t::assign(const DYTnPEff_t &e, TString tag)
{
  h2Eff_RecoID_Data= cloneHisto(e.h2Eff_RecoID_Data, "h2Eff_RecoID_Data"+tag,
				"h2Eff_RecoID_Data"+tag);
  h2Eff_Iso_Data= cloneHisto(e.h2Eff_Iso_Data, "h2Eff_Iso_Data"+tag,
			     "h2Eff_Iso_Data"+tag);
  h2Eff_HLT4p2_Data= cloneHisto(e.h2Eff_HLT4p2_Data, "h2Eff_HLT4p2_Data"+tag,
				"h2Eff_HLT4p2_Data"+tag);
  h2Eff_HLT4p3_Data= cloneHisto(e.h2Eff_HLT4p3_Data, "h2Eff_HLT4p3_Data"+tag,
				"h2Eff_HLT4p3_Data"+tag);
  h2Eff_RecoID_MC= cloneHisto(e.h2Eff_RecoID_MC, "h2Eff_RecoID_MC"+tag,
			      "h2Eff_RecoID_MC"+tag);
  h2Eff_Iso_MC= cloneHisto(e.h2Eff_Iso_MC, "h2Eff_Iso_MC"+tag,
			   "h2Eff_Iso_MC"+tag);
  h2Eff_HLT4p2_MC= cloneHisto(e.h2Eff_HLT4p2_MC, "h2Eff_HLT4p2_MC"+tag,
			      "h2Eff_HLT4p2_MC"+tag);
  h2Eff_HLT4p3_MC= cloneHisto(e.h2Eff_HLT4p3_MC, "h2Eff_HLT4p3_MC"+tag,
			      "h2Eff_HLT4p3_MC"+tag);
  return updateVectors();
}

// -------------------------------------------------------------

int DYTnPEff_t::randomize(const DYTnPEff_t &e, TString tag)
{
  if (!assign(e,tag)) {
    std::cout << "error in DYTnPEff_t::randomize with tag=" << tag << "\n";
    return 0;
  }
  int nonNegative=0;
  randomizeWithinErr(e.h2Eff_RecoID_Data, h2Eff_RecoID_Data, nonNegative);
  randomizeWithinErr(e.h2Eff_Iso_Data, h2Eff_Iso_Data, nonNegative);
  randomizeWithinErr(e.h2Eff_HLT4p2_Data, h2Eff_HLT4p2_Data, nonNegative);
  randomizeWithinErr(e.h2Eff_HLT4p3_Data, h2Eff_HLT4p3_Data, nonNegative);
  randomizeWithinErr(e.h2Eff_RecoID_MC, h2Eff_RecoID_MC, nonNegative);
  randomizeWithinErr(e.h2Eff_Iso_MC, h2Eff_Iso_MC, nonNegative);
  randomizeWithinErr(e.h2Eff_HLT4p2_MC, h2Eff_HLT4p2_MC, nonNegative);
  randomizeWithinErr(e.h2Eff_HLT4p3_MC, h2Eff_HLT4p3_MC, nonNegative);
  return 1;
}

// -------------------------------------------------------------

int DYTnPEff_t::save(TFile &fout)
{
  if (!fout.IsOpen()) {
    std::cout << "DYTnPEff_t::save(fout): file is not open\n";
    return 0;
  }
  if (!ptrsOk()) {
    std::cout << "DYTnPEff_T::save(fout): object is not ready\n";
    return 0;
  }
  fout.cd();
  fout.mkdir("DYTnPEff");
  fout.cd("DYTnPEff");
  for (unsigned int i=0; i<h2VReco.size(); i++) h2VReco[i]->Write();
  for (unsigned int i=0; i<h2VIso.size(); i++) h2VIso[i]->Write();
  for (unsigned int i=0; i<h2VHLT.size(); i++) h2VHLT[i]->Write();

  //TVectorD defVals(2);
  //defVals[0]= fDefVal;
  //defVals[1]= fDefErr;
  //defVals.Write("defaultValues");

  fout.cd();
  return 1;
}

// -------------------------------------------------------------

int DYTnPEff_t::load(TFile &fin, TString subdir, TString tag)
{
  if (!fin.IsOpen()) {
    std::cout << "DYTnPEff_t::load(fin): file is not open\n";
    return 0;
  }
  if (subdir.Length() && (subdir[subdir.Length()-1]!='/')) subdir.Append("/");

  h2Eff_RecoID_Data= loadHisto(fin,subdir+"h2Eff_RecoID_Data"+tag,"h2Eff_RecoID_Data"+tag,1,h2dummy);
  h2Eff_Iso_Data= loadHisto(fin,subdir+"h2Eff_Iso_Data"+tag,"h2Eff_Iso_Data"+tag,1,h2dummy);
  h2Eff_HLT4p2_Data= loadHisto(fin,subdir+"h2Eff_HLT4p2_Data"+tag,"h2Eff_HLT4p2_Data"+tag,1,h2dummy);
  h2Eff_HLT4p3_Data= loadHisto(fin,subdir+"h2Eff_HLT4p3_Data"+tag,"h2Eff_HLT4p3_Data"+tag,1,h2dummy);

  h2Eff_RecoID_MC= loadHisto(fin,subdir+"h2Eff_RecoID_MC"+tag,"h2Eff_RecoID_MC"+tag,1,h2dummy);
  h2Eff_Iso_MC= loadHisto(fin,subdir+"h2Eff_Iso_MC"+tag,"h2Eff_Iso_MC"+tag,1,h2dummy);
  h2Eff_HLT4p2_MC= loadHisto(fin,subdir+"h2Eff_HLT4p2_MC"+tag,"h2Eff_HLT4p2_MC"+tag,1,h2dummy);
  h2Eff_HLT4p3_MC= loadHisto(fin,subdir+"h2Eff_HLT4p3_MC"+tag,"h2Eff_HLT4p3_MC"+tag,1,h2dummy);

  if (subdir.Length()) fin.cd();

  return updateVectors();
}

// -------------------------------------------------------------
// -------------------------------------------------------------

EventSpace_t::EventSpace_t(TString setName, TH2D *set_effBinDef) :
  fName(setName), fh2EffBinDef(set_effBinDef),
  fh2ESV()
{
  if (!this->prepareVector()) std::cout << "error in EventSpace_t constructor\n";
}

// -------------------------------------------------------------

EventSpace_t::EventSpace_t(TString setName, const EventSpace_t &es) :
  fName(setName), fh2EffBinDef(NULL),
  fh2ESV()
{
  if (!assign(es)) std::cout << "error in EventSpace_t(ES) constructor\n";
  }

// -------------------------------------------------------------

int EventSpace_t::checkPtrs() const {
  if (!fh2EffBinDef) return 0;
  int res=1;
  for (unsigned int i=0; res && (i<fh2ESV.size()); i++) {
    if (!fh2ESV[i]) res=0;
  }
  return res;
}

// -------------------------------------------------------------

int EventSpace_t::assign(const EventSpace_t &es) {
  TString tagN="_" + fName;
  TString tagT=" " + fName;
  fh2EffBinDef= cloneHisto(es.fh2EffBinDef, es.fh2EffBinDef->GetName() + tagN, es.fh2EffBinDef->GetTitle() + tagT);
  if (fh2ESV.size()) fh2ESV.clear();
  for (unsigned int i=0; i<es.fh2ESV.size(); i++) {
    TH2D *h2= cloneHisto(es.fh2ESV[i], es.fh2ESV[i]->GetName() + tagN, es.fh2ESV[i]->GetTitle() + tagT);
    fh2ESV.push_back(h2);
  }
  int res=checkPtrs();
  if (!res) std::cout << "EventSpace_t::assign: checkPtrs=" << res << "\n";
  return res;
}


// -------------------------------------------------------------

TH1D* EventSpace_t::calculateScaleFactor(const DYTnPEff_t &eff, int hlt4p3,
					 TString hName, TString hTitle) const
{
  TH1D* h1rho= new TH1D(hName,hTitle, DYtools::nMassBins, DYtools::massBinEdges);
  if (!h1rho) {
    std::cout << "cannot create h1rho\n";
    return NULL;
  }
  h1rho->SetDirectory(0);
  h1rho->Sumw2();
  for (unsigned int im=0; im<fh2ESV.size(); im++) {
    std::cout << "im=" << im << "\n";
    const TH2D* h2sp= fh2ESV[im];
    if (!h2sp) HERE("h2sp is null");
    double sum=0, sumRho=0;
    
    for (int ibin1=1; ibin1<=fh2EffBinDef->GetNbinsX(); ibin1++) {
      const double eta1= fh2EffBinDef->GetXaxis()->GetBinCenter(ibin1);
      for (int jbin1=1; jbin1<=fh2EffBinDef->GetNbinsY(); jbin1++) {
	const double pt1= fh2EffBinDef->GetYaxis()->GetBinCenter(jbin1);
	const int fi1= DYtools::FlatIndex(fh2EffBinDef,eta1,pt1,1) + 1;
	std::cout << "ibin1=" << ibin1 << ", eta1=" << eta1 << ", jbin1=" << jbin1 << ", pt1=" << pt1 << ", fi1=" << fi1 << "\n";
	for (int ibin2=1; ibin2<=fh2EffBinDef->GetNbinsX(); ibin2++) {
	  const double eta2= fh2EffBinDef->GetXaxis()->GetBinCenter(ibin2);
	  for (int jbin2=1; jbin2<=fh2EffBinDef->GetNbinsY(); jbin2++) {
	    const double pt2= fh2EffBinDef->GetYaxis()->GetBinCenter(jbin2);
	    const int fi2= DYtools::FlatIndex(fh2EffBinDef,eta2,pt2,1) + 1;
	    std::cout << "ibin2=" << ibin2 << ", eta2=" << eta2 << ", jbin2=" << jbin2 << ", pt2=" << pt2 << ", fi2=" << fi2 << "\n";
	    if (!h2sp) HERE("h2sp is null"); else HERE("h2sp ok");
	    std::cout << "  bin content=" << h2sp->GetBinContent(fi1,fi2) << std::endl;
	    if (h2sp->GetBinContent(fi1,fi2) == 0) continue;
	    std::cout << "requesting scale factor\n";
	    double rho= eff.scaleFactor(eta1,pt1,eta2,pt2,hlt4p3);
	    std::cout << "got rho=" << rho << "\n";
	    sumRho += rho* h2sp->GetBinContent(fi1,fi2);
	    sum += h2sp->GetBinContent(fi1,fi2);
	  }
	}
      }
    }
    double avgRho= (sum==0.) ? 0. : sumRho/sum;
    std::cout << "im=" << im << ", avgRho=" << avgRho << "\n";
    h1rho->SetBinContent( im+1, avgRho );
    h1rho->SetBinError( im+1, 0. );
  }
  return h1rho;
}

// -------------------------------------------------------------

int EventSpace_t::save(TFile &fout)
{
  fout.cd();
  fout.mkdir(fName);
  fout.cd(fName);
  fh2EffBinDef->Write("h2EffBinDef");
  for (unsigned int im=0; im<fh2ESV.size(); im++) {
    TString mStr= DYtools::massStr(im,1);
    fh2ESV[im]->Write("h2ES_" + mStr);
  }
  fout.cd();
  return 1;
}

// -------------------------------------------------------------

int EventSpace_t::load(TFile &fin, TString subdir)
{
  TString tagN="_" + fName;
  TString tagT=" " + fName;
  fh2ESV.clear();
  fh2ESV.reserve(DYtools::nMassBins);

  fin.cd();
  //fin.cd(subdir);
  if (subdir.Length() && (subdir[subdir.Length()-1]!='/')) subdir.Append("/");
  fh2EffBinDef=loadHisto(fin,subdir+"h2EffBinDef","h2EffBinDef" + tagN,1,h2dummy);
  for (int im=0; im<DYtools::nMassBins; im++) {
    TString mStr= DYtools::massStr(im,1);
    TH2D* h2=loadHisto(fin, subdir+"h2ES_" + mStr, "h2ES_" + mStr + tagN,1,h2dummy);
    fh2ESV.push_back(h2);
  }
  fin.cd();
  return checkPtrs();
}


// -------------------------------------------------------------

int EventSpace_t::prepareVector() {
  if (fh2ESV.size()) fh2ESV.clear();
  fh2ESV.reserve(DYtools::nMassBins);

  TString tag=fName;
  tag.ReplaceAll(" ","_");
  int res=1;
  for (int im=0; res && (im<DYtools::nMassBins); im++) {
    TString mStr= DYtools::massStr(im);
    TString hTitle=TString("Event space for ") + mStr + TString(" GeV (") + fName + (");(eta,pt) flat index 1;(eta,pt) flat index 2");
    mStr.ReplaceAll("-","_");
    TH2D* h2= new TH2D("h2ES"+mStr+tag,hTitle,
	   DYtools::EtaPtFIMax+1, -0.5, DYtools::EtaPtFIMax+0.5,
	   DYtools::EtaPtFIMax+1, -0.5, DYtools::EtaPtFIMax+0.5);
    if (!h2) { res=0; continue; }
    h2->SetDirectory(0);
    h2->Sumw2();
    fh2ESV.push_back(h2);
  }
  return res;
}

// -------------------------------------------------------------
// -------------------------------------------------------------
