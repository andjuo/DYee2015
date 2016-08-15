#include "DYmm13TeV_eff.h"
#include <TVectorD.h>

// -------------------------------------------------------------
// -------------------------------------------------------------

DYTnPEff_t::DYTnPEff_t(int setElChannel) :
	     h2Eff_RecoID_Data(NULL), h2Eff_Iso_Data(NULL),
	     h2Eff_HLT4p2_Data(NULL), h2Eff_HLT4p3_Data(NULL),
	     h2Eff_RecoID_MC(NULL), h2Eff_Iso_MC(NULL),
	     h2Eff_HLT4p2_MC(NULL), h2Eff_HLT4p3_MC(NULL),
	     h2VReco(), h2VIso(), h2VHLT(),
	     //fDefVal(1.), fDefErr(0.),
	     fEffReco(0.), fEffIso(0.), fEffHlt1(0.), fEffHlt2(0.),
	     fElChannel(setElChannel)
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
	     fEffReco(0.), fEffIso(0.), fEffHlt1(0.), fEffHlt2(0.),
	     fElChannel(0)
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
  fElChannel= e.fElChannel;
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

int DYTnPEff_t::randomize(const DYTnPEff_t &e_errPos,
			  const DYTnPEff_t &e_errNeg, TString tag)
{
  if (!assign(e_errPos,tag)) {
    std::cout << "error in DYTnPEff_t::randomize with tag=" << tag << "\n";
    return 0;
  }
  int bound01=1;
  int res=1;
  res=res && randomizeWithinPosNegErr(e_errPos.h2Eff_RecoID_Data,e_errNeg.h2Eff_RecoID_Data, this->h2Eff_RecoID_Data, bound01);
  res=res && randomizeWithinPosNegErr(e_errPos.h2Eff_Iso_Data,e_errNeg.h2Eff_Iso_Data, this->h2Eff_Iso_Data, bound01);
  res=res && randomizeWithinPosNegErr(e_errPos.h2Eff_HLT4p2_Data,e_errNeg.h2Eff_HLT4p2_Data, this->h2Eff_HLT4p2_Data, bound01);
  res=res && randomizeWithinPosNegErr(e_errPos.h2Eff_HLT4p3_Data,e_errNeg.h2Eff_HLT4p3_Data, this->h2Eff_HLT4p3_Data, bound01);

  res=res && randomizeWithinPosNegErr(e_errPos.h2Eff_RecoID_MC,e_errNeg.h2Eff_RecoID_MC, this->h2Eff_RecoID_MC, bound01);
  res=res && randomizeWithinPosNegErr(e_errPos.h2Eff_Iso_MC,e_errNeg.h2Eff_Iso_MC, this->h2Eff_Iso_MC, bound01);
  res=res && randomizeWithinPosNegErr(e_errPos.h2Eff_HLT4p2_MC,e_errNeg.h2Eff_HLT4p2_MC, this->h2Eff_HLT4p2_MC, bound01);
  res=res && randomizeWithinPosNegErr(e_errPos.h2Eff_HLT4p3_MC,e_errNeg.h2Eff_HLT4p3_MC, this->h2Eff_HLT4p3_MC, bound01);

  return res;
}

// -------------------------------------------------------------

int DYTnPEff_helper_compareArrays(const TArrayD *xa, const TArrayD *ya,
				  const TArrayD *xb, const TArrayD *yb,
				  int axis_difference_detected=0)
{
  int ok=1;

  if ((xa->GetSize()!=xb->GetSize()) || (ya->GetSize()!=yb->GetSize())) {
    ok=0;
  }
  if (axis_difference_detected) {
    std::cout << "  " << xa->GetSize() << "x" << ya->GetSize() << ", "
	      << xb->GetSize() << "x" << yb->GetSize() << "\n";
    if (xa->GetSize()!=xb->GetSize())
      std::cout << " x size difference: " << xa->GetSize() << " vs "
		<< xb->GetSize() << "\n";
    if (ya->GetSize()!=yb->GetSize())
      std::cout << " y size difference: " << ya->GetSize() << " vs "
		<< yb->GetSize() << "\n";
  }

  for (int iOrd=0; iOrd<2; iOrd++) {
    const TArrayD *oa = (iOrd==0) ? xa : ya;
    const TArrayD *ob = (iOrd==0) ? xb : yb;
    if (axis_difference_detected) std::cout << " axis iOrd=" << iOrd << "\n";
    for (int i=0; i < oa->GetSize(); i++) {
      if ( fabs(oa->At(i) - ob->At(i)) > 1e-3 ) {
	ok=0;
	if (axis_difference_detected) {
	  std::cout << "value different: " << oa->At(i) << " vs "
		    << ob->At(i) << "\n";
	}
      }
    }
  }

  if (!ok && !axis_difference_detected) {
    DYTnPEff_helper_compareArrays(xa,ya,xb,yb,1);
  }

  return ok;
}

// -------------------------------------------------------------

int DYTnPEff_t::checkBinning() const
{
  std::vector<const TH2D*> h2Va, h2Vb;
  h2Va.push_back(h2Eff_RecoID_Data);
  h2Va.push_back(h2Eff_Iso_Data);
  h2Va.push_back(h2Eff_HLT4p2_Data);
  h2Va.push_back(h2Eff_HLT4p3_Data);
  h2Vb.push_back(h2Eff_RecoID_MC);
  h2Vb.push_back(h2Eff_Iso_MC);
  h2Vb.push_back(h2Eff_HLT4p2_MC);
  h2Vb.push_back(h2Eff_HLT4p3_MC);

  int ok=1;
  for (unsigned int ih=0; ih < h2Va.size(); ih++) {
    const TArrayD *xa= h2Va[ih]->GetXaxis()->GetXbins();
    const TArrayD *ya= h2Va[ih]->GetYaxis()->GetXbins();
    const TArrayD *xb= h2Vb[ih]->GetXaxis()->GetXbins();
    const TArrayD *yb= h2Vb[ih]->GetYaxis()->GetXbins();
    std::cout << "ih=" << ih << ", h2Va->GetName=" << h2Va[ih]->GetName() << "\n";
    if (!DYTnPEff_helper_compareArrays(xa,ya,xb,yb)) ok=0;
  }

  for (unsigned int ih1=0; ih1<h2Va.size(); ih1++) {
    for (unsigned int ih2=ih1+1; ih2<h2Va.size(); ih2++) {
      const TArrayD *xa= h2Va[ih1]->GetXaxis()->GetXbins();
      const TArrayD *ya= h2Va[ih1]->GetYaxis()->GetXbins();
      const TArrayD *xb= h2Va[ih2]->GetXaxis()->GetXbins();
      const TArrayD *yb= h2Va[ih2]->GetYaxis()->GetXbins();

      std::cout << "ih1=" << ih1 << ", h2Va->GetName=" << h2Va[ih1]->GetName()
		<< " vs ih2=" << ih2 << ", h2Va->GetName="
		<< h2Va[ih2]->GetName() << "\n";
      if (!DYTnPEff_helper_compareArrays(xa,ya,xb,yb)) ok=0;
      else std::cout << " .. compareArrays ok\n";
    }
  }
  return ok;
}
// -------------------------------------------------------------

void DYTnPEff_t::printNumbers() const
{
  std::vector<const TH2D*> h2Va, h2Vb;
  std::vector<TString> labelV;
  h2Va.push_back(h2Eff_RecoID_Data);
  h2Va.push_back(h2Eff_Iso_Data);
  h2Va.push_back(h2Eff_HLT4p2_Data);
  if (!fElChannel) h2Va.push_back(h2Eff_HLT4p3_Data);
  h2Vb.push_back(h2Eff_RecoID_MC);
  h2Vb.push_back(h2Eff_Iso_MC);
  h2Vb.push_back(h2Eff_HLT4p2_MC);
  if (!fElChannel) h2Vb.push_back(h2Eff_HLT4p3_MC);
  if (!fElChannel) {
    labelV.push_back("RECO+ID");
    labelV.push_back("ISO");
    labelV.push_back("HLT v.4p2");
    labelV.push_back("HLT v.4p3");
  }
  else {
    labelV.push_back("RECO");
    labelV.push_back("ID");
    labelV.push_back("HLT");
  }

  std::ofstream fout("tnpEff.tex");
  if (!fout.is_open()) {
    std::cout << "Failed to open latex file\n";
    return;
  }

  fout << "\\documentclass[12pt]{article}\n";
  fout << "\\usepackage{rotating}\n";
  fout << "\\begin{document}\n";

  for (unsigned int ih=0; ih<h2Va.size(); ih++) {
    const TArrayD *xa= h2Va[ih]->GetXaxis()->GetXbins();
    const TArrayD *ya= h2Va[ih]->GetYaxis()->GetXbins();

    TString h2rhoName=Form("h2rho_%d",ih);
    TH2D *h2rho= cloneHisto(h2Va[ih],h2rhoName,h2rhoName);
    h2rho->Divide(h2Vb[ih]);

    //printHisto(h2Va[ih]);

    fout << "%\\begin{sidewaystable}\n";
    fout << "\\centering\n";
    fout << "\\begin{tabular}{|l|";
    for (int ix=0; ix<xa->GetSize()-1; ix++) fout << "c|";
    fout << "}\n";
    fout << "\\multicolumn{" << xa->GetSize() << "}{l}{" << labelV[ih]
	 << "}\\\\\n";
    fout << "\\hline\n";
    fout << " & \\multicolumn{" << (xa->GetSize()-1) << "}{c|}{$\\eta$}\\\\\n";
    fout << Form("\\cline{2-%d}", xa->GetSize()) << "\n";
    fout << "$p_{\\mathrm{T}}$ [GeV] ";
    for (int ix=0; ix<xa->GetSize()-1; ix++ ) {
      fout << Form(" & $[%2.1lf, %2.1lf]$", xa->At(ix),xa->At(ix+1));
    }
    fout << "\n";
    fout << "\\\\\\hline\n";
    for (int iy=0; iy<ya->GetSize()-1; iy++) {
      fout << Form("%2.0lf-%2.0lf", ya->At(iy), ya->At(iy+1));
      for (int ix=0; ix<xa->GetSize()-1; ix++) {
	fout << Form(" & $%4.3lf \\pm %4.3lf$",
		     h2Va[ih]->GetBinContent(ix+1,iy+1),
		     h2Va[ih]->GetBinError(ix+1,iy+1));
	std::cout << " ix=" << ix << " "
		  << h2Va[ih]->GetBinContent(ix+1,iy+1) << "\n";
	printHisto(h2Va[ih]);
      }
      fout << "\\\\\n";

      fout << std::string(5,' ');
      for (int ix=0; ix<xa->GetSize()-1; ix++) {
	fout << Form(" & $%4.3lf \\pm %4.3lf$",
		     h2Vb[ih]->GetBinContent(ix+1,iy+1),
		     h2Vb[ih]->GetBinError(ix+1,iy+1));
      }
      fout << "\\\\\\hline\n";

      fout << std::string(5,' ');
      for (int ix=0; ix<xa->GetSize()-1; ix++) {
	fout << Form(" & $%4.3lf \\pm %4.3lf$",
		     h2rho->GetBinContent(ix+1,iy+1),
		     h2rho->GetBinError(ix+1,iy+1));
      }
      fout << "\\\\\\hline\n";
    }
    //fout << "\\hline\n";
    fout << "\\end{tabular}\n";
    fout << "%\\end{sidewaystable}\n";
    fout << "\n\\vspace{0.5cm}\n\n";
  }

  fout << "\\end{document}\n";
  fout.close();

  return;
}

// -------------------------------------------------------------

void DYTnPEff_t::printEffRatios(const DYTnPEff_t &e, int compareErrs) const
{
  if (!h2Eff_RecoID_Data || !e.h2Eff_RecoID_Data) {
    std::cout << "printEffRatios: initial check failed\n";
    return;
  }
  std::cout << "comparing Effs : e.g. "
	    << this->h2Eff_RecoID_Data->GetName() << " vs "
	    << e.h2Eff_RecoID_Data->GetName() << "\n";
  for (int iv=0; iv<3; iv++) {
    const std::vector<TH2D*> *h2va, *h2vb;
    switch(iv) {
    case 0: h2va= &h2VReco; h2vb= &e.h2VReco; break;
    case 1: h2va= &h2VIso; h2vb= &e.h2VIso; break;
    case 2: h2va= &h2VHLT; h2vb= &e.h2VHLT; break;
    default: std::cout << "code error\n"; return;
    }
    if (!h2va || !h2vb) {
      std::cout << "h2va or h2vb are null\n";
      return;
    }
    for (unsigned int i=0; i<h2va->size(); i++) {
      printRatio(h2va->at(i),h2vb->at(i),0,compareErrs);
    }
  }
}

// -------------------------------------------------------------

int DYTnPEff_t::save(TFile &fout, TString subdirTag)
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
  fout.mkdir("DYTnPEff"+subdirTag);
  fout.cd("DYTnPEff"+subdirTag);
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

  const TString mmNames[4]= { "h2Eff_RecoID", "h2Eff_Iso", "h2Eff_HLT4p2", "h2Eff_HLT4p3" };
  const TString eeNames[4] = { "h2Eff_Reco", "h2Eff_ID", "h2Eff_Trig", "h2Eff_TrigUnneeded" };
  const TString *hnames= (!fElChannel) ? mmNames : eeNames;

  h2Eff_RecoID_Data= loadHisto(fin,subdir+hnames[0]+"_Data"+tag,hnames[0]+"_Data"+tag,1,h2dummy);
  h2Eff_Iso_Data= loadHisto(fin,subdir+hnames[1]+"_Data"+tag,hnames[1]+"_Data"+tag,1,h2dummy);
  h2Eff_HLT4p2_Data= loadHisto(fin,subdir+hnames[2]+"_Data"+tag,hnames[2]+"_Data"+tag,1,h2dummy);
  h2Eff_HLT4p3_Data= loadHisto(fin,subdir+hnames[3]+"_Data"+tag,hnames[3]+"_Data"+tag,1,h2dummy);

  h2Eff_RecoID_MC= loadHisto(fin,subdir+hnames[0]+"_MC"+tag,hnames[0]+"_MC"+tag,1,h2dummy);
  h2Eff_Iso_MC= loadHisto(fin,subdir+hnames[1]+"_MC"+tag,hnames[1]+"_MC"+tag,1,h2dummy);
  h2Eff_HLT4p2_MC= loadHisto(fin,subdir+hnames[2]+"_MC"+tag,hnames[2]+"_MC"+tag,1,h2dummy);
  h2Eff_HLT4p3_MC= loadHisto(fin,subdir+hnames[3]+"_MC"+tag,hnames[3]+"_MC"+tag,1,h2dummy);

  if (subdir.Length()) fin.cd();

  return updateVectors();
}

// -------------------------------------------------------------

// Load Kyeongpil Lee files

int DYTnPEff_t::load_DYMM_TnP(TFile &fin, int version, TString tag)
{
  std::cout << "DYTnPEff_T::load_DYMM_TnP(" << fin.GetName() << ", version="
	    << version << "\n";
  if (!fin.IsOpen()) {
    std::cout << "DYTnPEff_t::load_DYMM_TnP: file is not open\n";
    return 0;
  }
  const int nH2Names=8;
  const TString h2names_ver1[nH2Names] = {
    "h_2D_Eff_RecoID_Data", "h_2D_Eff_Iso_Data",
    "h_2D_Eff_HLTv4p2_Data", "h_2D_Eff_HLTv4p3_Data",
    "h_2D_Eff_RecoID_MC", "h_2D_Eff_Iso_MC",
    "h_2D_Eff_HLTv4p2_MC", "h_2D_Eff_HLTv4p3_MC" };
  const TString h2intNames[nH2Names] = {
    "h2Eff_RecoID_Data", "h2Eff_Iso_Data",
    "h2Eff_HLT4p2_Data", "h2Eff_HLT4p3_Data",
    "h2Eff_RecoID_MC", "h2Eff_Iso_MC",
    "h2Eff_HLT4p2_MC", "h2Eff_HLT4p3_MC" };
  const TString *h2names=h2names_ver1;

  h2Eff_RecoID_Data= loadHisto(fin, h2names[0], h2intNames[0]+tag,1,h2dummy);
  h2Eff_Iso_Data= loadHisto(fin, h2names[1], h2intNames[1]+tag,1,h2dummy);
  h2Eff_HLT4p2_Data=loadHisto(fin, h2names[2], h2intNames[2]+tag,1,h2dummy);
  h2Eff_HLT4p3_Data=loadHisto(fin, h2names[3], h2intNames[3]+tag,1,h2dummy);
  h2Eff_RecoID_MC= loadHisto(fin, h2names[4], h2intNames[4]+tag,1,h2dummy);
  h2Eff_Iso_MC= loadHisto(fin, h2names[5], h2intNames[5]+tag,1,h2dummy);
  h2Eff_HLT4p2_MC= loadHisto(fin, h2names[6], h2intNames[6]+tag,1,h2dummy);
  h2Eff_HLT4p3_MC= loadHisto(fin, h2names[7], h2intNames[7]+tag,1,h2dummy);

  if (!this->updateVectors()) {
    std::cout << "updateVectors failed\n";
    return 0;
  }

  if (!this->ptrsOk()) {
    std::cout << "pointers failed the check\n";
    return 0;
  }

  if (!this->checkBinning()) {
    std::cout << "unforeseen binning problem\n";
    return 0;
  }

  return 1;
}

// -------------------------------------------------------------

TH2D *LoadAsymmErrGraphAsHisto_ver1(TFile &fin, TString base, int posErr,
				    TString tag)
{
  TH2D* h2=(TH2D*)fin.Get("h_2D_" + base);
  h2->Reset();
  h2->SetDirectory(0);
  if (tag.Length()) {
    TString oldName=h2->GetName();
    h2->SetName(oldName+tag);
  }
  for (int iEta=0; iEta<5; iEta++) {
    TString grName= "g_" + base + Form("_EtaBin%d",iEta);
    TGraphAsymmErrors *gr=(TGraphAsymmErrors*)fin.Get(grName);
    if (!gr) {
      std::cout << "failed to load graph " << grName << "\n";
      return NULL;
    }
    TString h1name= grName + tag;
    h1name.ReplaceAll("g_","h1_");
    const int plotIt=0;
    TH1D *h1= convert(gr,h1name,h1name,plotIt,(posErr==1) ? 1 : -1);
    if (!h1) {
      std::cout << "failed to convert gr to h1\n";
      return NULL;
    }
    for (int pTbin=1; pTbin<=h1->GetNbinsX(); pTbin++) {
      h2->SetBinContent(iEta+1, pTbin, h1->GetBinContent(pTbin));
      h2->SetBinError  (iEta+1, pTbin, h1->GetBinError(pTbin));
    }
  }
  return h2;
}

// -------------------------------------------------------------

int DYTnPEff_t::load_DYMM_TnP_asymmEff(TFile &fin, int version, int posErr,
				       TString tag)
{
  int ok=0;
  if (version==1) {
    h2Eff_RecoID_Data= LoadAsymmErrGraphAsHisto_ver1(fin,"Eff_RecoID_Data",posErr,tag);
    h2Eff_Iso_Data= LoadAsymmErrGraphAsHisto_ver1(fin,"Eff_Iso_Data",posErr,tag);
    h2Eff_HLT4p2_Data= LoadAsymmErrGraphAsHisto_ver1(fin,"Eff_HLTv4p2_Data",posErr,tag);
    h2Eff_HLT4p3_Data= LoadAsymmErrGraphAsHisto_ver1(fin,"Eff_HLTv4p3_Data",posErr,tag);

    h2Eff_RecoID_MC= LoadAsymmErrGraphAsHisto_ver1(fin,"Eff_RecoID_MC",posErr,tag);
    h2Eff_Iso_MC= LoadAsymmErrGraphAsHisto_ver1(fin,"Eff_Iso_MC",posErr,tag);
    h2Eff_HLT4p2_MC= LoadAsymmErrGraphAsHisto_ver1(fin,"Eff_HLTv4p2_MC",posErr,tag);
    h2Eff_HLT4p3_MC= LoadAsymmErrGraphAsHisto_ver1(fin,"Eff_HLTv4p3_MC",posErr,tag);
    ok=1;
  }

  if (ok) {
    if (!this->updateVectors()) {
      std::cout << "updateVectors failed\n";
      return 0;
    }
    if (!this->ptrsOk()) {
      std::cout << "pointers failed the check\n";
      return 0;
    }
    if (!this->checkBinning()) {
      std::cout << "unforeseen binning problem\n";
      return 0;
    }
  }
  return ok;
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
  h1rho->GetXaxis()->SetMoreLogLabels();
  h1rho->GetXaxis()->SetNoExponent();
  h1rho->GetXaxis()->SetTitle("M [Gev]");
  //h1rho->GetYaxis()->SetTitle("average event scale factor");
  for (unsigned int im=0; im<fh2ESV.size(); im++) {
    //std::cout << "im=" << im << "\n";
    const TH2D* h2sp= fh2ESV[im];
    if (!h2sp) { HERE("h2sp is null"); return NULL; }
    double sum=0, sumRho=0;

    if (0) {
    for (int ibin1=1; ibin1<=fh2EffBinDef->GetNbinsX(); ibin1++) {
      const double eta1= fh2EffBinDef->GetXaxis()->GetBinCenter(ibin1);
      for (int jbin1=1; jbin1<=fh2EffBinDef->GetNbinsY(); jbin1++) {
	const double pt1= fh2EffBinDef->GetYaxis()->GetBinCenter(jbin1);
	const int fi1= DYtools::FlatIndex(fh2EffBinDef,eta1,pt1,1) + 1;
	//std::cout << "ibin1=" << ibin1 << ", eta1=" << eta1 << ", jbin1=" << jbin1 << ", pt1=" << pt1 << ", fi1=" << fi1 << "\n";
	for (int ibin2=1; ibin2<=fh2EffBinDef->GetNbinsX(); ibin2++) {
	  const double eta2= fh2EffBinDef->GetXaxis()->GetBinCenter(ibin2);
	  for (int jbin2=1; jbin2<=fh2EffBinDef->GetNbinsY(); jbin2++) {
	    const double pt2= fh2EffBinDef->GetYaxis()->GetBinCenter(jbin2);
	    const int fi2= DYtools::FlatIndex(fh2EffBinDef,eta2,pt2,1) + 1;
	    //std::cout << "ibin2=" << ibin2 << ", eta2=" << eta2 << ", jbin2=" << jbin2 << ", pt2=" << pt2 << ", fi2=" << fi2 << "\n";
	    //if (!h2sp) HERE("h2sp is null"); else HERE("h2sp ok");
	    //std::cout << " eSpace bin content=" << h2sp->GetBinContent(fi1,fi2) << std::endl;
	    if (h2sp->GetBinContent(fi1,fi2) == 0) continue;
	    //std::cout << "requesting scale factor\n";
	    double rho= eff.scaleFactor(eta1,pt1,eta2,pt2,hlt4p3);
	    //std::cout << "got rho=" << rho << "\n";
	    sumRho += rho* h2sp->GetBinContent(fi1,fi2);
	    sum += h2sp->GetBinContent(fi1,fi2);
	  }
	}
      }
    }
    }
    else {
    for (int ibin1=1; ibin1<=fh2EffBinDef->GetNbinsX(); ibin1++) {
      for (int jbin1=1; jbin1<=fh2EffBinDef->GetNbinsY(); jbin1++) {
	const int fi1= DYtools::FlatIndex(fh2EffBinDef,ibin1,jbin1,1) + 1;
	//std::cout << "ibin1=" << ibin1 << ", jbin1=" << jbin1 << ", fi1=" << fi1 << "\n";
	for (int ibin2=1; ibin2<=fh2EffBinDef->GetNbinsX(); ibin2++) {
	  for (int jbin2=1; jbin2<=fh2EffBinDef->GetNbinsY(); jbin2++) {
	    const int fi2= DYtools::FlatIndex(fh2EffBinDef,ibin2,jbin2,1) + 1;
	    //std::cout << "ibin2=" << ibin2 << ", jbin2=" << jbin2 << ", fi2=" << fi2 << "\n";
	    //if (!h2sp) HERE("h2sp is null"); else HERE("h2sp ok");
	    //std::cout << " eSpace bin content=" << h2sp->GetBinContent(fi1,fi2) << std::endl;
	    if (h2sp->GetBinContent(fi1,fi2) == 0) continue;
	    //std::cout << "requesting scale factor\n";
	    double rho= eff.scaleFactorIdx(ibin1,jbin1,ibin2,jbin2, hlt4p3);
	    //std::cout << "got rho=" << rho << "\n";
	    sumRho += rho* h2sp->GetBinContent(fi1,fi2);
	    sum += h2sp->GetBinContent(fi1,fi2);
	  }
	}
      }
    }
    }

    double avgRho= (sum==0.) ? 0. : sumRho/sum;
    //std::cout << "im=" << im << ", avgRho=" << avgRho << "\n";
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

std::vector<TH1D*> EventSpace_t::avgAxisValues
        (std::vector<TH2D*> *h2ptSpace, std::vector<TH2D*> *h2etaSpace) const
{

  std::vector<TH1D*> h1V;
  h1V.reserve(4);

  if (int(fh2ESV.size())!=DYtools::nMassBins) {
    std::cout << "EventSpace::avgAxisValues -- nMassBins are different:"
	      << " object has "
	      << fh2ESV.size() << " bins, the code uses nMassBins="
	      << DYtools::nMassBins << "\n";
    return h1V;
  }

 if (h2ptSpace) {
    const TArrayD* ptb= fh2EffBinDef->GetYaxis()->GetXbins();
    if (ptb->GetSize()==0) {
      std::cout << "Yaxis has equal binning. This case is not implemented\n";
      return h1V;
    }
    const Double_t *pt= ptb->GetArray();

    for (int im=0; im<DYtools::nMassBins; im++) {
      TString mStr= DYtools::massStr(im);
      TString hname="h2ptSpace_" + mStr;
      TH2D *h2pt = new TH2D(hname,hname+";p_{1,T} [GeV];p_{2,T} [GeV]",
			    ptb->GetSize()-1,pt,ptb->GetSize()-1,pt);
      h2pt->SetDirectory(0);
      h2pt->SetStats(0);
      h2pt->Sumw2();
      h2ptSpace->push_back(h2pt);
    }
  }

  if (h2etaSpace) {
    const TArrayD* etab= fh2EffBinDef->GetXaxis()->GetXbins();
    if (etab->GetSize()==0) {
      std::cout << "Xaxis has equal binning. This case is not implemented\n";
      return h1V;
    }
    const Double_t *eta= etab->GetArray();

    for (int im=0; im<DYtools::nMassBins; im++) {
      TString mStr= DYtools::massStr(im);
      TString hname="h2etaSpace_" + mStr;
      TH2D *h2eta = new TH2D(hname,hname+";#eta_{1};#eta_{2,T}",
			     etab->GetSize()-1,eta,etab->GetSize()-1,eta);
      h2eta->SetDirectory(0);
      h2eta->SetStats(0);
      h2eta->Sumw2();
      h2etaSpace->push_back(h2eta);
    }
  }

  for (int idxLepton=0; idxLepton<2; idxLepton++) {
    for (int isPt=0; isPt<2; isPt++) {
      TString hname=Form((isPt) ? "avgPt%d" : "avgAbsEta%d", idxLepton);
      TString htitle= hname + TString(";M [GeV];") +
	TString((isPt) ? "<p_{T}> [GeV]" : "|#eta|");
      TH1D* h1avg= new TH1D(hname,htitle,
			    DYtools::nMassBins,DYtools::massBinEdges);
      h1V.push_back(h1avg);
    }
  }


  for (unsigned int im=0; im<fh2ESV.size(); im++) {
    const TH2D* h2sp= fh2ESV[im];
    double sumWPt[2], sumWAbsEta[2];
    double sumWPtSqr[2], sumWAbsEtaSqr[2];
    sumWPt[0]=0.; sumWPt[1]=0.;
    sumWAbsEta[0]=0.; sumWAbsEta[1]=0.;
    sumWPtSqr[0]=0.; sumWPtSqr[1]=0.;
    sumWAbsEtaSqr[0]=0.; sumWAbsEtaSqr[1]=0.;
    double sumW=0.;
    for (int ibin1=1; ibin1<=fh2EffBinDef->GetNbinsX(); ibin1++) {
      const double eta1= fh2EffBinDef->GetXaxis()->GetBinCenter(ibin1);
      for (int jbin1=1; jbin1<=fh2EffBinDef->GetNbinsY(); jbin1++) {
	const double pt1= fh2EffBinDef->GetYaxis()->GetBinCenter(jbin1);
	const int fi1= DYtools::FlatIndex(fh2EffBinDef, ibin1,jbin1,1) + 1;
	for (int ibin2=1; ibin2<=fh2EffBinDef->GetNbinsX(); ibin2++) {
	  const double eta2= fh2EffBinDef->GetXaxis()->GetBinCenter(ibin2);
	  for (int jbin2=1; jbin2<=fh2EffBinDef->GetNbinsY(); jbin2++) {
	    const double pt2= fh2EffBinDef->GetYaxis()->GetBinCenter(jbin2);
	    const int fi2= DYtools::FlatIndex(fh2EffBinDef, ibin2,jbin2,1) + 1;
	    const double w=h2sp->GetBinContent(fi1,fi2);
	    sumW += w;
	    sumWPt[0] += w * pt1;
	    sumWPt[1] += w * pt2;
	    sumWAbsEta[0] += w * fabs(eta1);
	    sumWAbsEta[1] += w * fabs(eta2);
	    sumWPtSqr[0] += w * pt1 * pt1;
	    sumWPtSqr[1] += w * pt2 * pt2;
	    sumWAbsEtaSqr[0] += w * eta1 * eta1;
	    sumWAbsEtaSqr[1] += w * eta2 * eta2;
	    if (w!=double(0)) {
	      std::cout << Form("(%lf,%lf), (%lf,%lf) ",pt1,eta1,pt2,eta2)
			<< " w=" << w << "\n";
	    }
	    if (h2ptSpace) h2ptSpace->at(im)->Fill(pt1,pt2, w);
	    if (h2etaSpace) h2etaSpace->at(im)->Fill(eta1,eta2, w);
	  }
	}
      }
    }

    for (int idxLepton=0, ih=0; idxLepton<2; idxLepton++) {
      for (int isPt=0; isPt<2; isPt++, ih++) {
	double avgVal=0, avgSqrVal=0, avgValVariance=0;
	if (sumW!=double(0)) {
	  if (isPt==1) {
	    avgVal= sumWPt[idxLepton]/sumW;
	    avgSqrVal= sumWPtSqr[idxLepton]/sumW;
	  }
	  else {
	    avgVal= sumWAbsEta[idxLepton]/sumW;
	    avgSqrVal= sumWAbsEtaSqr[idxLepton]/sumW;
	  }
	}
	if (avgSqrVal < avgVal*avgVal) {
	  std::cout << "cannot evaluate variance\n";
	}
	else {
	  avgValVariance = sqrt( avgSqrVal - avgVal * avgVal );
	}

	h1V[ih]->SetBinContent(im+1, avgVal);
	h1V[ih]->SetBinError(im+1, avgValVariance);
      }
    }
  }
  return h1V;
}

// -------------------------------------------------------------

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
