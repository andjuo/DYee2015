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

DYTnPEff_t::~DYTnPEff_t()
{
  if (h2Eff_RecoID_Data) delete h2Eff_RecoID_Data;
  if (h2Eff_Iso_Data) delete h2Eff_Iso_Data;
  if (h2Eff_HLT4p2_Data) delete h2Eff_HLT4p2_Data;
  if (h2Eff_HLT4p3_Data) delete h2Eff_HLT4p3_Data;
  if (h2Eff_RecoID_MC) delete h2Eff_RecoID_MC;
  if (h2Eff_Iso_MC) delete h2Eff_Iso_MC;
  if (h2Eff_HLT4p2_MC) delete h2Eff_HLT4p2_MC;
  if (h2Eff_HLT4p3_MC) delete h2Eff_HLT4p3_MC;
  h2VReco.clear();
  h2VIso.clear();
  h2VHLT.clear();
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

void DYTnPEff_t::removeError()
{
  for (int ikind=0; ikind<3; ikind++) {
    std::vector<TH2D*> *h2V= this->h2vecPtr(ikind);
    for (unsigned int i=0; i<h2V->size(); i++) {
      removeErrorH2((*h2V)[i]);
    }
  }
}

// -------------------------------------------------------------

void DYTnPEff_t::setError(const DYTnPEff_t &e)
{
  for (int ikind=0; ikind<3; ikind++) {
    std::vector<TH2D*> *h2V= this->h2vecPtr(ikind);
    const std::vector<TH2D*> *h2srcV= e.h2vecPtr(ikind);
    for (unsigned int i=0; i<h2V->size(); i++) {
      ::setError((*h2V)[i],(*h2srcV)[i]);
    }
  }
}

// -------------------------------------------------------------

void DYTnPEff_t::resetAll()
{
  for (int ikind=0; ikind<3; ikind++) {
    std::vector<TH2D*> *h2V= this->h2vecPtr(ikind);
    for (unsigned int i=0; i<h2V->size(); i++) {
      (*h2V)[i]->Reset();
    }
  }
}

// -------------------------------------------------------------

double DYTnPEff_t::totalEffIdx_misc(int ibin1, int jbin1,
				    int ibin2, int jbin2,
				    int mc, int miscFlag) const
{
  int hlt4p3=1;
  if (miscFlag==_misc_hlt4p2_plain) hlt4p3=0;
  else if (miscFlag==_misc_hlt4p3_plain) hlt4p3=1;
  else if ((miscFlag==_misc_hlt4p2_leadLep) ||
	   (miscFlag==_misc_hlt4p3_leadLep)) {
    // only leading electron determines the efficiency
    fEffReco= h2VReco[mc]->GetBinContent(ibin1,jbin1);
    fEffIso=  h2VIso[mc]->GetBinContent(ibin1,jbin1);
    hlt4p3= (miscFlag==_misc_hlt4p3_leadLep) ? 1 : 0;
    fEffHlt1= h2VHLT[mc+2*hlt4p3]->GetBinContent(ibin1,jbin1);
    fEffHlt2= 1.;
    double eff_lead= fEffReco *fEffIso * fEffHlt1;
    if (0) {
      std::cout << "totalEffIdx_misc(" << ibin1 << "," << jbin1 << "; "
		<< ibin2 << "," << jbin2 << "; miscFlag=" << miscFlag
		<< ")= "
		<< h2VReco[mc]->GetBinContent(ibin1,jbin1) << " x "
		<< h2VIso[mc]->GetBinContent(ibin1,jbin1) << " x "
		<< h2VHLT[mc+2*hlt4p3]->GetBinContent(ibin1,jbin1) << " x "
		<< "tot=" << eff_lead << "\n";
    }
    return eff_lead;
  }
  else {
    std::cout << "totalEffIdx_misc is not ready for miscFlag="
	      << int(miscFlag) << "\n";
    return 0.;
  }

  fEffReco=
    h2VReco[mc]->GetBinContent(ibin1,jbin1) *
    h2VReco[mc]->GetBinContent(ibin2,jbin2);
  fEffIso=
    h2VIso[mc]->GetBinContent(ibin1,jbin1) *
    h2VIso[mc]->GetBinContent(ibin2,jbin2);
  fEffHlt1= h2VHLT[mc+2*hlt4p3]->GetBinContent(ibin1,jbin1);
  fEffHlt2= h2VHLT[mc+2*hlt4p3]->GetBinContent(ibin2,jbin2);
  double effHlt= fEffHlt1*fEffHlt2;  /// Plain HLT
  double eff= fEffReco *fEffIso * effHlt;
  if (0) {
    std::cout << "totalEffIdx_misc(" << ibin1 << "," << jbin1 << "; "
	      << ibin2 << "," << jbin2 << "; "
	      << h2VReco[mc]->GetBinContent(ibin1,jbin1) << " x "
	      << h2VReco[mc]->GetBinContent(ibin2,jbin2) << "; "
	      << h2VIso[mc]->GetBinContent(ibin1,jbin1) << " x "
	      << h2VIso[mc]->GetBinContent(ibin2,jbin2) << "; "
	      << h2VHLT[mc+2*hlt4p3]->GetBinContent(ibin1,jbin1) << " x "
	      << h2VHLT[mc+2*hlt4p3]->GetBinContent(ibin2,jbin2) << "; "
	      << "tot=" << eff << "\n";
  }
  return eff;
}

// -------------------------------------------------------------

double DYTnPEff_t::totalEff_misc(double eta1, double pt1,
				 double eta2, double pt2,
				 int mc, int miscFlag) const
{
  int ibin1=h2Eff_RecoID_Data->GetXaxis()->FindBin(eta1);
  int ibin2=h2Eff_RecoID_Data->GetXaxis()->FindBin(eta2);
  int jbin1=h2Eff_RecoID_Data->GetYaxis()->FindBin(pt1);
  int jbin2=h2Eff_RecoID_Data->GetYaxis()->FindBin(pt2);
  if ((ibin1==0) || (jbin1==0) || (ibin2==0) || (jbin2==0)) return 0.;
  if (ibin1==h2Eff_RecoID_Data->GetNbinsX()+1) ibin1--;
  if (ibin2==h2Eff_RecoID_Data->GetNbinsX()+1) ibin2--;
  if (jbin1==h2Eff_RecoID_Data->GetNbinsY()+1) jbin1--;
  if (jbin2==h2Eff_RecoID_Data->GetNbinsY()+1) jbin2--;
  return totalEffIdx_misc(ibin1,jbin1,ibin2,jbin2, mc, miscFlag);
}

// -------------------------------------------------------------

double DYTnPEff_t::scaleFactor_misc(double eta1, double pt1,
				    double eta2, double pt2,
				    int miscFlag) const
{
  //std::cout << "scaleFactor_misc( " << Form(" (%2.2lf,%2.0lf)", eta1,pt1)
  //	      << Form(" (%2.2lf,%2.0lf)", eta2,pt2)
  //	      << " miscFlag=" << miscFlag << "\n";
  double effData= totalEff_misc(eta1,pt1,eta2,pt2, 0,miscFlag);
  //std::cout << " totalEffData=" << effData << "\n";
  double effMC= totalEff_misc(eta1,pt1,eta2,pt2, 1,miscFlag);
  //std::cout << " totalEffMC=" << effMC << "\n";
  double rho= (effMC==0.) ? 0. : effData/effMC;
  //std::cout << " rho=" << rho << "\n";
  return rho;
}

// -------------------------------------------------------------

double DYTnPEff_t::scaleFactorIdx_misc(int ibin1, int jbin1,
				       int ibin2, int jbin2,
				       int miscFlag) const
{
  //std::cout << "scaleFactorIdx_misc( " << Form(" (%2d,%2d)", ibin1,jbin1)
  //	      << Form(" (%2d,%2d)", ibin2,jbin2)
  //	      << " miscFlag=" << miscFlag << "\n";
  double effData= totalEffIdx_misc(ibin1,jbin1,ibin2,jbin2, 0,miscFlag);
  //std::cout << " totalEffData=" << effData << "\n";
  double effMC= totalEffIdx_misc(ibin1,jbin1,ibin2,jbin2, 1,miscFlag);
  //std::cout << " totalEffMC=" << effMC << "\n";
  double rho= (effMC==0.) ? 0. : effData/effMC;
  //std::cout << " rho=" << rho << "\n";
  return rho;
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

int DYTnPEff_t::assign_DiffAsUnc(const DYTnPEff_t &e, int relative)
{
  int res=assignDiffAsUnc(h2Eff_RecoID_Data, e.h2Eff_RecoID_Data, relative);
  if (res) res=assignDiffAsUnc(h2Eff_Iso_Data, e.h2Eff_Iso_Data, relative);
  if (res) res=assignDiffAsUnc(h2Eff_HLT4p2_Data, e.h2Eff_HLT4p2_Data, relative);
  if (res) res=assignDiffAsUnc(h2Eff_HLT4p3_Data, e.h2Eff_HLT4p3_Data, relative);
  if (res) res=assignDiffAsUnc(h2Eff_RecoID_MC, e.h2Eff_RecoID_MC, relative);
  if (res) res=assignDiffAsUnc(h2Eff_Iso_MC, e.h2Eff_Iso_MC, relative);
  if (res) res=assignDiffAsUnc(h2Eff_HLT4p2_MC, e.h2Eff_HLT4p2_MC, relative);
  if (res) res=assignDiffAsUnc(h2Eff_HLT4p3_MC, e.h2Eff_HLT4p3_MC, relative);
  fElChannel= e.fElChannel;
  if (res) res= updateVectors();
  return res;
}

// -------------------------------------------------------------

int DYTnPEff_t::randomize(const DYTnPEff_t &e, TString tag, int systematic)
{
  if (!assign(e,tag)) {
    std::cout << "error in DYTnPEff_t::randomize with tag=" << tag << "\n";
    return 0;
  }
  int nonNegative=0;
  //std::cout << "randomize within error : "; printHisto(e.h2Eff_RecoID_Data);
  randomizeWithinErr(e.h2Eff_RecoID_Data, h2Eff_RecoID_Data, nonNegative,0,systematic);
  randomizeWithinErr(e.h2Eff_Iso_Data, h2Eff_Iso_Data, nonNegative,0,systematic);
  randomizeWithinErr(e.h2Eff_HLT4p2_Data, h2Eff_HLT4p2_Data, nonNegative,0,systematic);
  randomizeWithinErr(e.h2Eff_HLT4p3_Data, h2Eff_HLT4p3_Data, nonNegative,0,systematic);
  randomizeWithinErr(e.h2Eff_RecoID_MC, h2Eff_RecoID_MC, nonNegative,0,systematic);
  randomizeWithinErr(e.h2Eff_Iso_MC, h2Eff_Iso_MC, nonNegative,0,systematic);
  randomizeWithinErr(e.h2Eff_HLT4p2_MC, h2Eff_HLT4p2_MC, nonNegative,0,systematic);
  randomizeWithinErr(e.h2Eff_HLT4p3_MC, h2Eff_HLT4p3_MC, nonNegative,0,systematic);
  return 1;
}

// -------------------------------------------------------------

int DYTnPEff_t::randomizeRelErr(const DYTnPEff_t &e0, const DYTnPEff_t &e,
				TString tag, int systematic)
{
  if (!assign(e,tag)) {
    std::cout << "error in DYTnPEff_t::randomize with tag=" << tag << "\n";
    return 0;
  }
  int nonNegative=0;
  randomizeWithinRelErrSetCentral(e0.h2Eff_RecoID_Data, e.h2Eff_RecoID_Data, h2Eff_RecoID_Data, nonNegative,systematic);
  randomizeWithinRelErrSetCentral(e0.h2Eff_Iso_Data, e.h2Eff_Iso_Data, h2Eff_Iso_Data, nonNegative,systematic);
  randomizeWithinRelErrSetCentral(e0.h2Eff_HLT4p2_Data, e.h2Eff_HLT4p2_Data, h2Eff_HLT4p2_Data, nonNegative,systematic);
  randomizeWithinRelErrSetCentral(e0.h2Eff_HLT4p3_Data, e.h2Eff_HLT4p3_Data, h2Eff_HLT4p3_Data, nonNegative,systematic);
  randomizeWithinRelErrSetCentral(e0.h2Eff_RecoID_MC, e.h2Eff_RecoID_MC, h2Eff_RecoID_MC, nonNegative,systematic);
  randomizeWithinRelErrSetCentral(e0.h2Eff_Iso_MC, e.h2Eff_Iso_MC, h2Eff_Iso_MC, nonNegative,systematic);
  randomizeWithinRelErrSetCentral(e0.h2Eff_HLT4p2_MC, e.h2Eff_HLT4p2_MC, h2Eff_HLT4p2_MC, nonNegative,systematic);
  randomizeWithinRelErrSetCentral(e0.h2Eff_HLT4p3_MC, e.h2Eff_HLT4p3_MC, h2Eff_HLT4p3_MC, nonNegative,systematic);
  return 1;
}

// -------------------------------------------------------------

int DYTnPEff_t::randomize(const DYTnPEff_t &e_errPos,
			  const DYTnPEff_t &e_errNeg, TString tag,
			  int systematic)
{
  if (!assign(e_errPos,tag)) {
    std::cout << "error in DYTnPEff_t::randomize with tag=" << tag << "\n";
    return 0;
  }
  int bound01=1;
  int res=1;
  res=res && randomizeWithinPosNegErr(e_errPos.h2Eff_RecoID_Data,e_errNeg.h2Eff_RecoID_Data, this->h2Eff_RecoID_Data, bound01, systematic);
  res=res && randomizeWithinPosNegErr(e_errPos.h2Eff_Iso_Data,e_errNeg.h2Eff_Iso_Data, this->h2Eff_Iso_Data, bound01, systematic);
  res=res && randomizeWithinPosNegErr(e_errPos.h2Eff_HLT4p2_Data,e_errNeg.h2Eff_HLT4p2_Data, this->h2Eff_HLT4p2_Data, bound01, systematic);
  res=res && randomizeWithinPosNegErr(e_errPos.h2Eff_HLT4p3_Data,e_errNeg.h2Eff_HLT4p3_Data, this->h2Eff_HLT4p3_Data, bound01, systematic);

  res=res && randomizeWithinPosNegErr(e_errPos.h2Eff_RecoID_MC,e_errNeg.h2Eff_RecoID_MC, this->h2Eff_RecoID_MC, bound01, systematic);
  res=res && randomizeWithinPosNegErr(e_errPos.h2Eff_Iso_MC,e_errNeg.h2Eff_Iso_MC, this->h2Eff_Iso_MC, bound01, systematic);
  res=res && randomizeWithinPosNegErr(e_errPos.h2Eff_HLT4p2_MC,e_errNeg.h2Eff_HLT4p2_MC, this->h2Eff_HLT4p2_MC, bound01, systematic);
  res=res && randomizeWithinPosNegErr(e_errPos.h2Eff_HLT4p3_MC,e_errNeg.h2Eff_HLT4p3_MC, this->h2Eff_HLT4p3_MC, bound01, systematic);

  return res;
}

// -------------------------------------------------------------

int DYTnPEff_t::multiply(const DYTnPEff_t &e)
{
  h2Eff_RecoID_Data->Multiply(e.h2Eff_RecoID_Data);
  h2Eff_Iso_Data->Multiply(e.h2Eff_Iso_Data);
  h2Eff_HLT4p2_Data->Multiply(e.h2Eff_HLT4p2_Data);
  h2Eff_HLT4p3_Data->Multiply(e.h2Eff_HLT4p3_Data);
  h2Eff_RecoID_MC->Multiply(e.h2Eff_RecoID_MC);
  h2Eff_Iso_MC->Multiply(e.h2Eff_Iso_MC);
  h2Eff_HLT4p2_MC->Multiply(e.h2Eff_HLT4p2_MC);
  h2Eff_HLT4p3_MC->Multiply(e.h2Eff_HLT4p3_MC);
  return 1;
}

// -------------------------------------------------------------

int DYTnPEff_t::add(const DYTnPEff_t &e)
{
  h2Eff_RecoID_Data->Add(e.h2Eff_RecoID_Data);
  h2Eff_Iso_Data->Add(e.h2Eff_Iso_Data);
  h2Eff_HLT4p2_Data->Add(e.h2Eff_HLT4p2_Data);
  h2Eff_HLT4p3_Data->Add(e.h2Eff_HLT4p3_Data);
  h2Eff_RecoID_MC->Add(e.h2Eff_RecoID_MC);
  h2Eff_Iso_MC->Add(e.h2Eff_Iso_MC);
  h2Eff_HLT4p2_MC->Add(e.h2Eff_HLT4p2_MC);
  h2Eff_HLT4p3_MC->Add(e.h2Eff_HLT4p3_MC);
  return 1;
}

// -------------------------------------------------------------

int DYTnPEff_helper_includeRelUnc(TH2D *h2dest,
				  const TH2D *h2central, const TH2D *h2relUnc)
{
  for (int ibin=1; ibin<=h2dest->GetNbinsX(); ibin++) {
    for (int jbin=1; jbin<=h2dest->GetNbinsY(); jbin++) {
      double e= h2dest->GetBinError(ibin,jbin);
      double ve= h2central->GetBinContent(ibin,jbin) *
	h2relUnc->GetBinError(ibin,jbin);
      h2dest->SetBinError(ibin,jbin, sqrt(e*e + ve*ve));
    }
  }
  return 1;
}

// -------------------------------------------------------------

int DYTnPEff_t::includeRelUnc(const DYTnPEff_t &e0, const DYTnPEff_t &eRelErr)
{
  DYTnPEff_helper_includeRelUnc(h2Eff_RecoID_Data,e0.h2Eff_RecoID_Data,eRelErr.h2Eff_RecoID_Data);
  DYTnPEff_helper_includeRelUnc(h2Eff_Iso_Data,e0.h2Eff_Iso_Data,eRelErr.h2Eff_Iso_Data);
  DYTnPEff_helper_includeRelUnc(h2Eff_HLT4p2_Data,e0.h2Eff_HLT4p2_Data,eRelErr.h2Eff_HLT4p2_Data);
  DYTnPEff_helper_includeRelUnc(h2Eff_HLT4p3_Data,e0.h2Eff_HLT4p3_Data,eRelErr.h2Eff_HLT4p3_Data);
  DYTnPEff_helper_includeRelUnc(h2Eff_RecoID_MC,e0.h2Eff_RecoID_MC,eRelErr.h2Eff_RecoID_MC);
  DYTnPEff_helper_includeRelUnc(h2Eff_Iso_MC,e0.h2Eff_Iso_MC,eRelErr.h2Eff_Iso_MC);
  DYTnPEff_helper_includeRelUnc(h2Eff_HLT4p2_MC,e0.h2Eff_HLT4p2_MC,eRelErr.h2Eff_HLT4p2_MC);
  DYTnPEff_helper_includeRelUnc(h2Eff_HLT4p3_MC,e0.h2Eff_HLT4p3_MC,eRelErr.h2Eff_HLT4p3_MC);
  return 1;
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

void DYTnPEff_t::printNumbers(TString fileTag) const
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

  TString fname="tnpEff.tex";
  if (fileTag.Length()) fname.ReplaceAll(".tex",fileTag+".tex");
  std::ofstream fout(fname.Data());
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

void DYTnPEff_t::displayAll() const
{
  int logScale=1;
#ifdef def_fewMassBins
  logScale=0;
#endif
  PlotCovCorrOpt_t opt(0,1,logScale);

  if (0) {
    for (int iv=0; iv<3; iv++) {
      const std::vector<TH2D*> *h2=
	(iv==0) ? &h2VReco : ((iv==1) ? &h2VIso : &h2VHLT);
      for (unsigned int ih=0; ih<h2->size(); ih++) {
	TString cName= TString("canv_") + (*h2)[ih]->GetName();
	plotHisto((*h2)[ih],cName,opt);
	TH2D *h2err= errorAsCentral((*h2)[ih],0);
	plotHisto(h2err,cName + "_err",opt);
      }
    }
  }
  else {
    TString cName1= TString("canv_") + h2Eff_RecoID_Data->GetName();
    plotHisto(h2Eff_RecoID_Data,cName1,opt);
    TH2D *h2err1= errorAsCentral(h2Eff_RecoID_Data,0);
    plotHisto(h2err1,cName1 + "_err",opt);

    TString cName2= TString("canv_") + h2Eff_Iso_MC->GetName();
    plotHisto(h2Eff_Iso_MC,cName2,opt);
    TH2D *h2err2= errorAsCentral(h2Eff_Iso_MC,0);
    plotHisto(h2err2,cName2 + "_err",opt);
  }
}

// -------------------------------------------------------------

void DYTnPEff_t::displayEffTot(int data_or_mc, TString tag, int hlt4p3,
			       int miscFlag) const
{
  if ((hlt4p3==-1) && (miscFlag==_misc_empty)) {
    std::cout << "displayEffTot: changed hlt4p3 from -1 to 0\n";
    hlt4p3=0;
  }
  const TH2D *h2def= h2Eff_RecoID_Data;
  int fiMax=
    DYtools::FlatIndex(h2def,h2def->GetNbinsX(),h2def->GetNbinsY(),0) + 1;
  std::cout << "fiMax=" << fiMax << "\n";
  TH2D *h2defTot= new TH2D("h2defTot"+tag,"h2defTot;fi1;fi2",
			   fiMax,-0.5,fiMax-0.5, fiMax,-0.5,fiMax-0.5);
  h2defTot->SetDirectory(0);
  h2defTot->SetStats(0);
  TH2D *h2MC=NULL, *h2Data=NULL;
  for (int isMC=0; isMC<2; isMC++) {
    if ( ((1<<isMC) & data_or_mc) != 0 ) {
      TString hName= h2defTot->GetName() + tag + Form("_isMC%d",isMC);
      if (miscFlag!=_misc_empty) hName.Append(Form("_misc%d",miscFlag));
      TH2D *h2= cloneHisto(h2defTot, hName,hName);
      if (isMC) h2MC=h2; else h2Data=h2;
      for (int ibin1=1; ibin1<=h2def->GetNbinsX(); ibin1++) {
	for (int jbin1=1; jbin1<=h2def->GetNbinsY(); jbin1++) {
	  int fi1= DYtools::FlatIndex(h2def,ibin1,jbin1,0);
	  for (int ibin2=1; ibin2<=h2def->GetNbinsX(); ibin2++) {
	    for (int jbin2=1; jbin2<=h2def->GetNbinsY(); jbin2++) {
	      int fi2= DYtools::FlatIndex(h2def,ibin2,jbin2,0);
	      double eff= (miscFlag!=_misc_empty) ?
		totalEffIdx_misc(ibin1,jbin1,ibin2,jbin2,isMC,miscFlag) :
		totalEffIdx(ibin1,jbin1,ibin2,jbin2,isMC,hlt4p3);
	      h2->SetBinContent(fi1+1,fi2+1, eff);
	      h2->SetBinError  (fi1+1,fi2+1, 0.);
	    }
	  }
	}
      }
    }
  }
  TString cStart= "canv_";
  if (h2MC) plotHisto(h2MC, cStart + h2MC->GetName());
  if (h2Data) plotHisto(h2Data, cStart + h2Data->GetName());
  if ( (data_or_mc & 4) != 0 ) {
    TString hName= h2Data->GetName() + TString("__DIV__") + h2MC->GetName();
    TH2D* h2SF= cloneHisto(h2Data,hName,hName);
    h2SF->Divide(h2MC);
    plotHisto(h2SF,cStart+h2SF->GetName());
  }
}

// -------------------------------------------------------------

void DYTnPEff_t::listNumbers() const
{
  std::cout << "\n" << string(20,'-') << "\n";
  printHisto(h2Eff_RecoID_Data);
  /*
  printHisto(h2Eff_RecoID_MC);
  printHisto(h2Eff_Iso_Data);
  */
  printHisto(h2Eff_Iso_MC);
  //printHisto(h2Eff_HLT4p2_Data);
  //printHisto(h2Eff_HLT4p2_MC);
  /*
  printHisto(h2Eff_HLT4p3_Data);
  printHisto(h2Eff_HLT4p3_MC);
  */
  std::cout << string(20,'-') << "\n";
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
  const TString eeNamesV2[4]= { "h_2D_Eff_RECO", "h_2D_Eff_ID", "h_2D_Eff_Trig", "h_2D_Eff_TrigUnneeded" };
  const TString *hnames = mmNames;
  if (fElChannel==1) hnames=eeNames;
  else if (fElChannel==2) hnames=eeNamesV2;

  h2Eff_RecoID_Data= loadHisto(fin,subdir+hnames[0]+"_Data"+tag,hnames[0]+"_Data"+tag,1,h2dummy);
  h2Eff_Iso_Data= loadHisto(fin,subdir+hnames[1]+"_Data"+tag,hnames[1]+"_Data"+tag,1,h2dummy);
  h2Eff_HLT4p2_Data= loadHisto(fin,subdir+hnames[2]+"_Data"+tag,hnames[2]+"_Data"+tag,1,h2dummy);
  if (!fElChannel) h2Eff_HLT4p3_Data= loadHisto(fin,subdir+hnames[3]+"_Data"+tag,hnames[3]+"_Data"+tag,1,h2dummy);
  else h2Eff_HLT4p3_Data=cloneHisto(h2Eff_HLT4p2_Data,hnames[3]+"_Data"+tag,h2Eff_HLT4p2_Data->GetTitle());

  h2Eff_RecoID_MC= loadHisto(fin,subdir+hnames[0]+"_MC"+tag,hnames[0]+"_MC"+tag,1,h2dummy);
  h2Eff_Iso_MC= loadHisto(fin,subdir+hnames[1]+"_MC"+tag,hnames[1]+"_MC"+tag,1,h2dummy);
  h2Eff_HLT4p2_MC= loadHisto(fin,subdir+hnames[2]+"_MC"+tag,hnames[2]+"_MC"+tag,1,h2dummy);
  if (!fElChannel) h2Eff_HLT4p3_MC= loadHisto(fin,subdir+hnames[3]+"_MC"+tag,hnames[3]+"_MC"+tag,1,h2dummy);
  else h2Eff_HLT4p3_MC= cloneHisto(h2Eff_HLT4p2_MC,hnames[3]+"_MC"+tag,h2Eff_HLT4p2_MC->GetTitle());

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

DYTnPEffColl_t::DYTnPEffColl_t(int setElChannel) :
  fTnPEff(setElChannel),
  fTnPEffSrcV(), fTnPEffSystV(),
  fElChannel(setElChannel)
{}


// -------------------------------------------------------------

DYTnPEffColl_t::DYTnPEffColl_t(const DYTnPEffColl_t &e, TString tag) :
  fTnPEff(e.isElChannel()),
  fTnPEffSrcV(), fTnPEffSystV(),
  fElChannel(e.isElChannel())
{
  if (!this->assign(e,tag)) {
    std::cout << "DYTnPEffColl_t:: failure in constructor\n";
  }
}

// -------------------------------------------------------------

DYTnPEff_t* DYTnPEffColl_t::getTnPWithSystUnc(int idx) const
{
  if ((idx<0) || (idx>=int(fTnPEffSystV.size()))) {
    std::cout << "getTnPWithSystUnc(idx=" << idx << ") idx is too large\n";
    return NULL;
  }
  DYTnPEff_t *e= new DYTnPEff_t(fTnPEff,Form("_systUnc%d",idx));
  e->setError(*fTnPEffSystV[idx]);
  return e;
}

// -------------------------------------------------------------

DYTnPEff_t* DYTnPEffColl_t::getTnPWithTotUnc() const
{
  DYTnPEff_t *e= new DYTnPEff_t(fTnPEff,"_totUnc");
  for (unsigned int i=0; i<fTnPEffSystV.size(); i++) {
    e->add(*fTnPEffSystV[i]);
  }
  return e;
}

// -------------------------------------------------------------

int DYTnPEffColl_t::ptrsOk() const
{
  int ok= (fTnPEff.ptrsOk() && (fTnPEffSrcV.size()==fTnPEffSystV.size())) ? 1:0;
  for (unsigned int i=0; ok && (i<fTnPEffSrcV.size()); i++) {
    ok= fTnPEffSrcV[i]->ptrsOk();
  }
  for (unsigned int i=0; ok && (i<fTnPEffSystV.size()); i++) {
    ok= fTnPEffSystV[i]->ptrsOk();
  }
  return ok;
}

// -------------------------------------------------------------

int DYTnPEffColl_t::assign(const DYTnPEffColl_t &coll, TString tag)
{
  int ok= fTnPEff.assign(coll.fTnPEff,tag);
  for (unsigned int i=0; ok && (i<coll.fTnPEffSrcV.size()); i++) {
    DYTnPEff_t *e=new DYTnPEff_t(*coll.fTnPEffSrcV[i],tag);
    if (!e) ok=0;
    else { ok= e->ptrsOk(); fTnPEffSrcV.push_back(e); }
  }
  for (unsigned int i=0; ok && (i<coll.fTnPEffSystV.size()); i++) {
    DYTnPEff_t *e=new DYTnPEff_t(*coll.fTnPEffSystV[i],tag);
    if (!e) ok=0;
    else { ok=e->ptrsOk(); fTnPEffSystV.push_back(e); }
  }
  fElChannel= coll.fElChannel;
  return ok;
}

// -------------------------------------------------------------

int DYTnPEffColl_t::assignStatErr(const DYTnPEff_t &e, TString tag)
{ return fTnPEff.assign(e,tag); }

// -------------------------------------------------------------

int DYTnPEffColl_t::addSystErrSource(const DYTnPEff_t &e, TString tag)
{
  DYTnPEff_t *e0= new DYTnPEff_t(e,tag);
  if (!e0) return 0;
  fTnPEffSrcV.push_back(e0);
  DYTnPEff_t *e1= new DYTnPEff_t(e,tag + "_syst");
  int relative=0; // not relative!! It should be absolute
  if (!e1 || !e1->assign_DiffAsUnc(fTnPEff,relative)) {
    std::cout << "failed to create syst unc\n";
    return 0;
  }
  fTnPEffSystV.push_back(e1);
  return 1;
}

// -------------------------------------------------------------

DYTnPEff_t* DYTnPEffColl_t::randomize(int srcIdx, TString tag) const
{
  DYTnPEff_t *e= NULL;
  int res=1;
  if (srcIdx==9999) {
    e= getTnPWithTotUnc();
  }
  else if (srcIdx==0) {
    e= new DYTnPEff_t(fElChannel);
    // randomize stat unc
    res= e->randomize(fTnPEff,tag,0);
  }
  else if (srcIdx==-1111) { // full randomization
    e= new DYTnPEff_t(fElChannel);
    res= e->randomize(fTnPEff,tag,0);
    for (unsigned int i=0; res && (i<fTnPEffSystV.size()); i++) {
      DYTnPEff_t *eRnd= new DYTnPEff_t(fElChannel);
      res= eRnd->randomize(*fTnPEffSystV[i],tag+Form("_iSyst%d",i),1);
      e->add(*eRnd);
      delete eRnd;
    }
  }
  else if (srcIdx<0) {
    // randomization of the systematic error
    srcIdx=-srcIdx-1; // -1 - correction for index
    if (srcIdx<int(fTnPEffSystV.size())) {
      e= new DYTnPEff_t(fElChannel);
      res= e->randomize(*fTnPEffSystV[srcIdx],tag,1);
    }
    else res=0;
  }
  else if (srcIdx>0) {
    // randomized for a particular source
    srcIdx--; // correction for index
    if (srcIdx<int(fTnPEffSystV.size())) {
      DYTnPEff_t *src= getTnPWithSystUnc(srcIdx);
      e= new DYTnPEff_t(fElChannel);
      res= e->randomize(*src,tag,1);
      delete src;
    }
    else res=0;
  }
  if (!res) {
    if (e) { delete e; e=NULL; }
    std::cout << "DYTnPEffColl_t::randomize failed\n";
  }
  return e;
}

// -------------------------------------------------------------

void DYTnPEffColl_t::displayAll() const
{
  for (unsigned int i=0; i<fTnPEffSrcV.size(); i++)
    fTnPEffSrcV[i]->displayAll();
  fTnPEff.displayAll();
  for (unsigned int i=0; i<fTnPEffSystV.size(); i++)
    fTnPEffSystV[i]->displayAll();
}

// -------------------------------------------------------------

void DYTnPEffColl_t::listNumbers() const
{
  std::cout << "\n" << std::string(20,'=') << "\n";
  fTnPEff.listNumbers();
  std::cout << "\n" << std::string(10,'=') << " SOURCES " << std::string(10,'=')
	    << "\n";
  for (unsigned int i=0; i<fTnPEffSrcV.size(); i++) {
    fTnPEffSrcV[i]->listNumbers();
  }
  std::cout << "\n" << std::string(10,'=') << " ABSOLUTE SYST "
	    << std::string(10,'=') << "\n";
  for (unsigned int i=0; i<fTnPEffSrcV.size(); i++) {
    fTnPEffSystV[i]->listNumbers();
  }
  std::cout << "\n" << std::string(20,'=') << "\n";
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

    // evaluate uncertainty
    double sumDRhoSqr=0.;
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
	    sumDRhoSqr += pow( (rho - sumRho/sum)/sum *
			       h2sp->GetBinError(fi1,fi2) , 2);
	  }
	}
      }
    }


    double avgRho= (sum==0.) ? 0. : sumRho/sum;
    //std::cout << "im=" << im << ", avgRho=" << avgRho << "\n";
    h1rho->SetBinContent( im+1, avgRho );
    double dAvgRho=(sumDRhoSqr<0) ? -sqrt(-sumDRhoSqr) : sqrt(sumDRhoSqr);
    h1rho->SetBinError( im+1, dAvgRho );
  }
  return h1rho;
}

// -------------------------------------------------------------

TH1D* EventSpace_t::calculateScaleFactor_misc(const DYTnPEff_t &eff,
			  int hlt4p3, TString hName, TString hTitle) const
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
	    double rho= eff.scaleFactorIdx_misc(ibin1,jbin1,ibin2,jbin2, hlt4p3);
	    //std::cout << "got rho=" << rho << "\n";
	    sumRho += rho* h2sp->GetBinContent(fi1,fi2);
	    sum += h2sp->GetBinContent(fi1,fi2);
	  }
	}
      }
    }

    // evaluate uncertainty
    double sumDRhoSqr=0.;
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
	    double rho= eff.scaleFactorIdx_misc(ibin1,jbin1,ibin2,jbin2, hlt4p3);
	    //std::cout << "got rho=" << rho << "\n";
	    sumDRhoSqr += pow( (rho - sumRho/sum)/sum *
			       h2sp->GetBinError(fi1,fi2) , 2);
	  }
	}
      }
    }

    double avgRho= (sum==0.) ? 0. : sumRho/sum;
    //std::cout << "im=" << im << ", avgRho=" << avgRho << "\n";
    h1rho->SetBinContent( im+1, avgRho );
    double dAvgRho=(sumDRhoSqr<0) ? -sqrt(-sumDRhoSqr) : sqrt(sumDRhoSqr);
    h1rho->SetBinError( im+1, dAvgRho );
  }
  return h1rho;
}

// -------------------------------------------------------------

void EventSpace_t::displayAll() const
{
  for (unsigned int im=0; im<fh2ESV.size(); im++) {
    TString mStr= DYtools::massStr(im,1);
    plotHisto(fh2ESV[im],"cES_"+ fName + "_" + mStr,
	      PlotCovCorrOpt_t(0,1,0,1.5,0.15,0.15));
  }
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
    h2->SetStats(0);
    fh2ESV.push_back(h2);
  }
  return res;
}

// -------------------------------------------------------------
// -------------------------------------------------------------
