#include "studyCSSteps.h"

// -----------------------------------------------------------
// -----------------------------------------------------------

int StudyCSSteps_t::addInfo(const CrossSection_t &cs, TString tag,
			    const TH1D* h1style)
{
  fh1YieldV.push_back(clone_local(cs.h1Yield(),tag,h1style));
  fh1SigV.push_back(clone_local(cs.h1Signal(),tag,h1style));
  fh1UnfV.push_back(clone_local(cs.h1Unf(),tag,h1style));
  fh1UnfRhoCorrV.push_back(clone_local(cs.h1UnfRhoCorr(),tag,h1style));
  fh1UnfRhoEffCorrV.push_back(clone_local(cs.h1UnfRhoEffCorr(),tag,h1style));
  fh1UnfRhoEffAccCorrV.push_back(clone_local(cs.h1UnfRhoEffAccCorr(),tag,h1style));
  fh1PostFsrV.push_back(clone_local(cs.h1UnfRhoEffAccCorr(),tag,h1style));
  fh1PreFsrV.push_back(clone_local(cs.h1PreFsr(),tag,h1style));
  fh1CSV.push_back(clone_local(cs.h1PreFsrCS(),tag,h1style));
  return 1;
}

// -----------------------------------------------------------

int StudyCSSteps_t::addInfo(const MuonCrossSection_t &muCS, TString tag,
			    const TH1D* h1style)
{
  int res=1;
  if (fTopLevel) {
    if (!fCSa) {
      fCSa= new StudyCSSteps_t();
      fCSa->fTopLevel=0;
    }
    if (!fCSb) {
      fCSb= new StudyCSSteps_t();
      fCSb->fTopLevel=0;
    }
    if (!fCSa || !fCSb) {
      std::cout << "StudyCSSteps_t::addInfo: null fCSa or fCSb\n";
      return 0;
    }
    res= (fCSa->addInfo(muCS.csA(),tag+"_csA") &&
	  fCSb->addInfo(muCS.csB(),tag+"_csB")) ? 1:0;
  }
  if (res) {
    fh1PostFsrV.push_back(clone_local(muCS.h1PostFsr(),tag,h1style));
    fh1PreFsrV.push_back(clone_local(muCS.h1PreFsr(),tag,h1style));
    fh1CSV.push_back(clone_local(muCS.h1CS(),tag,h1style));
  }
  if (!res) {
    std::cout << "error in StudyCSSteps_t::addInfo\n";
  }
  return res;
}

// -----------------------------------------------------------

TH1D* StudyCSSteps_t::clone_local(const TH1D* h1, TString tag, const TH1D *h1style)
{
  if (!h1) {
    std::cout << "creating empty histo h1" << tag << "\n";
    return (new TH1D("h1"+tag,"h1 empty histo",1,0,1));
  }
  TH1D *h1c= cloneHisto(h1,h1->GetName()+tag, h1->GetTitle()+tag);
  if (h1style) copyStyle(h1c,h1style);
  return h1c;
}

// -----------------------------------------------------------
