#ifndef studyCSSteps_H
#define studyCSSteps_H

#include "crossSection.h"

// -----------------------------------------------------------

struct StudyCSSteps_t
{
  std::vector<TH1D*> fh1YieldV,fh1SigV,fh1UnfV,fh1UnfRhoCorrV,fh1UnfRhoEffCorrV;
  std::vector<TH1D*> fh1UnfRhoEffAccCorrV; // postFSR
  std::vector<TH1D*> fh1PostFsrV, fh1PreFsrV, fh1CSV; // top level
  StudyCSSteps_t *fCSa,*fCSb;
  int fTopLevel;
public:
StudyCSSteps_t() : fh1YieldV(),fh1SigV(),fh1UnfV(),fh1UnfRhoCorrV(),
    fh1UnfRhoEffCorrV(), fh1UnfRhoEffAccCorrV(),
    fh1PostFsrV(), fh1PreFsrV(), fh1CSV(),
    fCSa(NULL), fCSb(NULL),
    fTopLevel(1)
  {}

  std::vector<TH1D*>& h1YieldV() { return fh1YieldV; }
  std::vector<TH1D*>& h1SigV() { return fh1SigV; }
  std::vector<TH1D*>& h1UnfV() { return fh1UnfV; }
  std::vector<TH1D*>& h1UnfRhoCorrV() { return fh1UnfRhoCorrV; }
  std::vector<TH1D*>& h1UnfRhoEffCorrV() { return fh1UnfRhoEffCorrV; }
  std::vector<TH1D*>& h1UnfRhoEffAccCorrV() { return fh1UnfRhoEffAccCorrV; }
  std::vector<TH1D*>& h1PostFsrV() { return fh1PostFsrV; }
  std::vector<TH1D*>& h1PreFsrV() { return fh1PreFsrV; }
  std::vector<TH1D*>& h1CSV() { return fh1CSV; }
  StudyCSSteps_t* csA() { return fCSa; }
  StudyCSSteps_t* csB() { return fCSb; }

  int addInfo(const CrossSection_t &cs, TString tag, const TH1D *h1style=NULL);
  int addInfo(const MuonCrossSection_t &muCS, TString tag, const TH1D *h1style=NULL);
  TH1D* clone_local(const TH1D* h1, TString tag, const TH1D *h1style=NULL);
};

// -----------------------------------------------------------


#endif
