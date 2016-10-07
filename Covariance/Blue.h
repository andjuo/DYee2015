#ifndef BLUE_H
#define BLUE_H

#include <TROOT.h>
#include <TMatrixD.h>
#include <iostream>
#include <sstream>
#include <TGraphErrors.h>
#include <TH1D.h>

// -------------------------------------------------------------------

#ifndef Inputs_H
inline void HERE(const char *msg)
{ if (msg) std::cout << msg; std::cout << std::endl; }

template<class type_t>
inline void HERE(const char *format, type_t x)
{ std::cout << Form(format,x) << std::endl; }
#endif

// -------------------------------------------------------------------
// -------------------------------------------------------------------

void print(TString msg, const TMatrixD &m);
void print(TString msg, const TMatrixD *m);
void print(TString msg, double x);
void printCov(TString msg, const TMatrixD &m);
void printCov(TString msg, const TMatrixD *m);
int assignValues(TMatrixD &m, int ir, TString str);
int assignValuesT(TMatrixD &m, int ic, TString str);
TMatrixD covToCorr(const TMatrixD &cov);
TMatrixD covToCorrPartial(const TMatrixD &cov, const TMatrixD &covTot);
TMatrixD reduceCorrelations(const TMatrixD &cov, double by_fraction);

TMatrixD removeNaNs(const TMatrixD &m);

TGraphErrors* createGraph(TString setTitle,
			  const TMatrixD &meas, const TMatrixD &cov,
			  int nXbins, const double *xBins,
			  int shift);

TH1D *createErrHisto(TString setName, TString setTitle,
		     const TMatrixD &cov,
		     int nXbins, const double *xBins);

// -------------------------------------------------------------------

inline
TMatrixD operator*(const TMatrixD &m1, const TMatrixD *m2)
{ return (m1 * (*m2)); }

inline
TMatrixD operator*(const TMatrixD *m1, const TMatrixD &m2)
{ return ((*m1) * m2); }

inline
TMatrixD operator-(const TMatrixD &m1, const TMatrixD *m2)
{ return (m1 - (*m2)); }

inline
TMatrixD operator-(const TMatrixD *m1, const TMatrixD &m2)
{ return ((*m1) - m2); }

inline
TMatrixD operator+(const TMatrixD &m1, const TMatrixD *m2)
{ return (m1 + (*m2)); }

inline
TMatrixD operator+(const TMatrixD *m1, const TMatrixD &m2)
{ return ((*m1) + m2); }

inline
TMatrixD mult(const TMatrixD *m1, const TMatrixD *m2)
{ return ((*m1) * (*m2)); }

// -------------------------------------------------------------------
// -------------------------------------------------------------------

struct Info_t {
  TString description;
  std::vector<TString> info;
public:
  Info_t(TString set_description="study") : description(set_description),info()  {}

  Info_t(TString set_description, const std::vector<TString> &set_info) :
    description(set_description), info(set_info) {}

  unsigned int size() const { return info.size(); }
  void add(TString line) { info.push_back(line); }

  template<class idx_t>
  TString at(idx_t i) const { return info[i]; }

  void print() const;
};

// -------------------------------------------------------------------
// -------------------------------------------------------------------

struct BLUEResult_t {
  TMatrixD *meas; // total measurement
  TMatrixD *U; // relational matrix - measurement--its-result
  TMatrixD *covInp; // input covariance matrix
  TMatrixD *lambda;  // weights
  TMatrixD *est; // best estimate
  TMatrixD *covOut; // output covariance matrix
  TMatrixD *corrOut; // output correlation matrix
  double chi2;
public:
  BLUEResult_t();
  BLUEResult_t(const BLUEResult_t &r);
  BLUEResult_t(const TMatrixD &set_meas, const TMatrixD &set_relations,
	       const TMatrixD &set_covariance);
  ~BLUEResult_t();
  const TMatrixD* getMeas() const { return meas; }
  double getMeas(int ir, int ic=0) const { return (*meas)(ir,ic); }
  const TMatrixD* getU() const { return U; }
  double getU(int ir, int ic) const { return (*U)(ir,ic); }
  const TMatrixD* getInpCov() const { return covInp; }
  double getInpCov(int ir, int ic) const { return (*covInp)(ir,ic); }
  double getInpErr(int i) const;
  const TMatrixD* getLambda() const { return lambda; }
  double getLambda(int ir, int ic) const { return (*lambda)(ir,ic); }
  const TMatrixD* getEstimate() const { return est; }
  const TMatrixD* getEst() const { return est; }
  double getEst(int ir, int ic=0) const { return (*est)(ir,ic); }
  const TMatrixD* getCov() const { return covOut; }
  double getCov(int ir, int ic) const { return (*covOut)(ir,ic); }
  double getErr(int i) const;
  const TMatrixD* getCorr() const { return corrOut; }

  double getChi2() const { return chi2; }
  double getChi2(int iMeas) const;

  BLUEResult_t& operator=(const BLUEResult_t &r);
  int estimate();
  int estimate(const TMatrixD &set_meas, const TMatrixD &set_relations,
	       const TMatrixD &set_cov);
  int estimate(const TMatrixD &set_meas, const TMatrixD &set_cov);
  int estimate(const TMatrixD &set_measA, const TMatrixD &set_covA,
	       const TMatrixD &set_measB, const TMatrixD &set_covB,
	       const TMatrixD &set_corrAB);


  // arrange according to Valassi paper 2003
  int estimateTest(const TMatrixD &set_measA, const TMatrixD &set_covA,
		   const TMatrixD &set_measB, const TMatrixD &set_covB,
		   const TMatrixD &set_corrAB);

  // more studies
  TMatrixD combinedLambda(int iMeas1, int iMeas2, double factor) const;
  TMatrixD combinedEst(int iMeas1, int iMeas2, double factor) const;
  TMatrixD combinedCov(int iMeas1, int iMeas2, double factor) const;
  //double getCombinedChi2(int iMeas1, int iMeas2, double factor) const;

  // covariance by source, eq.(18)
  static TMatrixD measCovMatrix(const TMatrixD &inpCov, int isMeasA);
  TMatrixD contributedCov(const TMatrixD &covSrc) const
  { return (*lambda) * covSrc * TMatrixD(TMatrixD::kTransposed,*lambda); }

  // checks
  static int compareDimRC(const TMatrixD &a, const TMatrixD &b)
  { return ((a.GetNrows()==b.GetNrows()) && (a.GetNcols()==b.GetNcols()))?1:0; }
  static int compareDimR(const TMatrixD &a, const TMatrixD &b)
  { return (a.GetNrows()==b.GetNrows()) ? 1:0; }
  static int compareDimC(const TMatrixD &a, const TMatrixD &b)
  { return (a.GetNcols()==b.GetNcols()) ? 1:0; }
  static int isSquare(const TMatrixD &a)
  { return (a.GetNrows()==a.GetNcols()) ? 1:0; }

};

// -------------------------------------------------------------------

#endif
