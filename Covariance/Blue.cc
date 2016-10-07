#include "Blue.h"

// -------------------------------------------------------------------
// -------------------------------------------------------------------

void print(TString msg, const TMatrixD &m)
{
  std::cout << msg << Form(" [%d][%d]",m.GetNrows(),m.GetNcols()) << "\n";
  for (int ir=0; ir<m.GetNrows(); ir++) {
    std::cout << "ir=" << ir << " : ";
    for (int ic=0; ic<m.GetNcols(); ic++) {
      std::cout << " " << Form("%6.3lf",m(ir,ic));
    }
    std::cout << "\n";
  }
}

// -------------------------------------------------------------------

void print(TString msg, const TMatrixD *m)
{
  if (m==NULL) {
    std::cout << "error: print(msg=" << msg << ",NULL TMatrixD) is called\n";
    return;
  }
  print(msg,*m);
}

// -------------------------------------------------------------------

void print(TString msg, double x)
{
  std::cout << msg << "  value=" << x << "\n";
}

// -------------------------------------------------------------------

void printCov(TString msg, const TMatrixD &m)
{
  TMatrixD corr( covToCorr(m) );
  std::cout << msg << Form(" [%d][%d]",m.GetNrows(),m.GetNcols()) << "\n";
  for (int ir=0; ir<m.GetNrows(); ir++) {
    std::cout << "ir=" << ir << " : ";
    for (int ic=0; ic<m.GetNcols(); ic++) {
      std::cout << " " << Form("%6.3lf",m(ir,ic));
    }
    std::cout << "   ";
    for (int ic=0; ic<corr.GetNcols(); ic++) {
      std::cout << " " << Form("%6.3lf",corr(ir,ic));
    }
    std::cout << "\n";
  }
}

// -------------------------------------------------------------------

void printCov(TString msg, const TMatrixD *m)
{
  if (m==NULL) {
    std::cout << "error: printCov(msg=" << msg << ",NULL TMatrixD) is called\n";
    return;
  }
  printCov(msg,*m);
}

// -------------------------------------------------------------------

int assignValues(TMatrixD &m, int ir, TString str)
{
  std::stringstream ss(str.Data());
  int ic=0;
  double x;
  while (!ss.eof() && (ic<m.GetNcols())) {
    ss >> x;
    m(ir,ic) = x;
    ic++;
  }
  if (ic!=m.GetNcols()) {
    std::cout << "assignValues ir=" << ir << ", str=<" << str
	      << "> did not contain enough numbers\n";
    return 0;
  }
  return 1;
}

// -------------------------------------------------------------------

int assignValuesT(TMatrixD &m, int ic, TString str)
{
  std::stringstream ss(str.Data());
  int ir=0;
  double x;
  while (!ss.eof() && (ir<m.GetNrows())) {
    ss >> x;
    m(ir,ic) = x;
    ir++;
  }
  if (ir!=m.GetNrows()) {
    std::cout << "assignValuesT ic=" << ic << ", str=<" << str
	      << "> did not contain enough numbers\n";
    return 0;
  }
  return 1;
}

// -------------------------------------------------------------------

TMatrixD covToCorr(const TMatrixD &m)
{
  TMatrixD corr(TMatrixD::kZero,m);
  for (int ir=0; ir<m.GetNrows(); ir++) {
    for (int ic=0; ic<m.GetNcols(); ic++) {
      corr(ir,ic) = m(ir,ic)/ sqrt(m(ir,ir)*m(ic,ic));
    }
  }
  //print("matrix for corr",m);  print("corr",corr);
  return corr;
}

// -------------------------------------------------------------------

TMatrixD covToCorrPartial(const TMatrixD &m, const TMatrixD &covTot)
{
  TMatrixD corr(TMatrixD::kZero,m);
  for (int ir=0; ir<m.GetNrows(); ir++) {
    for (int ic=0; ic<m.GetNcols(); ic++) {
      corr(ir,ic) = m(ir,ic)/ sqrt(covTot(ir,ir)*covTot(ic,ic));
    }
  }
  //print("matrix for corr partial",m);  print("corr",corr);
  return corr;
}

// -------------------------------------------------------------------

TMatrixD reduceCorrelations(const TMatrixD &cov, double by_fraction)
{
  TMatrixD m(cov);
  for (int ir=0; ir<m.GetNrows(); ir++) {
    for (int ic=0; ic<m.GetNcols(); ic++) {
      double f= (ir==ic) ? 1 : (1.-by_fraction);
      m(ir,ic) = f * cov(ir,ic);
    }
  }
  return m;
}

// -------------------------------------------------------------------

TMatrixD removeNaNs(const TMatrixD &m)
{
  TMatrixD mNew(m);
  for (int ir=0; ir<m.GetNrows(); ir++)
    for (int ic=0; ic<m.GetNcols(); ic++) {
      if (m(ir,ic)!=m(ir,ic)) mNew(ir,ic)=0.;
    }
  return mNew;
}

// -------------------------------------------------------------------

TGraphErrors* createGraph(TString setTitle,
			  const TMatrixD &meas, const TMatrixD &cov,
			  int nXbins, const double *xBins,
			  int shift)
{
  if (meas.GetNrows()!=nXbins) {
    std::cout << "measurement does not match the number of bins\n";
    return NULL;
  }
  double *xc= new double[nXbins];
  double *yc= new double[nXbins];
  double *ex= new double[nXbins];
  double *ey= new double[nXbins];
  for (int i=0; i<nXbins; i++) {
    double dx=0;
    if (xBins[nXbins]>1000.) {
      if (xBins[i]<=60) dx=1.;
      else if (xBins[i]<=76) dx=0.8;
      else if (xBins[i]<=141) dx=1.25;
      else if (xBins[i]<=200) dx=2;
      else if (xBins[i]<=250) dx=5;
      else if (xBins[i]<=380) dx=8;
      else if (xBins[i]<=600) dx=10;
      else if (xBins[i]<=1000) dx=20;
      else dx=50;
    }
    else {
      dx=0.25;
    }
    dx*=shift;
    xc[i]= 0.5*(xBins[i+1]+xBins[i]) + dx;
    yc[i]= meas(i,0);
    ex[i]= (shift==0) ?  0.5*(xBins[i+1]-xBins[i]) : 0.;
    ey[i]= sqrt(cov(i,i));
    if (0) {
      std::cout << "i=" << i << " x range " << xBins[i] << " .. "
		<< xBins[i+1] << "\n";
      std::cout << "i=" << i << " x=" << xc[i] << " +- " << ex[i]
		<< ", y=" << yc[i] << " +- " << ey[i] << "\n";
    }
  }
  TGraphErrors *gr= new TGraphErrors(nXbins,xc,yc,ex,ey);
  gr->SetTitle(setTitle);
  delete xc; delete yc; delete ex; delete ey;
  return gr;
}

// -------------------------------------------------------------------

TH1D *createErrHisto(TString setName, TString setTitle,
		     const TMatrixD &cov,
		     int nXbins, const double *xBins)
{
  TH1D *h1= new TH1D(setName,setTitle,nXbins,xBins);
  h1->SetDirectory(0);
  for (int i=0; i<nXbins; i++) {
    double err=(cov(i,i)>=0) ? sqrt(cov(i,i)) : -sqrt(-cov(i,i));
    h1->SetBinContent(i+1, err);
    h1->SetBinError  (i+1, 0.);
  }
  return h1;
}

// -------------------------------------------------------------------
// -------------------------------------------------------------------

void Info_t::print() const
{
  std::cout << "Description: " << description << "\n";
  std::cout << "\n";
  for (unsigned int i=0; i<info.size(); i++) {
    std::cout << info[i] << "\n";
  }
  std::cout << "\n";
}

// -------------------------------------------------------------------
// -------------------------------------------------------------------

BLUEResult_t::BLUEResult_t() :
  meas(NULL), U(NULL), covInp(NULL), lambda(NULL), est(NULL), covOut(NULL),
  corrOut(NULL),
  chi2(0)
{
}

// -------------------------------------------------------------------

BLUEResult_t::BLUEResult_t(const BLUEResult_t &r) :
  meas(NULL), U(NULL), covInp(NULL), lambda(NULL), est(NULL), covOut(NULL),
  corrOut(NULL),
  chi2(0)
{
  this->operator=(r);
}

// -------------------------------------------------------------------

BLUEResult_t::BLUEResult_t(const TMatrixD &set_meas,
			  const TMatrixD &set_relations,
			  const TMatrixD &set_covariance) :
  meas(new TMatrixD(set_meas)),
  U(new TMatrixD(set_relations)),
  covInp(new TMatrixD(set_covariance)),
  lambda(NULL),
  est(NULL),
  covOut(NULL),
  corrOut(NULL),
  chi2(0)
{
}


// -------------------------------------------------------------------

BLUEResult_t::~BLUEResult_t()
{
  if (meas) delete meas;
  if (U) delete U;
  if (covInp) delete covInp;
  if (lambda) delete lambda;
  if (est) delete est;
  if (covOut) delete covOut;
  if (corrOut) delete corrOut;
}

// -------------------------------------------------------------------

double BLUEResult_t::getInpErr(int i) const
{
  double c= (*covInp)(i,i);
  return (c>=0) ? sqrt(c) : -sqrt(c);
}

// -------------------------------------------------------------------

double BLUEResult_t::getErr(int i) const
{
  double c= (*covOut)(i,i);
  return (c>=0) ? sqrt(c) : -sqrt(c);
}

// -------------------------------------------------------------------

double BLUEResult_t::getChi2(int iMeas) const
{
  if (!est) return -999.;
  TMatrixD sel(TMatrixD::kZero,*covInp);
  for (int ir=0; ir<sel.GetNrows(); ir++) {
    if ((*U)(ir,iMeas)==1) sel(ir,ir)=1;
  }
  //print("sel",sel);
  TMatrixD covInpInv(TMatrixD::kInverted,*covInp);
  TMatrixD meas_diff_est( meas - mult(U,est) );
  TMatrixD meas_diff_est_tr(TMatrixD::kTransposed, meas_diff_est);
  TMatrixD chiSel( meas_diff_est_tr * covInpInv * sel * meas_diff_est );
  return chiSel(0,0);
}

// -------------------------------------------------------------------

BLUEResult_t& BLUEResult_t::operator=(const BLUEResult_t &r)
{
  if (this!=&r) {
    if (meas) delete meas;
    meas= new TMatrixD(*r.meas);
    if (U) delete U;
    U= new TMatrixD(*r.U);
    if (covInp) delete covInp;
    covInp= new TMatrixD(*r.covInp);
    if (lambda) delete lambda;
    lambda= new TMatrixD(*r.lambda);
    if (est) delete est;
    est= new TMatrixD(*r.est);
    if (covOut) delete covOut;
    covOut= new TMatrixD(*r.covOut);
    if (corrOut) delete corrOut;
    corrOut= new TMatrixD(*r.corrOut);
    chi2= r.chi2;
  }
  return *this;
}

// -------------------------------------------------------------------

int BLUEResult_t::estimate()
{
  if (!meas || !U || !covInp) {
    std::cout << "BLUEResult_t::estimate: input is not assigned\n";
    return 0;
  }
  if (lambda) delete lambda;
  if (est) delete est;
  if (covOut) delete covOut;
  if (corrOut) delete corrOut;
  chi2=0;

  std::cout << "\nestimate()\n";
  print("meas",meas);
  print("covInp",covInp);

  TMatrixD Ut(TMatrixD::kTransposed,*U);

  print("Ut*meas",Ut*meas);

  TMatrixD covInpInv(*covInp);
  double det=0;
  covInpInv.Invert(&det);
  print(Form("covTotInv [det=%g]",det),covInpInv);

  TMatrixD nom = Ut * covInpInv; // eq.8
  //print("nom=Ut*covInpInv",Ut*covInpInv);

  TMatrixD denomM= Ut * (covInpInv * U);
  //print("denomM",denomM);
  denomM.Invert(&det);
  print(Form("1/denomM [det=%g]",det),denomM);

  lambda= new TMatrixD(denomM * nom);
  print("lambda",lambda);
  lambda->Print();

  print("lambda*U should be ones", mult(lambda,U));

  est= new TMatrixD( mult(lambda,meas) );
  print("estimate",est);
  print("U*est",mult(U,est));

  TMatrixD lambdaTr(TMatrixD::kTransposed,*lambda);
  //print("lambdaTr",lambdaTr);
  //print("lambdaTr*covInp*lambda",lambdaTr * mult(covInp,lambda));
  covOut= new TMatrixD( mult(lambda,covInp) * lambdaTr );
  print("finCov",covOut);

  print("denomM_minus_covOut", denomM - covOut );

  corrOut= new TMatrixD(covToCorr(*covOut));
  print("finCorr",corrOut);

  // calculate chi2
  TMatrixD meas_diff_est ( meas - mult(U,est) );
  //print("meas_diff_est",meas_diff_est);

  TMatrixD meas_diff_est_tr(TMatrixD::kTransposed,meas_diff_est);

  TMatrixD chi2_asM= meas_diff_est_tr * covInpInv * meas_diff_est;
  print("chi2",chi2_asM);
  chi2= chi2_asM(0,0);
  return 1;
}

// -------------------------------------------------------------------

int BLUEResult_t::estimate(const TMatrixD &set_meas, const TMatrixD &set_U,
			    const TMatrixD &set_cov)
{
  if (meas) delete meas;
  meas= new TMatrixD(set_meas);
  if (U) delete U;
  U= new TMatrixD(set_U);
  if (covInp) delete covInp;
  covInp= new TMatrixD(set_cov);
  return estimate();
}

// -------------------------------------------------------------------

int BLUEResult_t::estimate(const TMatrixD &set_meas, const TMatrixD &set_cov)
{
  if (set_meas.GetNcols()!=1) {
    std::cout << "estimate(meas,cov) expects meas(N,1)\n";
    return 0;
  }
  if (meas) delete meas;
  meas= new TMatrixD(set_meas);
  if (U) delete U;
  U= new TMatrixD(set_meas.GetNrows(),1);
  for (int ir=0; ir<U->GetNrows(); ir++) (*U)(ir,0)=1;
  if (covInp) delete covInp;
  covInp= new TMatrixD(set_cov);
  return estimate();
}

// -------------------------------------------------------------------

int BLUEResult_t::estimate(const TMatrixD &set_measA, const TMatrixD &set_covA,
			   const TMatrixD &set_measB, const TMatrixD &set_covB,
			   const TMatrixD &set_covAB)
{
  if (!compareDimRC(set_measA,set_measB)) {
    std::cout << "estimate A+B is not ready for different number of measurements\n";
    return 0;
  }

  if (meas) delete meas;
  if (U) delete U;
  if (covInp) delete covInp;

  int dim=set_measA.GetNrows();

  U= new TMatrixD(2*dim,dim);
  for (int ir=0; ir<U->GetNrows(); ir++) {
    int ic= (ir < dim) ? ir : (ir-dim);
    (*U)(ir,ic) = 1;
  }
  print("U",U);

  if (0) {
    meas= new TMatrixD(2*dim,dim);
    for (int ir=0; ir<dim; ir++) {
      (*meas)(ir    ,ir) = set_measA(ir,0);
      (*meas)(ir+dim,ir) = set_measB(ir,0);
    }
  }
  else {
    meas= new TMatrixD(2*dim,1);
    for (int ir=0; ir<dim; ir++) {
      (*meas)(ir    ,0) = set_measA(ir,0);
      (*meas)(ir+dim,0) = set_measB(ir,0);
    }
  }

  covInp= new TMatrixD(2*dim,2*dim);
  if ((set_covA.GetNrows()==0) && (set_covB.GetNrows()==0)) {
    if (!compareDimRC(*covInp,set_covAB)) {
      std::cout << "compare(A,B): set_covAB is wrong size\n";
      return 0;
    }
    (*covInp)=set_covAB;
  }
  else {
    for (int ir=0; ir<set_covA.GetNrows(); ir++) {
      for (int ic=0; ic<set_covA.GetNcols(); ic++) {
	(*covInp)(ir,ic)= set_covA(ir,ic);
      }
    }
    for (int ir=0; ir<set_covB.GetNrows(); ir++) {
      for (int ic=0; ic<set_covB.GetNcols(); ic++) {
	(*covInp)(ir+dim,ic+dim) = set_covB(ir,ic);
      }
    }
    for (int ir=0; ir<set_covAB.GetNrows(); ir++) {
      for (int ic=0; ic<set_covAB.GetNcols(); ic++) {
	(*covInp)(ir    ,ic+dim) = set_covAB(ir,ic);
	(*covInp)(ic+dim,ir)= set_covAB(ir,ic);
      }
    }
  }
  print("covInp",covInp);
  return estimate();
}

// -------------------------------------------------------------------

int BLUEResult_t::estimateTest(
		     const TMatrixD &set_measA, const TMatrixD &set_covA,
		     const TMatrixD &set_measB, const TMatrixD &set_covB,
		     const TMatrixD &set_covAB)
{
  if (!compareDimRC(set_measA,set_measB)) {
    std::cout << "estimate A+B is not ready for different number of measurements\n";
    return 0;
  }

  if (meas) delete meas;
  if (U) delete U;
  if (covInp) delete covInp;

  int dim=set_measA.GetNrows();

  U= new TMatrixD(2*dim,2);
  for (int ir=0; ir<U->GetNrows(); ir++) {
    int ic= (ir < dim) ? 0 : 1;
    (*U)(ir,ic) = 1;
  }
  print("U",U);

  meas= new TMatrixD(2*dim,1);
  for (int ir=0; ir<dim; ir++) {
    (*meas)(2*ir  ,0) = set_measA(ir,0);
    (*meas)(2*ir+1,0) = set_measB(ir,0);
  }

  covInp= new TMatrixD(2*dim,2*dim);
  if ((set_covA.GetNrows()==0) && (set_covB.GetNrows()==0)) {
    if (!compareDimRC(*covInp,set_covAB)) {
      std::cout << "compare(A,B): set_covAB is wrong size\n";
      return 0;
    }
    (*covInp)=set_covAB;
  }
  else {
    for (int ir=0; ir<set_covA.GetNrows(); ir++) {
      for (int ic=0; ic<set_covA.GetNcols(); ic++) {
	(*covInp)(2*ir,2*ic)= set_covA(ir,ic);
      }
    }
    for (int ir=0; ir<set_covB.GetNrows(); ir++) {
      for (int ic=0; ic<set_covB.GetNcols(); ic++) {
	(*covInp)(2*ir+1,2*ic+1) = set_covB(ir,ic);
      }
    }
    for (int ir=0; ir<set_covAB.GetNrows(); ir++) {
      for (int ic=0; ic<set_covAB.GetNcols(); ic++) {
	(*covInp)(2*ir,2*ic+1) = set_covAB(ir,ic);
	(*covInp)(2*ic+1, 2*ir)= set_covAB(ir,ic);
      }
    }
  }
  print("covInp",covInp);
  return estimate();
}

// -------------------------------------------------------------------

TMatrixD BLUEResult_t::combinedLambda(int iMeas1, int iMeas2, double factor) const
{
  TMatrixD combL(1,lambda->GetNcols());
  for (int ic=0; ic<combL.GetNcols(); ic++) {
    combL(0,ic) = (*lambda)(iMeas1,ic) + factor * (*lambda)(iMeas2,ic);
  }
  return combL;
}

// -------------------------------------------------------------------

TMatrixD BLUEResult_t::combinedEst(int iMeas1, int iMeas2, double factor) const
{
  TMatrixD diff= combinedLambda(iMeas1,iMeas2,factor);
  //print("combinedLambda",diff);
  TMatrixD combEst (diff * meas);
  //print("combEst",combEst);
  return combEst;
}

// -------------------------------------------------------------------

TMatrixD BLUEResult_t::combinedCov(int iMeas1, int iMeas2, double factor) const
{
  TMatrixD diff=  combinedLambda(iMeas1,iMeas2,factor);
  TMatrixD diff_tr (TMatrixD::kTransposed, diff);
  TMatrixD combCov( diff * covInp * diff_tr );
  //print("combCov",combCov);
  return combCov;
}

// -------------------------------------------------------------------

TMatrixD BLUEResult_t::measCovMatrix(const TMatrixD &inpCov, int isMeasA)
{
  int dim=inpCov.GetNrows();
  TMatrixD cov(2*dim,2*dim);
  int irShift=(isMeasA!=0) ? 0 : dim; // dim, 0, 0
  int icShift=(isMeasA==1) ? 0 : dim; // dim, 0, dim
  for (int ir=0; ir<inpCov.GetNrows(); ir++) {
    for (int ic=0; ic<inpCov.GetNcols(); ic++) {
      cov(ir+irShift,ic+icShift) = inpCov(ir,ic);
      if (isMeasA==2) {
	cov(ic+icShift,ir+irShift) = inpCov(ir,ic);
      }
    }
  }
  return cov;
}

// -------------------------------------------------------------------

/*
double BLUEResult_t::getCombinedChi2(int iMeas1, int iMeas2, double factor) const
{
  TMatrixD combEst= combinedEst(iMeas1,iMeas2,factor);
  if (0) combEst.Print();
  return 0.;
}
*/

// -------------------------------------------------------------------
// -------------------------------------------------------------------
