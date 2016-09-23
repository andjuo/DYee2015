// Examples from Valassi 2003

#include <TMatrixD.h>
#include <iostream>
#include <math.h>
#include <sstream>
#include <TString.h>

// ---------------------------------------------------------------

// ---------------------------------------------------------------

inline void HERE(const char *msg)
{
  if (msg) std::cout << msg;
  std::cout << std::endl;
}

// ------------------------------------------

struct Info_t {
  TString caseName;
  std::vector<TString> info;
public:
  Info_t(TString set_caseName="study") : caseName(set_caseName), info() {}
  Info_t(TString set_caseName, const std::vector<TString> &set_info) :
    caseName(set_caseName), info(set_info) {}
  unsigned int size() const { return info.size(); }
  void add(TString line) { info.push_back(line); }
  template<class idx_t>
  TString at(idx_t i) const { return info[i]; }
};

// ------------------------------------------

// ------------------------------------------

void print(TString msg, const TMatrixD &m);
//void print(TString msg, const TMatrixD &m, const TMatrixD &cov);
void print(TString msg, double x);
int assignValues(TMatrixD &m, int ir, TString str);
int assignValuesT(TMatrixD &m, int ic, TString str);

// combine uncorrelated measurements
void combine(const TMatrixD &measA_inp, const TMatrixD &measB_inp,
	     const TMatrixD &covTot_inp, const Info_t *info=NULL);
		   
// ---------------------------------------------------------------
// ---------------------------------------------------------------
// ---------------------------------------------------------------

void myBlue(int theCase=1, int printInfo=0)
{

  TMatrixD *measAinp=NULL, *measBinp=NULL;
  TMatrixD *covInp=NULL;
  int nMeas=2;
  TString caseName="unknown";
  Info_t info;
  info.add(Form("Case number in the macro theCase=%d",theCase));

  if (theCase==1) {
    caseName="4.1 no correlations";
    nMeas=2;
    measAinp=new TMatrixD(2,1);
    measBinp=new TMatrixD(2,1);
    assignValuesT(*measAinp,0,"10.50 9.50");
    assignValuesT(*measBinp,0,"13.50 14.00");
    covInp=new TMatrixD(4,4);
    assignValues(*covInp,0,"1.00 0    0    0   ");
    assignValues(*covInp,1,"0    9.00 0    0   ");
    assignValues(*covInp,2,"0    0    9.00 0   ");
    assignValues(*covInp,3,"0    0    0    9.00");
    info.add("Eq.33 defines the measurements.");
    info.add("Eq.34(&36)-35 describe totCov and U matrix.");
    info.add("Results: eqs.37,38\n");
  }
  else if (theCase==2) {
    caseName="4.1 no correlations. Lepton universality";
    nMeas=1;
    measAinp=new TMatrixD(4,1);
    measBinp=new TMatrixD(0,0);
    assignValuesT(*measAinp,0,"10.50 13.50   9.50 14.00");
    covInp=new TMatrixD(4,4);
    assignValues(*covInp,0,"1.00 0    0    0   ");
    assignValues(*covInp,1,"0    9.00 0    0   ");
    assignValues(*covInp,2,"0    0    9.00 0   ");
    assignValues(*covInp,3,"0    0    0    9.00");
    info.add("Results: eqs.39 and 40");
  }
  else if (theCase==31) {
    caseName="4.2 Correlations between B_e";
    nMeas=2;
    measAinp=new TMatrixD(2,1);
    measBinp=new TMatrixD(2,1);
    assignValuesT(*measAinp,0,"10.50 9.50");
    assignValuesT(*measBinp,0,"13.50 14.00");
    covInp=new TMatrixD(4,4);
    assignValues(*covInp,0,"1.00 0.45 0    0   ");
    assignValues(*covInp,1,"0.45 9.00 0    0   ");
    assignValues(*covInp,2,"0    0    9.00 0   ");
    assignValues(*covInp,3,"0    0    0    9.00");
    info.add("Tot cov from eq.41, results eq.42");
  }
  else if ((theCase==3) || (theCase==4) || (theCase==5) || (theCase==6)) {
    double corr=(theCase==4) ? 1/3. : 1/2.;
    if (theCase==3) corr=0.15;
    else if (theCase==6) corr=-0.15;
    caseName=TString("4.2 Correlations between B_e measurements") +
      Form(" (correlation %4.1lf%%)",corr*100.);
    nMeas=2;
    measAinp=new TMatrixD(2,1);
    measBinp=new TMatrixD(2,1);
    assignValuesT(*measAinp,0,"10.50 9.50");
    assignValuesT(*measBinp,0,"13.50 14.00");
    covInp=new TMatrixD(4,4);
    assignValues(*covInp,0,Form("1.00 %lf 0    0   ",corr*sqrt(1*9)));
    assignValues(*covInp,1,Form("%lf 9.00 0    0   ",corr*sqrt(1*9)));
    assignValues(*covInp,2,"0    0    9.00 0   ");
    assignValues(*covInp,3,"0    0    0    9.00");
    if (theCase==3) info.add("Tot cov from eq.41, results eq.42");
    else if (theCase==4) info.add("Tot cov similar to eq.41, results discussed below eq.42");
    else if (theCase==5) info.add("Tot cov similar to eq.41, results probe the decrease of uncertainty for correlations larger than 1/3");
    else if (theCase==6) info.add("Tot cov similar to eq.41, results discussed on p.401 top left");
  }
  else if ((theCase==7) || (theCase==8)) {
    double corr= 0.15;
    if (theCase==8) corr=-0.15;
    caseName=TString("4.2 Correlations between B_e measurements") +
      Form(" (correlation %4.1lf%%)",corr*100.);
    caseName+= TString(". Lepton universality");
    nMeas=1;
    measAinp=new TMatrixD(4,1);
    measBinp=new TMatrixD(0,0);
    assignValuesT(*measAinp,0,"10.50 13.50  9.50 14.00");
    covInp=new TMatrixD(4,4);
    assignValues(*covInp,0,Form("1.00 %lf 0    0",corr*sqrt(1*9)));
    assignValues(*covInp,1,Form("%lf 9.00 0    0",corr*sqrt(1*9)));
    assignValues(*covInp,2,"0    0    9.00 0   ");
    assignValues(*covInp,3,"0    0    0    9.00");
    info.add("Tot cov similar to eq.41. Results discussed in the last paragraph of 4.2");
  }
  else if ((theCase==9) || (theCase==11)) {
    double corr=(theCase==9) ? 0.995 : -0.995;
    caseName=TString("4.3 Correlations between B_e and B_tau measurementB") +
      Form(" (correlation %4.1lf%%)",corr*100.);
    nMeas=2;
    measAinp=new TMatrixD(2,1);
    measBinp=new TMatrixD(2,1);
    assignValuesT(*measAinp,0,"10.50 9.50");
    assignValuesT(*measBinp,0,"13.50 14.00");
    covInp=new TMatrixD(4,4);
    assignValues(*covInp,0,"1.00 0    0    0");
    assignValues(*covInp,1,Form("0 9.00 0    %lf",corr*sqrt(9*9)));
    assignValues(*covInp,2,"0    0    9.00 0");
    assignValues(*covInp,3,Form("0 %lf 0    9.00",corr*sqrt(9*9)));
    if (theCase==9) info.add("Tot cov eq.44. Results discussed eqs.45-46");
    else if (theCase==11) info.add("Tot cov similar to eq.44. Results discussed in Eq.49");
  }
  else if ((theCase==10) || (theCase==12)) {
    double corr=(theCase==10) ? 0.995 : -0.995;
    caseName=TString("4.3 Correlations between B_e and B_tau measurementB") +
      Form(" (correlation %4.1lf%%)",corr*100.);
    caseName+=". Lepton universality";
    nMeas=1;
    measAinp=new TMatrixD(4,1);
    measBinp=new TMatrixD(0,00);
    assignValuesT(*measAinp,0,"10.50 13.50  9.50 14.00");
    covInp=new TMatrixD(4,4);
    assignValues(*covInp,0,"1.00 0    0    0");
    assignValues(*covInp,1,Form("0 9.00 0    %lf",corr*sqrt(9*9)));
    assignValues(*covInp,2,"0    0    9.00 0");
    assignValues(*covInp,3,Form("0 %lf 0    9.00",corr*sqrt(9*9)));
    if (theCase==10) info.add("Tot cov eq.44. Results discussed eq.48");
    else info.add("Tot cov. similar to eq.44. Results discussed eq.50");
  }
  else {
    std::cout << "not ready for case=" << theCase << "\n";
    return;
  }
  info.caseName=caseName;

  TMatrixD measA(*measAinp);
  TMatrixD measB(*measBinp);
  TMatrixD covTot(*covInp);
  print("measA",measA);
  print("measB",measB);

  combine(measA,measB,covTot,(printInfo) ? &info : NULL);
}

// -------------------------------------------------------------
// -------------------------------------------------------------

// combine uncorrelated measurements
void combine(const TMatrixD &measA_inp, const TMatrixD &measB_inp,
	     const TMatrixD &covTot_inp, const Info_t *info)
{
  TMatrixD measA(measA_inp);
  TMatrixD measB(measB_inp);
  TMatrixD covTot(covTot_inp);

  int nMeas= (measB.GetNrows()==0) ? 1 : 2;

  print("measA",measA);
  print("measB",measB);


  int dim1= measA.GetNrows();
  int dim2= measB.GetNrows();
  int dimTot=dim1+dim2;

  TMatrixD U(dimTot,nMeas);
  if (nMeas==1) {
    for (int ir=0; ir<U.GetNrows(); ir++) U(ir,0)=1;
  }
  else {
    for (int ir=0; ir<U.GetNrows(); ir++) {
      int i=(ir < U.GetNrows()/2) ? 0 : 1;
      U(ir,i)=1;
    }
  }
  print("U",U);
  TMatrixD Ut(TMatrixD::kTransposed,U);
  //Ut.T(); // transpose
  print("Ut",Ut);

  TMatrixD measTot(dimTot,1); //rows,columns
  if (nMeas==1) {
    measTot=measA;
  }
  else {
    for (int ir=0; ir<dim1; ir++) {
      measTot(2*ir  ,0)= measA(ir,0);
      measTot(2*ir+1,0)= measB(ir,0);
    }
  }

  print("covTot",covTot);
  print("measTot",measTot);

  TMatrixD covTotU= covTot*U;
  print("covTotU",covTotU);

  TMatrixD covTotInv= covTot;
  double det=0;
  covTotInv.Invert(&det);
  print(Form("covTotInv [det=%g]",det),covTotInv);

  TMatrixD nom= Ut*covTotInv; // eq.8
  //print("nom=Ut*covTotInv",Ut*covTotInv);

  TMatrixD denomM= Ut * (covTotInv * U);
  //print("denomM",denomM);
  denomM.Invert(&det);
  print(Form("1/denomM [det=%g]",det),denomM);
  //double denom=denomM(0,0);
  
  /*
  TMatrixD nom_v2 = covTotInv * U;
  //print("nom_v2",nom_v2);

  TMatrixD lambda_v2= denom * nom_v2;
  print("lambda_v2",lambda_v2);

  TMatrixD lambda_v2_tr(TMatrixD::kTransposed,lambda_v2);
  TMatrixD est_v2= lambda_v2_tr * measTot;
  print("estimate_v2",est_v2);
  std::cout << "final_v2 = " << est_v2(0,0) << " +- " << sqrt(denom_v2) << "\n";
  */

  TMatrixD lambda= denomM * nom;
  print("lambda",lambda);

  TMatrixD lambdaU(lambda*U);
  print("(lambda * U) should be ones",lambdaU);

  TMatrixD est= lambda * measTot;
  print("estimate",est);

  if (nMeas==1) {
    std::cout << "measurements: ";
    for (int i=0; i<measA.GetNrows(); i++) {
      std::cout << Form("  %5.2lf +- %4.2lf",measA(i,0),sqrt(covTot(i,i)));
    }
    std::cout << "\n";

    for (int i=0; i<est.GetNrows(); i++) {
      double err=sqrt(denomM(i,i));
      std::cout << "final = " << est(i,0) << " +- " << Form("%4.2lf",err)
		<< "  " << Form("(relErr=%4.2lf%%)",100.*err/est(i,0)) << "\n";
    }
  }
  else {
    std::cout << "lept  measurementA   measurementB   finalValue\n";
    for (int i=0; i<est.GetNrows(); i++) {
      std::cout << ((i==0) ? "(ele)" : "(tau)") << "  ";
      std::cout << Form("%5.2lf +- %4.2lf",measA(i,0),sqrt(covTot(i,i)))
		<< "  "
		<< Form("%5.2lf +- %4.2lf",measB(i,0),sqrt(covTot(i+2,i+2)))
		<< "  "
		<< Form("%5.2lf +- %4.2lf",est(i,0),sqrt(denomM(i,i)))
		<< "\n";
    }
  }

  TMatrixD lambda_tr(lambda);
  lambda_tr.T();
  TMatrixD finCov = lambda * covTot * lambda_tr;
  print("finCov",finCov);

  TMatrixD finCorr(finCov);
  for (int ir=0; ir<finCorr.GetNrows(); ir++) {
    for (int ic=0; ic<finCorr.GetNcols(); ic++) {
      finCorr(ir,ic) = finCov(ir,ic)/ sqrt(finCov(ir,ir)*finCov(ic,ic));
    }
  }
  print("finCorr",finCorr);


  TMatrixD Uest = U*est;
  //print("U*est",Uest);

  TMatrixD meas_diff_est = measTot - U*est;
  print("meas_diff_est",meas_diff_est);

  TMatrixD meas_diff_est_tr(TMatrixD::kTransposed,meas_diff_est);

  TMatrixD chi2= meas_diff_est_tr * covTotInv * meas_diff_est;
  print("chi2",chi2);

  if (nMeas==2) {
    TMatrixD iMeasE(covTot), iMeasT(covTot);
    iMeasE.Zero(); iMeasT.Zero();
    iMeasE(0,0)=1; iMeasE(1,1)=1; // select ele channel
    iMeasT(2,2)=1; iMeasT(3,3)=1; // select tau channel
  
    TMatrixD chi2measE(meas_diff_est_tr*covTotInv*iMeasE*meas_diff_est);
    TMatrixD chi2measT(meas_diff_est_tr*covTotInv*iMeasT*meas_diff_est);
    //print("chi2measE",chi2measE);
    //print("chi2measT",chi2measT);
    std::cout << "chi2=" << chi2(0,0) << " = " << chi2measE(0,0) << " (e) + "
	      << chi2measT(0,0) << " (tau)\n";
  }

  if (nMeas==2) {
    TMatrixD Be_diff_Bt(1,lambda.GetNcols());
    for (int ic=0; ic<lambda.GetNcols(); ic++) {
      Be_diff_Bt(0,ic)=lambda(0,ic)-lambda(1,ic);
    }
    TMatrixD Be_diff_Bt_tr(Be_diff_Bt);
    Be_diff_Bt_tr.T();
    print("Be_diff_Bt_tr",Be_diff_Bt_tr);

    TMatrixD Be_diff_Bt_est(Be_diff_Bt*measTot);
    print("Be_diff_Bt_est",Be_diff_Bt_est);

    TMatrixD diff_unc(Be_diff_Bt*covTot*Be_diff_Bt_tr);
    print("dif_unc",diff_unc);
    std::cout << "the difference (Be-Bt)= " << Be_diff_Bt_est(0,0)
	      << " +- " << sqrt(diff_unc(0,0)) << "\n";

    TMatrixD Be_sum_Bt(1,lambda.GetNcols());
    for (int ic=0; ic<lambda.GetNcols(); ic++) {
      Be_sum_Bt(0,ic)= lambda(0,ic) + lambda(1,ic);
    }
    TMatrixD Be_sum_Bt_tr(Be_sum_Bt);
    Be_sum_Bt_tr.T();
    print("Be_sum_Bt_tr",Be_sum_Bt_tr);

    TMatrixD Be_sum_Bt_est(Be_sum_Bt*measTot);
    //print("Be_sum_Bt_est",Be_sum_Bt_est);

    TMatrixD sum_unc(Be_sum_Bt*covTot*Be_sum_Bt_tr);
    //print("sum_unc",sum_unc);
    std::cout << "the sum (Be+Bt)= " << Be_sum_Bt_est(0,0)
	      << " +- " << sqrt(sum_unc(0,0)) << "\n";
  }

  if (info) {
    std::cout << "\nCase: " << info->caseName << "\n";
    std::cout << "\n";
    for (unsigned int i=0; i<info->size(); i++)
      std::cout << info->at(i) << "\n";
    std::cout << "\n";
  }
}

// ---------------------------------------------------------------
// ---------------------------------------------------------------


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

// ---------------------------------------------------------------

void print(TString msg, double x)
{
  std::cout << msg << "  value=" << x << "\n";
}

// ---------------------------------------------------------------

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

// ---------------------------------------------------------------

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

// ---------------------------------------------------------------
// ---------------------------------------------------------------


