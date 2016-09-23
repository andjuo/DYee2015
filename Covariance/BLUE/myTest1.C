// Examples from Valassi 2003

#include <TMatrixD.h>
#include <iostream>
#include <math.h>
#include <sstream>
#include <TString.h>
#include "Blue.h"

// ---------------------------------------------------------------
// ---------------------------------------------------------------

// combine uncorrelated measurements
void combine(const TMatrixD &measA_inp, const TMatrixD &measB_inp,
	     const TMatrixD &covTot_inp, const Info_t *info);
// more relies on BLUEResult_t
void combineX(const TMatrixD &measA_inp, const TMatrixD &measB_inp,
	      const TMatrixD &covTot_inp, const Info_t *info);

// ---------------------------------------------------------------
// ---------------------------------------------------------------

void myTest1(int theCase=1, int printInfo=0)
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
    if (theCase==10) caseName+=". eq.48."; else caseName+=". eq.50.";
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
  info.description=caseName;

  TMatrixD measA(*measAinp);
  TMatrixD measB(*measBinp);
  TMatrixD covTot(*covInp);
  print("measA",measA);
  print("measB",measB);

  //combine(measA,measB,covTot,(printInfo) ? &info : NULL);
  combineX(measA,measB,covTot,(printInfo) ? &info : NULL);
  std::cout << "caseName=" << caseName << "\n";
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


  BLUEResult_t blue(measTot,U,covTot);
  if (!blue.estimate()) {
    std::cout << "combination failed\n";
    return;
  }


  if (nMeas==1) {
    std::cout << "measurements: ";
    for (int i=0; i<measA.GetNrows(); i++) {
      std::cout << Form("  %5.2lf +- %4.2lf",measA(i,0),sqrt(covTot(i,i)));
    }
    std::cout << "\n";

    for (int i=0; i<blue.getEst()->GetNrows(); i++) {
      double err=sqrt(blue.getCov(i,i));
      std::cout << "final = " << blue.getEst(i,0) << " +- "
		<< Form("%4.2lf",err) << "  "
		<< Form("(relErr=%4.2lf%%)",100.*err/blue.getEst(i,0)) << "\n";
    }
  }
  else {
    std::cout << "lept  measurementA   measurementB   finalValue\n";
    for (int i=0; i<blue.getEst()->GetNrows(); i++) {
      std::cout << ((i==0) ? "(ele)" : "(tau)") << "  ";
      std::cout << Form("%5.2lf +- %4.2lf",measA(i,0),sqrt(covTot(2*i,2*i)))
		<< "  "
		<< Form("%5.2lf +- %4.2lf",measB(i,0),sqrt(covTot(2*i+1,2*i+1)))
		<< "  "
		<< Form("%5.2lf +- %4.2lf",blue.getEst(i,0),sqrt(blue.getCov(i,i)))
		<< "\n";
    }
  }

  print("finCov",blue.getCov());
  print("finCorr",blue.getCorr());
  print("chi2",blue.getChi2());

  if (nMeas==2) {
    std::cout << "chi2=" << blue.getChi2() << " = " << blue.getChi2(0)
	      << " (e) + " << blue.getChi2(1) << " (tau)\n";
  }

  if (nMeas==2) {
    print("contributions to diff",blue.combinedLambda(0,1,-1.));

    TMatrixD Be_diff_Bt_est = blue.combinedEst(0,1, -1.);
    //print("Be_diff_Bt_est",Be_diff_Bt_est);

    TMatrixD diff_unc = blue.combinedCov(0,1, -1.);
    //print("dif_unc",diff_unc);
    std::cout << "the difference (Be-Bt)= " << Be_diff_Bt_est(0,0)
	      << " +- " << sqrt(diff_unc(0,0)) << "\n";

    print("contributions to sum",blue.combinedLambda(0,1,1.));

    TMatrixD Be_sum_Bt_est = blue.combinedEst(0,1, 1.);
    //print("Be_sum_Bt_est",Be_sum_Bt_est);

    TMatrixD sum_unc = blue.combinedCov(0,1, 1.);
    //print("sum_unc",sum_unc);
    std::cout << "the sum (Be+Bt)= " << Be_sum_Bt_est(0,0)
	      << " +- " << sqrt(sum_unc(0,0)) << "\n";
  }

  if (info) info->print();
}


// ---------------------------------------------------------------
// -------------------------------------------------------------

// combine uncorrelated measurements
void combineX(const TMatrixD &measA_inp, const TMatrixD &measB_inp,
	      const TMatrixD &covTot_inp, const Info_t *info)
{

  int nMeas= (measB_inp.GetNrows()==0) ? 1 : 2;
  BLUEResult_t blue;

  if (nMeas==1) {
    if (!blue.estimate(measA_inp,covTot_inp)) {
      std::cout << "estimation for 1 parameter combination failed\n";
      return;
    }
  }
  else {
    int dim1= measA_inp.GetNrows();
    TMatrixD covA(dim1,dim1);
    TMatrixD covB(dim1,dim1);
    TMatrixD covAB(dim1,dim1);
    for (int ir=0; ir<dim1; ir++) {
      for (int ic=0; ic<dim1; ic++) {
	covA(ir,ic)= covTot_inp(2*ir  , 2*ic);
	covB(ir,ic)= covTot_inp(2*ir+1, 2*ic+1);
	covAB(ir,ic)=covTot_inp(2*ir, 2*ic+1);
      }
    }
    print("covA",covA);
    print("covB",covB);
    print("covAB",covAB);

    if (!blue.estimateTest(measA_inp, covA, measB_inp, covB, covAB)) {
      std::cout <<" estimation of 2 measurement combination failed\n";
      return;
    }
  }

  if (!blue.getEst()) {
    std::cout << "estimation is null ptr\n";
    return;
  }

  TMatrixD measA(measA_inp);
  TMatrixD measB(measB_inp);
  

  if (nMeas==1) {
    std::cout << "measurements: ";
    for (int i=0; i<measA.GetNrows(); i++) {
      std::cout << Form("  %5.2lf +- %4.2lf",measA(i,0),sqrt(covTot_inp(i,i)));
    }
    std::cout << "\n";

    for (int i=0; i<blue.getEst()->GetNrows(); i++) {
      double err=sqrt(blue.getCov(i,i));
      std::cout << "final = " << blue.getEst(i,0) << " +- "
		<< Form("%4.2lf",err) << "  "
		<< Form("(relErr=%4.2lf%%)",100.*err/blue.getEst(i,0)) << "\n";
    }
  }
  else {
    std::cout << "lept  measurementA   measurementB   finalValue\n";
    for (int i=0; i<blue.getEst()->GetNrows(); i++) {
      std::cout << ((i==0) ? "(ele)" : "(tau)") << "  ";
      std::cout << Form("%5.2lf +- %4.2lf",measA(i,0),sqrt(covTot_inp(2*i,2*i)))
		<< "  "
		<< Form("%5.2lf +- %4.2lf",measB(i,0),sqrt(covTot_inp(2*i+1,2*i+1)))
		<< "  "
		<< Form("%5.2lf +- %4.2lf",blue.getEst(i,0),sqrt(blue.getCov(i,i)))
		<< "\n";
    }
  }

  print("finCov",blue.getCov());
  print("finCorr",blue.getCorr());
  print("chi2",blue.getChi2());

  if (nMeas==2) {
    std::cout << "chi2=" << blue.getChi2() << " = " << blue.getChi2(0)
	      << " (e) + " << blue.getChi2(1) << " (tau)\n";
  }

  if (nMeas==2) {
    print("contributions to diff",blue.combinedLambda(0,1,-1.));

    TMatrixD Be_diff_Bt_est = blue.combinedEst(0,1, -1.);
    //print("Be_diff_Bt_est",Be_diff_Bt_est);

    TMatrixD diff_unc = blue.combinedCov(0,1, -1.);
    //print("dif_unc",diff_unc);
    std::cout << "the difference (Be-Bt)= " << Be_diff_Bt_est(0,0)
	      << " +- " << sqrt(diff_unc(0,0)) << "\n";

    print("contributions to sum",blue.combinedLambda(0,1,1.));

    TMatrixD Be_sum_Bt_est = blue.combinedEst(0,1, 1.);
    //print("Be_sum_Bt_est",Be_sum_Bt_est);

    TMatrixD sum_unc = blue.combinedCov(0,1, 1.);
    //print("sum_unc",sum_unc);
    std::cout << "the sum (Be+Bt)= " << Be_sum_Bt_est(0,0)
	      << " +- " << sqrt(sum_unc(0,0)) << "\n";
  }

  if (info) info->print();
}


// ---------------------------------------------------------------
// ---------------------------------------------------------------


