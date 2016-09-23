#include "Blue.h"

// ---------------------------------------------------------------
// ---------------------------------------------------------------

// more relies on BLUEResult_t
void combineX(const TMatrixD &measA_inp, const TMatrixD &measB_inp,
	      const TMatrixD &covStat_inp, const TMatrixD &covSyst_inp,
	      const Info_t *info);

// ---------------------------------------------------------------
// ---------------------------------------------------------------


void myTest2_StatSyst(int theCase, int printInfo=0)
{
  TMatrixD *measAinp=NULL, *measBinp=NULL;
  TMatrixD *covStat=NULL, *covSyst=NULL;
  int nMeas=2;
  TString caseName="unknown";
  Info_t info;
  info.add(Form("myBlue_StatSyst: Case number in the macro theCase=%d",theCase));

  if ((theCase==1) || (theCase==2) || (theCase==3) || (theCase==4)
      ||
      (theCase==5) || (theCase==6) || (theCase==7) || (theCase==8)
      ) {
    caseName="4.3.3. Breakdown of error contributions";
    if (theCase==1) caseName.Append(" (uncorrelated, eq.54)");
    else if (theCase==2) caseName.Append(" (correlated, eq.55)");
    else if (theCase==3) caseName.Append(" (anti-correlated, eq.53)");
    else if (theCase==4) caseName.Append(". No syst.err.; eq.56");
    else if ((theCase>=5) && (theCase<=8)) {
      caseName.Append(". Lepton universality");
      if (theCase==5) caseName.Append(". Correlated, eq.57");
      else if (theCase==6) caseName.Append(". Uncorrelated, eq.58");
      else if (theCase==7) caseName.Append(". Unticorrelated, eq.59");
      else if (theCase==8) caseName.Append(". No systematic err, eq.60");
    }
    nMeas=2;
    measAinp= new TMatrixD(2,1);
    measBinp= new TMatrixD(2,1);
    assignValuesT(*measAinp,0,"10.50 9.50");
    assignValuesT(*measBinp,0,"13.50 14.00");
    covStat= new TMatrixD(4,4);
    covSyst= new TMatrixD(4,4);
    info.add("Eq.(51) values:");
    info.add("B_A^e=10.50 +- 1.00(stat)");
    info.add("B_B^e=13.50 +- 0.21(stat) +- 2.99(syst)");
    info.add("B_A^t= 9.50 +- 3.00(stat)");
    info.add("B_B^t=14.00 +- 0.21(stat) +- 2.99(syst)");
    info.add("To get the previous case, stat^2+syst^2=9");
    assignValues(*covStat,0,"1.00 0.     0.   0.   ");
    assignValues(*covStat,1,"0.   0.0441 0.   0.   ");
    assignValues(*covStat,2,"0.   0.     9.00 0.   ");
    assignValues(*covStat,3,"0.   0.     0.   0.0441");
    assignValues(*covSyst,0,"0.   0.      0.   0.   ");
    assignValues(*covSyst,1,"0.   8.9559  0.   8.9559");
    assignValues(*covSyst,2,"0.   0.      0.   0.   ");
    assignValues(*covSyst,3,"0.   8.9559  0.   8.9559");
    if (0)
    for (int ir=0; ir<covStat->GetNrows(); ir++) {
      for (int ic=0; ic<covStat->GetNcols(); ic++) {
	(*covStat)(ir,ic)= 0.01*trunc(100.*(*covStat)(ir,ic));
	(*covSyst)(ir,ic)= 0.01*trunc(100.*(*covSyst)(ir,ic));
      }
    }
    if ((theCase==1) || (theCase==6)) {
      (*covSyst)(1,3)=0;
      (*covSyst)(3,1)=0;
    }
    if ((theCase==3) || (theCase==7)) {
      (*covSyst)(1,3)=-(*covSyst)(1,3);
      (*covSyst)(3,1)=-(*covSyst)(3,1);
    }
    if ((theCase==4) || (theCase==8)) {
      covSyst->Zero();
    }
    if ((theCase>=5) && (theCase<=8)) {
      delete measAinp;
      delete measBinp;
      measAinp= new TMatrixD(4,1);
      measBinp= new TMatrixD(0,0);
      assignValuesT(*measAinp,0,"10.50 13.50  9.50 14.00");
    }
  }
  else {
    std::cout << "not ready for case=" << theCase << "\n";
    return;
  }

  info.description=caseName;

  combineX(*measAinp,*measBinp,*covStat,*covSyst,(printInfo) ? &info : NULL);
  std::cout << "caseName=" << caseName << "\n";
}

// ---------------------------------------------------------------
// -------------------------------------------------------------

// combine uncorrelated measurements
void combineX(const TMatrixD &measA_inp, const TMatrixD &measB_inp,
	      const TMatrixD &covStat_inp, const TMatrixD &covSyst_inp,
	      const Info_t *info)
{

  int nMeas= (measB_inp.GetNrows()==0) ? 1 : 2;
  TMatrixD covTot_inp(covStat_inp,TMatrixD::kPlus,covSyst_inp);

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
  
  TMatrixD finCovStat= blue.contributedCov(covStat_inp);
  TMatrixD finCovSyst= blue.contributedCov(covSyst_inp);
  printCov("finCovStat",finCovStat);
  printCov("finCovSyst",finCovSyst);


  if (nMeas==1) {
    std::cout << "measurements: ";
    for (int i=0; i<measA.GetNrows(); i++) {
      std::cout << Form("  %5.2lf +- %4.2lf",blue.getMeas(i),blue.getInpErr(i));
    }
    std::cout << "\n";

    for (int i=0; i<blue.getEst()->GetNrows(); i++) {
      double err=blue.getErr(i);
      std::cout << "final = " << blue.getEst(i) << " +- "
		<< Form("%4.2lf",err) << "  "
		<< " = " << blue.getEst(i)
		<< Form(" +- %4.2lf(stat) +- %4.2lf(syst)",sqrt(finCovStat(i,i)),sqrt(finCovSyst(i,i)))
		<< "\n";
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
		<< Form("%5.2lf +- %4.2lf",blue.getEst(i),blue.getErr(i))
		<< "  "
		<< Form("%5.2lf +- %4.2lf(stat) +- %4.2lf(syst)",blue.getEst(i),sqrt(finCovStat(i,i)),sqrt(finCovSyst(i,i)))
		<< "\n";
    }
  }

  printCov("finCov",blue.getCov());
  //print("finCorr",blue.getCorr());
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
