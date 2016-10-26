#include "inputs.h"
#include "crossSection.h"
#include "DYmm13TeV_eff.h"

  const int nEtaBins=2;
  const double etaArr[nEtaBins+1]= { -2.5, 0., 2.5 };
  //const int nPtBins=3;
  //const double ptArr[nPtBins+1]= { 0., 10., 20., 50. };
  const int nPtBins=1;
  const double ptArr[nPtBins+1] = { 0., 100. };

// ----------------------------------------------------------

void setValue(TH2D *h2, int ibin, int jbin, double v, double dv)
{
  h2->SetBinContent(ibin,jbin, v);
  h2->SetBinError  (ibin,jbin,dv);
}

// ----------------------------------------------------------

//std::vector<TH1D*>* randomizeSF(const DYTnPEffColl_t &eff, const EventSpace_t &es,
//				int nToys, TString nameBase);

DYTnPEff_t* randomizeEffs(const DYTnPEffColl_t &coll, int srcIdx, int nToys, TString tag);

// ----------------------------------------------------------

void testDYTnPEffColl(int nExps_input=100)
{

  TH2D *h2def= new TH2D("h2def","h2def;eta;pt", nEtaBins,etaArr,nPtBins,ptArr);
  h2def->SetStats(0);

  for (int ibin=1; ibin<=nEtaBins; ibin++) {
    for (int jbin=1; jbin<=nPtBins; jbin++) {
      //std::cout << "ibin=" << ibin << ", jbin=" << jbin << "\n";
      h2def->SetBinContent(ibin,jbin, 1.);
      h2def->SetBinError  (ibin,jbin, 0.);
    }
  }
  printHisto(h2def);

  DYTnPEff_t tnpDef;
  tnpDef.h2Eff_RecoID_Data= cloneHisto(h2def,"h2eff_recoID_data","h2eff_recoID_data");
  tnpDef.h2Eff_Iso_Data= cloneHisto(h2def,"h2eff_iso_data","h2eff_iso_data");
  tnpDef.h2Eff_HLT4p2_Data= cloneHisto(h2def,"h2eff_hlt4p2_data","h2eff_hlt4p2_data");
  tnpDef.h2Eff_HLT4p3_Data= cloneHisto(h2def,"h2eff_hlt4p3_data","h2eff_hlt4p3_data");
  tnpDef.h2Eff_RecoID_MC= cloneHisto(h2def,"h2eff_recoID_mc","h2eff_recoID_mc");
  tnpDef.h2Eff_Iso_MC= cloneHisto(h2def,"h2eff_iso_mc","h2eff_iso_mc");
  tnpDef.h2Eff_HLT4p2_MC= cloneHisto(h2def,"h2eff_hlt4p2_mc","h2eff_hlt4p2_mc");
  tnpDef.h2Eff_HLT4p3_MC= cloneHisto(h2def,"h2eff_hlt4p3_mc","h2eff_hlt4p3_mc");

  DYTnPEff_t eStat(tnpDef,"stat");
  DYTnPEff_t eSystA(tnpDef,"systA");
  DYTnPEff_t eSystB(tnpDef,"systB");

  setValue(eStat.h2Eff_RecoID_Data,1,1, 0.5,0.5);
  setValue(eStat.h2Eff_RecoID_Data,2,1, 1.25,0.75);
  setValue(eSystB.h2Eff_Iso_MC,2,1, 0.75, 100.);


  DYTnPEffColl_t coll(0);
  if (!coll.assignStatErr(eStat,"_coll") ||
      !coll.addSystErrSource(eSystA,"_collA") ||
      !coll.addSystErrSource(eSystB,"_collB")) {
    std::cout << "error adding stat and syst err\n";
    return;
  }

  coll.listNumbers();
  //coll.displayAll();

  //DYTnPEff_t *eTot= coll.randomize(9999,"_total");
  //eTot->listNumbers();
  //eTot->displayAll();
  //eTot->displayEffTot(1+2,"_total",0,1);

  if (0) {
    DYTnPEff_t *effSystA_static= coll.getTnPWithSystUnc(0);
    DYTnPEff_t *effSystB_static= coll.getTnPWithSystUnc(1);
    effSystA_static->listNumbers();
    effSystB_static->listNumbers();
    return;
  }

  int nToys=nExps_input;

  if (0) {
    coll.getTnPWithStatUnc().displayEffTot(1+2,"",-1,DYTnPEff_t::_misc_hlt4p2_leadLep);
    DYTnPEff_t *effStat=randomizeEffs(coll,0, nToys,"_statOnly");
    effStat->displayAll();
    effStat->listNumbers();
    return;
  }

  if (0) {
    //DYTnPEff_t *effSystA= coll.getTnPWithSystUnc(0);
    //effSystA->displayEffTot(1+2,"_systA",-1,DYTnPEff_t::_misc_hlt4p2_leadLep);
    //effSystA->listNumbers();
    DYTnPEff_t *effSystB_static= coll.getTnPWithSystUnc(1);
    effSystB_static->displayEffTot(1+2,"_systB_static",-1,DYTnPEff_t::_misc_hlt4p2_leadLep);
    effSystB_static->listNumbers();
    DYTnPEff_t *effSystB= randomizeEffs(coll,2, nToys,"_systBOnly");
    effSystB->displayAll();
    effSystB->listNumbers();
    return;
  }

  if (1) {
    DYTnPEff_t *effTot= coll.getTnPWithTotUnc();
    effTot->listNumbers();
    DYTnPEff_t *effTotRnd= randomizeEffs(coll,-1111, nToys,"_all");
    effTotRnd->listNumbers();
    effTotRnd->displayAll();
  }

    return;
  randomizeEffs(coll,-1,nToys,"_systAOnly");
  randomizeEffs(coll,-2,nToys,"_systBOnly");

  return;


  EventSpace_t esTest("esTest",h2def);
  TLorentzVector mass[4],v1[4];
  double m=DYtools::massBinCenter(0);
  mass[0].SetPtEtaPhiE(2.,2.2,1.3,m);
  v1[0].SetPtEtaPhiE(1.,0.8,0.4,0.5*m);
  m=DYtools::massBinCenter(1);
  mass[1].SetPtEtaPhiE(10.,-1.2,1.3,m);
  v1[1].SetPtEtaPhiE(1.,-0.8,0.4,0.5*m);
  m=DYtools::massBinCenter(2);
  mass[2].SetPtEtaPhiE(10.,-1.2,1.3,m);
  v1[2].SetPtEtaPhiE(1.,2.2,0.4,0.8*m);
  mass[3].SetPtEtaPhiE(10., 1.2,1.3,m);
  v1[3].SetPtEtaPhiE(1.,0.8,0.4,0.5*m);
  for (int im=0; im<4; im++) {
    TLorentzVector v2= mass[im] - v1[im];
    std::cout << "v1=" << v1[im] << ", v2=" << v2 << ", mass=" << mass[im]
	      << ", " << "(v1+v2)=" << (v1[im]+v2)
      //<< ", (v1+v2-m)=" << (v1[im]+v2-mass[im])
	      << "\n";
    esTest.fill(v1+im, &v2, 1.);
  }
  esTest.displayAll();

  //deriveUn
}

// ----------------------------------------------------------
// ----------------------------------------------------------

DYTnPEff_t* randomizeEffs(const DYTnPEffColl_t &coll, int srcIdx, int nToys, TString tagInp)
{
  std::vector<DYTnPEff_t*> rndV;
  rndV.reserve(nToys);
  TH2D *h2chk= NULL;
  if (0) {
    h2chk=new TH2D("h2chk","h2chk",2,0.,2.,100,-0.5,2.5);
    h2chk->SetDirectory(0);
    h2chk->SetStats(0);
  }
  for (int itoy=0; itoy<nToys; itoy++) {
    TString tag= tagInp + Form("_%d",itoy);
    DYTnPEff_t* eff= coll.randomize(srcIdx,tag);
    rndV.push_back(eff);

    //printHisto(eff->h2Eff_RecoID_Data);
    if (h2chk) {
      for (int i=1; i<=2; i++) {
	//std::cout << "fill @ " << i+0.2 << " : " << eff->h2Eff_RecoID_Data->GetBinContent(i,1) << "\n";
	h2chk->Fill(i-0.2, eff->h2Eff_RecoID_Data->GetBinContent(i,1));
      }
    }
  }
  if (h2chk) {
    plotHisto(h2chk,"ch2chk_recoId_data",0,0,"COLZ");
    printHisto(h2chk);
  }

  // calculate average
  DYTnPEff_t *effAvg= new DYTnPEff_t(coll.getTnPWithStatUnc(),
				     Form("effAvg_srcIdx%d",srcIdx));
  effAvg->resetAll();
  for (int ikind=0; ikind<3; ikind++) {
    std::vector<TH2D*> *h2avgV= effAvg->h2vecPtr(ikind);
    for (unsigned int itoy=0; itoy<rndV.size(); itoy++) {
      const std::vector<TH2D*> *h2rndV= rndV[itoy]->h2vecPtr(ikind);
      for (unsigned int i=0; i<h2rndV->size(); i++) {
	(*h2avgV)[i] -> Add( (*h2rndV)[i], 1/double(nToys) );
      }
    }
  }

  // calculate variance
  DYTnPEff_t effVar(coll.getTnPWithStatUnc(),Form("effVar_srcIdx%d",srcIdx));
  effVar.resetAll();
  for (int ikind=0; ikind<3; ikind++) {
    std::vector<TH2D*> *h2varV= effVar.h2vecPtr(ikind);
    // accumulate
    const std::vector<TH2D*> *h2avgV= effAvg->h2vecPtr(ikind);
    for (unsigned int itoy=0; itoy<rndV.size(); itoy++) {
      const std::vector<TH2D*> *h2rndV= rndV[itoy]->h2vecPtr(ikind);
      for (unsigned int i=0; i<h2varV->size(); i++) {
	TH2D *h2var= (*h2varV)[i];
	const TH2D *h2rnd= (*h2rndV)[i];
	const TH2D *h2avg= (*h2avgV)[i];
	for (int ibin=1; ibin<=h2var->GetNbinsX(); ibin++) {
	  for (int jbin=1; jbin<=h2var->GetNbinsY(); jbin++) {
	    double dev=
	      h2rnd->GetBinContent(ibin,jbin) - h2avg->GetBinContent(ibin,jbin);
	    h2var->SetBinContent(ibin,jbin,
		 h2var->GetBinContent(ibin,jbin) +  pow(dev,2)/double(nToys));
	  }
	}
      }
    }
  }

  // calculate and set the error
  for (int ikind=0; ikind<3; ikind++) {
    std::vector<TH2D*> *h2avgV= effAvg->h2vecPtr(ikind);
    const std::vector<TH2D*> *h2varV= effVar.h2vecPtr(ikind);
    for (unsigned int i=0; i<h2varV->size(); i++) {
      TH2D *h2avg= (*h2avgV)[i];
      const TH2D *h2var= (*h2varV)[i];
      for (int ibin=1; ibin<=h2avg->GetNbinsX(); ibin++) {
	for (int jbin=1; jbin<=h2avg->GetNbinsY(); jbin++) {
	  double errSqr= h2var->GetBinContent(ibin,jbin);
	  double err=0.;
	  if (errSqr<0) {
	    std::cout << "ibin=" << ibin << ", jbin=" << jbin
		      << ", negative err^2=" << errSqr << "\n";
	    err=sqrt(-errSqr);
	  }
	  else err=sqrt(errSqr);
	  h2avg->SetBinError(ibin,jbin, err);
	}
      }
    }
  }

  return effAvg;
}

// ----------------------------------------------------------
// ----------------------------------------------------------
/*
std::vector<TH1D*>* randomizeSF(const DYTnPEffColl_t &coll,
				int rndSrcIdx,
				const EventSpace_t &es,
				int nToys, TString nameBase)
{
  std::vector<TH1D*> *h1V= new std::vector<TH1D*>();
  for (int iToy=0; iToy<nToys; iToy++) {
    //DYTnPEff_t *rndEff= coll.randomize(srcIdx,
    //rndEff.randomize(
    //TH1D *h1toy= es.calculateScaleFactor(eff
  }
  return NULL;
}
*/
// ----------------------------------------------------------
