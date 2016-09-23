#ifndef compareVersions_C
#define compareVersions_C

#include "../Covariance/inputs.h"
#include "../Covariance/crossSection.h"
#include "../RooUnfold/src/RooUnfoldResponse.h"
#include <map>

typedef enum { _fn_main=0, _fn_cmp1, _fn_cmp2, _fn_last } TSource_t;

typedef enum { _iEff, _iRho, _iAcc, _iEffAcc,
	       _iSignal, _iUnf, _iUnfRho, _iUnfRhoEffAcc,
	       _iPostFsr, _iPreFsr,
	       _iEffPass, _iEffTot, _iAccPass, _iAccTot,
	       _iDetResFakes, _iFSRResFakes,
	       _iCmp1, _iCmp2,
	       _ih2Mig, _ih2FsrMig
	       //_iRRmeas, _iRRtrue
} THisto_t;

// ---------------------------------------------------------

struct InfoBundle_t {
  TString fileName;
  std::map<THisto_t,TString> *histoNames;
public:
  InfoBundle_t() : fileName(), histoNames(new std::map<THisto_t,TString>()) {}

  InfoBundle_t(TString set_fname, std::map<THisto_t,TString> *set_hNames) :
    fileName(set_fname), histoNames(set_hNames)
  {}

  void add(THisto_t ih, TString hname)
  { (*histoNames)[ih]=hname; }
};

// ---------------------------------------------------------

TString THistoName(THisto_t ih);
TH1D* loadHisto(const std::map<THisto_t,TString> &histoNames,
		THisto_t idxHisto,
		TString fileName, TString newHistoName, TH2D **h2ptr=NULL);

TH1D* loadHisto(THisto_t ih, TString histoDiskName,
		TString fileName, TString newHistoName,
		TH2D **h2ptr);

int drawComparison(const InfoBundle_t &main,
		   const InfoBundle_t &c1,
		   const InfoBundle_t &c2,
		   int noH2Ratios,
		   TString showOnly);

TString lumiScale=":scale2316.969";
TString lumiInvScale=":scale4.31598164844603134e-04";

// ---------------------------------------------------------

void compareVersions(int theCase=1, int version=2, int noH2Ratios_user=0, TString showOnly="")
{
  TString mainFName;
  TString cmpFName1, cmpFName2;
  std::map<THisto_t,TString> histoNames[_fn_last];

  std::vector<TString> vecMainFName,vecCmp1FName,vecCmp2FName;
  addToVector(vecMainFName,". . cs_DYee_13TeV_El2.root");
  addToVector(vecCmp1FName,". . dyee_test_El2skim2_excludeGap.root");
  addToVector(vecCmp2FName,". . dyee_test_RECO_El2skim.root");

  // version 3
  addToVector(vecMainFName,"cs_DYee_13TeV_El2.root");
  addToVector(vecCmp1FName,"dyee_test_dressed_El2skim2.root");
  addToVector(vecCmp2FName,".");

  // version 4
  addToVector(vecMainFName,"cs_DYee_13TeV_El3.root");
  if (1) addToVector(vecCmp1FName,"dyee_test_dressed_El2skim3.root");
  else addToVector(vecCmp1FName,"dyee_test_dressed_El2skim2-old20160916.root");
  addToVector(vecCmp2FName,".");

  // version 5
  addToVector(vecMainFName,"cs_DYee_13TeV_El3.root");
  addToVector(vecCmp1FName,"/mnt/sdc/andriusj/DY13TeV/EgammaWork-20160912/ElectronNtupler/work_76x/Unfolding/RespObj_detUnfolding.root");
  addToVector(vecCmp2FName,".");

  // version 6
  addToVector(vecMainFName,"/mnt/sdc/andriusj/DY13TeV/EgammaWork-20160912/ElectronNtupler/work_76x/Unfolding/dyee_unf_input.root");
  addToVector(vecCmp1FName,"/mnt/sdc/andriusj/DY13TeV/EgammaWork-20160912/ElectronNtupler/work_76x/Unfolding/RespObj_detUnfolding.root");
  addToVector(vecCmp2FName,".");


  if (theCase==1) {
    std::cout << "compare DYee corrections\n";
    mainFName=vecMainFName[version];
    histoNames[_fn_main][_iSignal]="RRmeas+detRes";
    histoNames[_fn_main][_iUnf]="RRtrue+detRes";
    histoNames[_fn_main][_ih2Mig]="RRmig+detRes";
    histoNames[_fn_main][_iDetResFakes]="RRfakes+detRes";
    histoNames[_fn_main][_iEff]="h1Eff";
    histoNames[_fn_main][_iRho]="h1Rho";
    histoNames[_fn_main][_iAcc]="h1Acc";
    histoNames[_fn_main][_iEffAcc]="h1EffAcc";
    histoNames[_fn_main][_iPostFsr]="RRmeas+FSRRes";
    histoNames[_fn_main][_iPreFsr]="RRtrue+FSRRes";
    histoNames[_fn_main][_ih2FsrMig]="RRmig+FSRRes";
    cmpFName1=vecCmp1FName[version];
    if (valueEquals(version,"3 4")) {
      histoNames[_fn_cmp1][_iSignal]="RRmeas+rooUnf_detResRespPU"+lumiScale;
      //histoNames[_fn_cmp1][_iSignal]="h1recoSelNoPostFsrAcc_MWPU"+lumiScale;
      histoNames[_fn_cmp1][_iUnf]="RRtrue+rooUnf_detResRespPU"+lumiScale;
      histoNames[_fn_cmp1][_ih2Mig]="RRmig+rooUnf_detResRespPU"+lumiScale;
      histoNames[_fn_cmp1][_iDetResFakes]="RRfakes+rooUnf_detResRespPU"+lumiScale;
    }
    //histoNames[_fn_cmp1][_iEff]="h1Eff";
    histoNames[_fn_cmp1][_iEff]="h1EffPU";
    if (valueEquals(version,"3 4")) {
      histoNames[_fn_cmp1][_iRho]="h1rho_inPostFsrAcc";
    }
    histoNames[_fn_cmp1][_iAcc]="h1Acc";
    histoNames[_fn_cmp1][_iEffAcc]="h1EffPUAcc";
    histoNames[_fn_cmp1][_iPostFsr]="h1_postFsr_Mweighted"+lumiScale;
    //histoNames[_fn_cmp1][_iPostFsr]="RRmeas+rooUnf_fsrResp"+lumiScale;
    histoNames[_fn_cmp1][_iPreFsr]="RRtrue+rooUnf_fsrResp"+lumiScale;
    histoNames[_fn_cmp1][_ih2FsrMig]="RRmig+rooUnf_fsrResp"+lumiScale;
    cmpFName2=vecCmp2FName[version];
    histoNames[_fn_cmp2][_iRho]="h1rho";
    if (version==2) {
      histoNames[_fn_cmp2][_iSignal]="RRmeas+rooUnf_detResRespPU"+lumiScale;
      histoNames[_fn_cmp2][_iUnf]="RRtrue+rooUnf_detResRespPU"+lumiScale;
    }
  }
  else if ((theCase==2) || (theCase==3)) {
    std::cout << "compare DYee CS closure test\n";
    mainFName="csClosure_DYee_13TeV_El2.root";
    if (theCase==3) {
      std::cout << " -- A x eff replaced by my distribution\n";
      mainFName="csClosure_DYee_13TeV_El2_AJeffAcc.root";
    }
    histoNames[_fn_main][_iSignal]="h1Signal";
    histoNames[_fn_main][_iUnf]="h1Unf";
    histoNames[_fn_main][_iUnfRho]="h1UnfRho";
    histoNames[_fn_main][_iUnfRhoEffAcc]="h1UnfRhoEffAcc";
    cmpFName1="cs_DYee_13TeV_El2.root";
    histoNames[_fn_cmp1][_iSignal]="RRmeas+detRes";
    histoNames[_fn_cmp1][_iUnf]="RRtrue+detRes";
    histoNames[_fn_cmp1][_iUnfRho]="RRtrue+detRes";
    histoNames[_fn_cmp1][_iUnfRhoEffAcc]="RRmeas+FSRRes";
  }
  else if (theCase==4) {
    std::cout << "compare DYee Rho corrections\n";
    mainFName=vecMainFName[version];
    histoNames[_fn_main][_iRho]="h1Rho";
    cmpFName1=vecCmp1FName[version];
    histoNames[_fn_cmp1][_iRho]="h1rho";
    cmpFName2=vecCmp1FName[version];
    histoNames[_fn_cmp2][_iRho]="h1rho_inPostFsrAcc";
  }
  else if (theCase==5) {
    std::cout << "compare DYee Rho corrections\n";
    mainFName=vecMainFName[version];
    histoNames[_fn_main][_iRho]="h1Rho";
    cmpFName1=vecCmp1FName[version];
    histoNames[_fn_cmp1][_iRho]="h1rho_reco";
    cmpFName2=vecCmp1FName[version];
    histoNames[_fn_cmp2][_iRho]="h1rho_recoInPostFsrAcc";
  }
  else if (theCase==6) {
    std::cout << "compare DYee Rho corrections\n";
    mainFName=vecCmp1FName[version];
    histoNames[_fn_main][_iSignal]="RRmeas+rooUnf_detResRespPU";
    cmpFName1=vecCmp1FName[version];
    histoNames[_fn_cmp1][_iSignal]="h1recoSel_MWPU";
  }
  else if (theCase==100) {
    std::cout << "compare detRes objects\n";
    mainFName=vecMainFName[version];
    histoNames[_fn_main][_iSignal]="RRmeas+detRes";
    histoNames[_fn_main][_iUnf]="RRtrue+detRes";
    histoNames[_fn_main][_iDetResFakes]="RRfakes+detRes";
    histoNames[_fn_main][_ih2Mig]="RRmig+detRes";
    cmpFName1=vecCmp1FName[version];
    histoNames[_fn_cmp1][_iSignal]="RRmeas+Unfold_DetectorRes1";
    histoNames[_fn_cmp1][_iUnf]="RRtrue+Unfold_DetectorRes1";
    histoNames[_fn_cmp1][_iDetResFakes]="RRfakes+Unfold_DetectorRes1";
    histoNames[_fn_cmp1][_ih2Mig]="RRmig+Unfold_DetectorRes1";
  }
  else if (theCase==200) {
    std::cout << "compare DYmm preFSR 4pi theory\n";
    mainFName="theory13TeVmm.root";
    histoNames[_fn_main][_iPreFsr]="h1cs_theory";
    lumiScale="";
    cmpFName1="dyee_test_El2skim2_excludeGap.root";
    histoNames[_fn_cmp1][_iPreFsr]="RRtrue+rooUnf_fsrResp"+lumiScale;
    histoNames[_fn_cmp1][_iPreFsr]="h1_preFsr_Mweighted";
    if (1) {
      cmpFName2="cs_DYee_13TeV_El2.root";
      histoNames[_fn_cmp2][_iPreFsr]="RRtrue+FSRRes"+lumiInvScale;
    }
    else {
      cmpFName2="dyee_test_El2skim_excludeGap.root";
      histoNames[_fn_cmp2][_iPreFsr]="RRtrue+rooUnf_fsrResp"+lumiScale;
      histoNames[_fn_cmp2][_iPreFsr]="h1_preFsr_Mweighted";
    }
  }
  else if (theCase==201) {
    std::cout << "compare postFSRInAccSel with and without recoSCGap\n";
    mainFName="dyee_test_dressed_El2skim2.root";
    histoNames[_fn_main][_iUnf]="h1_postFsrInAccSelNoRecoGap_MWPU";
    histoNames[_fn_main][_iUnf]="h1_postFsrInAccSel_MWPU";
    //histoNames[_fn_main][_iEff]="h1_postFsrInAccSelNoRecoGap_MWPU:ratio:h1_postFsrInAcc_MW";
    cmpFName1=mainFName;
    //histoNames[_fn_cmp1][_iUnf]="h1_postFsrInAccSel_MWPU";
    //histoNames[_fn_cmp1][_iEff]="h1EffPU";
    histoNames[_fn_cmp1][_iUnf]="RRtrue+rooUnf_detResRespPU";
  }
  //else if (theCase==4)
  else {
    std::cout << "not ready for theCase=" << theCase << "\n";
  }


  if (mainFName.Length()<2) {
    std::cout << "too short mainFName\n";
    return;
  }

  drawComparison(InfoBundle_t(mainFName,&histoNames[_fn_main]),
		 InfoBundle_t(cmpFName1,&histoNames[_fn_cmp1]),
		 InfoBundle_t(cmpFName2,&histoNames[_fn_cmp2]),
		 noH2Ratios_user,showOnly);

}

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

int drawComparison(const InfoBundle_t &m,
		   const InfoBundle_t &c1,
		   const InfoBundle_t &c2,
		   int noH2Ratios,
		   TString showOnly)
{
  for (std::map<THisto_t,TString>::const_iterator it= m.histoNames->begin();
       it!=m.histoNames->end(); it++) {
    std::cout << THistoName(it->first) << " " << it->second << "\n";    if ((showOnly.Length()>0) && (showOnly!=THistoName(it->first))) {
      std::cout << "   ... skipping by request\n";
      continue;
    }

    TString h1nameBase= "h1" + THistoName(it->first);
    TH2D *h2main=NULL, *h2cmp1=NULL, *h2cmp2=NULL;
    TH1D *h1main= loadHisto(*m.histoNames,it->first,m.fileName,h1nameBase+"_main",&h2main);
    TH1D *h1cmp1=loadHisto(*c1.histoNames,it->first,c1.fileName,h1nameBase+"_cmp1",&h2cmp1);
    TH1D* h1cmp2= loadHisto(*c2.histoNames,it->first,c2.fileName,h1nameBase+"_cmp2",&h2cmp2);

    if (h1main) {
      histoStyle(h1main,kBlue,24,1);
      if (h1cmp1) histoStyle(h1cmp1,kGreen+1,5,1);
      if (h1cmp2) histoStyle(h1cmp2,kRed,7,1);

      TString cName="c" + h1nameBase;
      plotHisto(h1main,cName,1,0,"LPE1","main");
      if (h1cmp1) plotHistoSame(h1cmp1,cName,"LPE1","cmp1");
      if (h1cmp2) plotHistoSame(h1cmp2,cName,"LPE1","cmp2");
      if (1) {
	if (h1cmp1) {
	  std::cout << "mainFile=" << m.fileName
		    << ", cmpFile1=" << c1.fileName << "\n";
	  printRatio(h1main,h1cmp1);
	}
	if (h1cmp2) {
	  std::cout << "mainFile=" << m.fileName
		    << ", cmpFile2=" << c2.fileName << "\n";
	  printRatio(h1main,h1cmp2);
	}
      }
    }
    else if (h2main) {
      int idx=0;
      TH2D *h2_cmp=NULL;
      if (h2cmp1!=NULL) { idx=1; h2_cmp=h2cmp1; }
      else if (h2cmp2!=NULL) { idx=2; h2_cmp=h2cmp2; }
      TString diffTag=Form("_diff%d",idx);
      TString cName="c" + h1nameBase + diffTag;
      TCanvas *cx= new TCanvas(cName,cName,900,300);
      cx->Divide(3,1);
      cx->cd(1);
      h2main->Draw("COLZ");
      if (idx!=0) {
	cx->cd(2);
	h2_cmp->Draw("COLZ");
	TString h2name_new=h2main->GetName()+TString("_diff");
	TH2D *h2main_diff_cmp= cloneHisto(h2main,h2name_new,h2name_new);
	h2main_diff_cmp->Add(h2_cmp,-1);
	cx->cd(3);
	h2main_diff_cmp->Draw("COLZ");
	if (noH2Ratios) std::cout << "skipping h2 ratio print out\n";
	else printHisto(h2main_diff_cmp,1,1, h2main,1e-3);
	std::cout << "h2main(1,1)=" << h2main->GetBinContent(1,1) << "\n";
      }
      cx->Update();
    }
  }
  return 1;
}

// ---------------------------------------------------------
// ---------------------------------------------------------

TString THistoName(THisto_t ih)
{
  TString name;
  switch(ih) {
    // corrections
  case _iEff: name="Eff"; break;
  case _iRho: name="Rho"; break;
  case _iAcc: name="Acc"; break;
  case _iEffAcc: name="EffAcc"; break;
    // distributions
  case _iSignal: name="Signal"; break;
  case _iUnf: name="Unf"; break;
  case _iUnfRho: name="UnfRho"; break;
  case _iUnfRhoEffAcc: name="UnfRhoEffAcc"; break;
  case _iPostFsr: name="postFSR"; break;
  case _iPreFsr: name="preFSR"; break;
  case _iEffPass: name="effPass"; break;
  case _iEffTot: name="effTotal"; break;
  case _iAccPass: name="accPass"; break;
  case _iAccTot: name="accTotal"; break;
  case _iDetResFakes: name="detResFakes"; break;
  case _iFSRResFakes: name="fsrResFakes"; break;
  case _iCmp1: name="cmp1"; break;
  case _iCmp2: name="cmp2"; break;
  case _ih2Mig: name="h2Mig"; break;
  case _ih2FsrMig: name="h2FsrMig"; break;
  default:
    name="UNKNOWN";
  }
  return name;
}

// ---------------------------------------------------------

TH1D* loadHisto(THisto_t ih, TString histoDiskName,
		TString fileName, TString newHistoName,
		TH2D **h2ptr)
{
  std::map<THisto_t,TString> hNames;
  hNames[ih]=histoDiskName;
  return loadHisto(hNames,ih,fileName,newHistoName,h2ptr);
}

// ---------------------------------------------------------

TH1D* loadHisto(const std::map<THisto_t,TString> &histoNames,
		THisto_t idxHisto,
		TString fileName, TString newHistoName,
		TH2D **h2ptr)
{
  TH1D *h1=NULL;
  int is2D=0;
  if (h2ptr) *h2ptr=NULL;
  if (fileName.Length()==0) return h1;
  if (fileName.Length()<2) return h1;

  for (std::map<THisto_t,TString>::const_iterator it=histoNames.begin();
       it!=histoNames.end(); it++) {
    if (it->first == idxHisto) {
      TString hNameDisk=it->second;
      double scale=1.;
      if (hNameDisk.Index(":scale")!=-1) {
	TString valueStr=hNameDisk(hNameDisk.Index(":scale")+6,hNameDisk.Length());
	//std::cout << "valueStr=" << valueStr << "\n";
	if (valueStr.Length()) scale=atof(valueStr.Data());
	hNameDisk.ReplaceAll(":scale" + valueStr,"");
      }
      if (hNameDisk.Index(":ratio:")!=-1) {
	int idx=hNameDisk.Index(":ratio:");
	TString hist1name= hNameDisk(0,idx);
	TString hist2name= hNameDisk(idx+7,hNameDisk.Length());
	std::cout << "histnames=<" << hist1name << "> and <"
		  << hist2name << ">\n";
	h1= loadHisto(fileName,hist1name, newHistoName,h1dummy);
	TH1D *denom= loadHisto(fileName,hist2name, newHistoName+TString("_denom"),h1dummy);
	if (!denom) {
	  std::cout << "failed to load " << hist2name << " from file "
		    << fileName << "\n";
	  if (h1) { delete h1; h1=NULL; }
	}
	else {
	  if (h1) h1->Divide(denom);
	}
      }
      else if (hNameDisk.Index(":unf")!=-1) {
	int idx=hNameDisk.Index(":unf");
	int nIters=atoi(hNameDisk.Data()+idx+4);
	TString histname= hNameDisk(0,idx);
	TString respName= hNameDisk(hNameDisk.Index(":",idx+1),hNameDisk.Length());
	std::cout << ":unf: histname=" << histname << ", respName=" << respName << ", nIters=" << nIters << "\n";
	return NULL;
      }
      else if (hNameDisk.Index("+")==-1) {
	h1= loadHisto(fileName,hNameDisk, newHistoName,h1dummy);
      }
      else if ((hNameDisk.Index("RRmeas+")!=-1) ||
	       (hNameDisk.Index("RRtrue+")!=-1) ||
	       (hNameDisk.Index("RRmig+")!=-1) ||
	       (hNameDisk.Index("RRfakes+")!=-1)) {
	TString respName=hNameDisk;
	respName.Remove(0,respName.Index("+")+1);
	std::cout << "respName=" << respName << "\n";
	RooUnfoldResponse *r= loadRooUnfoldResponse(fileName,respName,"resp"+newHistoName);
	if (!r) {
	  std::cout << "failed to get response " << respName << " from "
		    << fileName << "\n";
	}
	else {
	  if (hNameDisk.Index("RRmeas+")!=-1)
	    h1=cloneHisto((TH1D*)r->Hmeasured(),newHistoName,newHistoName);
	  else if (hNameDisk.Index("RRtrue+")!=-1)
	    h1=cloneHisto((TH1D*)r->Htruth(),newHistoName,newHistoName);
	  else if (hNameDisk.Index("RRfakes+")!=-1)
	    h1=cloneHisto((TH1D*)r->Hfakes(),newHistoName,newHistoName);
	  else if (hNameDisk.Index("RRmig+")!=-1) {
	    is2D=1;
	    if (!h2ptr) { std::cout << "cannot receive RRmig\n"; }
	    else *h2ptr=cloneHisto((TH2D*)r->Hresponse(),newHistoName,newHistoName);
	  }
	  delete r;
	}
      }
      if ((!h1 && !is2D) || (is2D && (!h2ptr || !(*h2ptr)))) {
	std::cout << "failed to get " << it->second << "\n";
      }
      else if (scale!=1.) {
	if (h1) h1->Scale(scale);
	else (*h2ptr)->Scale(scale);
      }
      break;
    }
  }
  return h1;
}

// ---------------------------------------------------------

#endif
