#include "inputs.h"
#include "DYbinning.h"

void convertUncTable(int elChannel=0)
{
  
  TString fname=(elChannel)?"dyee_uncertainties.txt" : "dymm_uncertainties.txt";
  std::ifstream fin(fname);
  if (!fin.is_open()) {
    std::cout << "failed to open the file <" << fname << ">\n";
    return;
  }

  TString foutName="dymm_unc_fromNote.root";
  if (elChannel) foutName.ReplaceAll("dymm","dyee");
  TString h1NameBase=(elChannel) ? "h1csEE" : "h1csMM";

  TH1D* h1cs= new TH1D(h1NameBase,h1NameBase,DYtools::nMassBins43,DYtools::massBinEdges43);
  TH1D* h1csStatUnc= new TH1D(h1NameBase + "_statUnc",h1NameBase + "_statUnc",DYtools::nMassBins43,DYtools::massBinEdges43);
  TH1D* h1csExpUnc= new TH1D(h1NameBase + "_expUnc",h1NameBase + "_expUnc",DYtools::nMassBins43,DYtools::massBinEdges43);
  TH1D* h1csTheoUnc= new TH1D(h1NameBase + "_theoUnc",h1NameBase + "_theoUnc",DYtools::nMassBins43,DYtools::massBinEdges43);
  TH1D* h1csTotUnc= new TH1D(h1NameBase + "_totUnc",h1NameBase + "_totUnc",DYtools::nMassBins43,DYtools::massBinEdges43);

  std::vector<TH1D*> h1V;
  h1V.push_back(h1cs);
  h1V.push_back(h1csStatUnc);
  h1V.push_back(h1csExpUnc);
  h1V.push_back(h1csTheoUnc);
  h1V.push_back(h1csTotUnc);


  std::string line;
  TString s;
  int im=-1;
  const int nVals=5;
  double value[2*nVals];
  while (!fin.eof() && getline(fin,line)) {
    std::cout << "got line <" << line << ">\n";
    std::stringstream ss(line);
    int idx=0, iv=0;
    while (!ss.eof()) {
      ss >> s;
      std::cout << "idx=" << idx << ", s[" << s.Length() << "]=<" << s << ">\n";
      switch(idx) {
      case 0: break; // range
      case 1: case 4: case 7: case 10: case 13:
	// mantisse
	value[iv]=atof(s.Data());
	break;
      case 2:  case 5: case 8: case 11: case 14:
	// times
	break;
      case 3: case 6: case 9: case 12: case 15:
	// exponent
	{
	  s.ReplaceAll("10??","");
	  double x=1.;
	  if (s=="102") x=100;
	  else if (s=="101") x=10;
	  else if (s=="100") x=1;
	  else if (s=="1") x=1e-1;
	  else if (s=="2") x=1e-2;
	  else if (s=="3") x=1e-3;
	  else if (s=="4") x=1e-4;
	  else if (s=="5") x=1e-5;
	  else if (s=="6") x=1e-6;
	  else if (s=="7") x=1e-7;
	  else if (s=="8") x=1e-8;
	  else if (s=="9") x=1e-9;
	  else if (s=="10") x=1e-10;
	  else {
	    std::cout << "not ready for s=" << s << "\n";
	  }
	  std::cout << "iv=" << iv << ", set value=" << value[iv]
		    << " * " << x << " = " << (value[iv]*x) << "\n";
	  value[iv]*=x;
	  iv++;
	}
	break;
      default: ; // nothing
      }
      idx++;
    }
    std::cout << "was line=<" << line << ">\n";
    std::cout << "im=" << im << ", numbers: ";
    for (int i=0; i<nVals; i++) std::cout << Form(" %7.2g",value[i]);
    std::cout << "\n";
    if (im>0) {
      for (unsigned int i=0; i<h1V.size(); i++) {
	h1V[i]->SetBinContent(im,value[i]);
      }
    }
    if (im==3) break;
    im++;
  }
  fin.close();
  return;

  //TString foutName="dyLL_unc_fromNote.root";
  TFile fout(foutName,"RECREATE");
  for (unsigned int i=0; i<h1V.size(); i++)
    h1V[i]->Write();
  writeTimeTag(&fout);
  fout.Close();
  std::cout << "file <" << fout.GetName() << "> created\n";
  return;
}
