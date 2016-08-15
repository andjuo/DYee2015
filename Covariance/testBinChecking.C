#include "DYbinning.h"
#include "inputs.h"
//#include <iostream>

void testBinChecking() {
  const int ny=4;
  const double ybins[ny+1] = { 0.1, 0.2, 0.4, 0.8, 1.6 };
  TH2D *h2a= new TH2D("h2a","h2a",10,0.4,10.4, ny,ybins);
  TH2D *h2b= new TH2D("h2b","h2b",10,0.4,10.4, ny,ybins);
  TH2D *h2c= new TH2D("h2c","h2c",11,0.4,10.4, ny,ybins);
  TH2D *h2d= new TH2D("h2d","h2d",10,0.4,10.4, ny,ybins);

  std::cout << DYtools::compareBinning(3,h2a,h2b,h2c) << "\n";
  std::cout << DYtools::compareBinning(3,h2a,h2b,h2d) << "\n";
}
