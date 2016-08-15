#include "inputs.h"

void testM()
{
  int dim=5;
  TMatrixD m(dim,dim);
  for (int i=0; i<dim; i++)
    for (int j=0; j<dim; j++) {
      m(i,j)= abs(i-dim/2) + abs(j-dim/2);
    }

  double *x= new double[dim+1];
  for (int i=0; i<=dim; i++) x[i]= (100.-0.)*i/double(dim);

  TH1D *h1= new TH1D("h1","h1;x label;y label",dim, x);
  h1->Fill(5.0);
  printHisto(h1);
  plotHisto(h1,"c1");

  const TArrayD *xb= h1->GetXaxis()->GetXbins();
  std::cout << "xb->GetSize=" << xb->GetSize() << "\n";
  for (int i=1; i<xb->GetSize(); i++) {
    std::cout << "i=" << i << ", val=" << xb->At(i) << "\n";
  }

  TH2D *h2=convert2histo(m,h1, "h2","h2");
  TH2D *h2a= new TH2D(m);
  plotHisto(h2,"c2",0,0);
  plotHisto(h2a,"c2a",0,0);
}
