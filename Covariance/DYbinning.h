#ifndef DYbinning
#define DYbinning

#include <TLorentzVector.h>

namespace DYtools {

const Int_t nMassBins = 45;
const Double_t massBinEdges[nMassBins+1] =
  {15, 20, 25, 30, 35, 40, 45, 50, 55, 60,
   64, 68, 72, 76, 81, 86, 91, 96, 101, 106,
   110, 115, 120, 126, 133, 141, 150, 160, 171, 185,
   200, 220, 243, 273, 320, 380, 440, 510, 600, 700,
   830, 1000, 1200, 1500, 2000, 3000};

};

// -------------------------------------------------------------

//int inAcceptance(const TLorentzVector *v1, const TLorentzVector *v2) {
//  return (v1->
//}


// -------------------------------------------------------------

#endif
