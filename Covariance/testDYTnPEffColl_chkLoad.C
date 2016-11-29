#include "DYmm13TeV_eff.h"

void testDYTnPEffColl_chkLoad()
{
  DYTnPEffColl_t collChk(0);
  if (!collChk.load("data_test_TnPEffColl.root","","_coll _collA _collB")) {
    std::cout << "loading failed\n";
  }
  else {
    std::cout << "loading ok\n";
    collChk.listNumbers(1);
  }
  return;
}
