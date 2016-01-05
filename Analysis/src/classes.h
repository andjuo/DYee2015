#include "DYee/Analysis/Include/TEventInfo.hh"
#include "DYee/Analysis/Include/TElectron.hh"
#include "DYee/Analysis/Include/TMuon.hh"

namespace {
  struct dictionary {
    mithep::TEventInfo dummyEventInfo;
    mithep::TElectron dummyElectron;
    mithep::TMuon dummyMuon;
  };
}
