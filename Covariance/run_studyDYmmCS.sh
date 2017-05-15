#!/bin/bash

vars=
vars="_varYieldPoisson"
vars="${vars} _varBkg _varDetRes _varFSRRes _varEff _varAcc"
vars="${vars} _varRhoFile"

nSample=2000
doSave=1

for v in ${vars} ; do
    root -l -q -b studyDYmmCS.C+"(${v},${nSample},${doSave})"
done
