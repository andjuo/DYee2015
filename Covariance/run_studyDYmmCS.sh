#!/bin/bash

vars=
vars="_varYieldPoisson"
vars="${vars} _varBkg _varDetRes _varFSRRes _varEff _varAcc"
vars="${vars} _varRhoFile"

nSample=2000
doSave=1

code="studyDYxxCS.C"
if [ ${#1} -gt 0 ] ; then
    if [ "$1" == "ee" ] ; then
	code="studyDYeeCS.C"
	echo -e "\n\tElectron code\n"
    elif [ "$1" == "mm" ] ; then
	code="studyDYmmCS.C"
	echo -e "\n\tMuon code\n"
    fi
else
    echo -e "\nUse:./run_studyDYmmCS.sh {ee|mm}"
    exit 1
fi

for v in ${vars} ; do
    echo "code=${code}, v=${v}"
    root -l -q -b ${code}+"(${v},${nSample},${doSave})"
done
