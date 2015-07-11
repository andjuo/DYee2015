Adapted from Kevin Sung's code
https://github.com/ksung25/UserCode
by A.Juodagalvis on July 11, 2015

The code to produce the MIT Bambu files were found in
the repository by Christoph Paus
https://github.com/cpausmit?tab=repositories


Ilya Kravchenko wrote an explanation about the code (in April 17, 2013):

I have compiled info on how the present DY->ee ntuples are created,
and documenting it in this email. What I understand is that you will
take a look at this and see how easy/difficult it is for you to put
together electron ntupling in a common format with the muon
ntupling. If it is a matter of 2-3 days, it may be a good idea to do
it, but it if is 3+ weeks, maybe it is not a good idea.

There are two aspects to this info: running the code so that all
sequences are included and all corrections are done, and the way the
AOD info is accessed.

First of all, there are two stages to ntuple making.
- stage 1: a full dataset (such as DYToEE_M-20, or DoubleElectron)
is processed AOD(SIM)->Bambu. The "Bambu" format is
analogous to PAT-tuples, it is rather complete, but reduced,
event record. This is a CMSSW job.
- stage 2: Bambu root files are processed by ntupling code that
produces specifically ntuples for DY->ee analysis that contain
(almost) only variables of interest.

I imagine if you want to see how everything is implemented, you will
need to
1) start with the final ntuple's data structures
2) see how they are filled in "stage 2" from Bambu files, with which accessors
3) go to the code of "stage 1" and see how Bambu data structures are filled
from AOD(SIM) data structures.

So let's start with the final data structures. The whole package for final
ntupling ("stage 2") is found here:
http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/ksung/DYAna/

The data structures are in DYAna/EWKAna/Ntupler/interface/.
The ntuples we use contain events that look like this:
TTree Events
--> Branch "Info" - single object of TEventInfo.hh type
--> Branch "Gen" - single object of TGenInfo.hh type
--> Branch "Electron" - TClonesArray of TElectron.hh objects
--> Branch "Dielectron" - TClonesArray of TDielectron.hh objects
--> Branch "Muon" - TClonesArray of TMuon.hh objects
--> Branch "Photon" - TClonesArray of TPhoton.hh objects
--> Branch "Jet" - TClonesArray of TJet.hh objects
--> Branch "PV" - TClonesArray of TVertex.hh objects
Additionally, interface/EWKAnaDefs.hh contains definitions of various
constants like trigger bits.

Each of these data structures has many fields. We do not use
explicitly the "Jet", but we use all other branches. The branch
"Dielectron" contains full info about two electrons, so the fields in
Dielectron and Electron branches are filled by related values (there
is some duplication, for convenience).  It would take me non-trivial
time to go over all fields and email you the info what is used and
what is not in our present analysis, some number of hours… I am
guessing that at least 60% of all variables are used. In principle I
could do the inventory, but I would need to find time for it.

The code that fills these ntuples is from this DYAna/EWKAna/Ntupler
area, src/NtuplerModes.cc and interface/NtuplerModes.hh. This module
is run as a simple ROOT script (not CMSSW) by macros/runNtupler.C.
It takes as input the Bambu format files.

The code that creates Bambu root files out of AOD(SIM) is found in several
packages in UserCode:

MitAna <-- defines needed classes
MitCommon <-- defines many useful functions and computational tools
MitPhysics <-- contains support material such as BDT weight files,
jet energy corrections, etc. Much of it isn't used by the DYAna code, but
there are most likely some dependencies.

Additionally, the code that actually does AOD->Bambu production:
MitEdm
MitProd <-- most of the CMSSW python configs and sequences can be found here

Now, regarding running CMSSW with all the right sequences and corrections:
the master python config file is this one
MitProd/Configuration/python/BAMBUProd_AODSIM.py

Among the included fragments and configs this one is important:
MitProd/BAMBUSequences/python/BambuFillAOD_cfi.py

If you would like to run side by side your ntupling and our present one,
I can, e.g., share our present 8 TeV ntuples in which you can make a cut
on runs and lumis for data, or event numbers for MC, and just dump data
structures.

Let me know how it goes. If you decide that you can replicate/improve
ntupling for electrons, I would be interested to discuss the ntuple
structure with you before you do most of the work. (I find the present
final electron ntuple structure quite convenient, by the way.)

Ilya



// ------------------------------------------------------------------
// Part II (April 18, 2013)
// ------------------------------------------------------------------


here
is the set of instructions. You need to work on a computer that
has CMSSW available, and you will need the original AOD or AODSIM
example file for your chosen run. (I have done it at UNL Tier2 from
start to end just now). It is fairly straightforward, if you follow the
instructions below exactly.

(I am also cc'ing this to Andrius and Kevin for their records. Thanks
to Kevin and Guillelmo for babysitting me through all steps).

#----------------------------------------------
# Get all code and build it:
#----------------------------------------------

scram project CMSSW CMSSW_5_3_3_patch3
cd CMSSW_5_3_3_patch3/
cd src
cmsenv
export MIT_TAG=Mit_029a
export MIT_VERS=029a
# The setup step below is necessary as it checks out a number
# of needed packages with specific versions.
./MitAna/scripts/setup.sh
cvs co -d MitCommon -r $MIT_TAG UserCode/MitCommon
cvs co -d MitEdm -r $MIT_TAG UserCode/MitEdm
cvs co -d MitAna -r $MIT_TAG UserCode/MitAna
cvs co -d MitProd -r $MIT_TAG UserCode/MitProd
# Note: MitPhysics is not necessary for Stage 1 below,
# but may be needed via some inclusions for Stage 2.
cvs co -d MitPhysics -r $MIT_TAG UserCode/MitPhysics
cvs co -d EWKAna UserCode/ksung/DYAna/EWKAna
# build everything with two parallel processes. This takes O(1 hour)
scram build -j 2

#------------------------------------------
# Stage 1: AOD(SIM) -> Bambu
#----------------------------------------------

# Assuming you are continuing the above and are in $CMSSW_BASE/src,
# and cmsenv is executed

… Edit the top level config script, change there the number
… of events and the input file name (in format 'file:localname.root').
… Depending on whether you run on data or MC, the script is:
… MitProd/Configuration/python/BAMBUProd_AOD.py
… or
… MitProd/Configuration/python/BAMBUProd_AODSIM.py

# run production of Bambu root file
cmsRun ./MitProd/Configuration/python/BAMBUProd_AOD.py

# locate the output of the above, if no changes are made the default is
# something like XX-MITDATASET-XX_000.root in the current directory.

#------------------------------------------
# Stage 2: Bambu -> DY ntuples
#----------------------------------------------

# Assuming you are continuing the above and are in $CMSSW_BASE/src,
# and cmsenv is executed

cd EWKAna/Ntupler/macros

… Find there the file runNtupler.C
… In this file, find THE SECOND instance of runNtupler(...)
… function, designed for running on a single file. Change
… the default inputs there in the function declaration (data/mc, 2011 or not,
… input file name, etc). The input file name is, obviously, the file from Stage 1 above.

# run ntuple making as
root -b -q runNtupler.C+

# by default, the result is found in the ntuple.root

-------------------------------------------------------------------------------

Ilya

