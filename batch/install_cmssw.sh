source /cvmfs/cms.cern.ch/cmsset_default.sh
#export SCRAM_ARCH=slc7_amd64_gcc10
thisdir=$(pwd)
release=CMSSW_12_6_0_patch1
if [ $# -gt 0 ]
then
    if [ ${1} == "2022central" ]
    then
	release=CMSSW_12_4_16
    fi
    if [ ${1} == "2023" ]
    then
	release=CMSSW_13_0_10
    fi
    if [ ${1} == "2023central" ]
    then
	release=CMSSW_13_0_14
    fi
    if [ ${1} == "2024" ]
    then
    release=CMSSW_14_0_18
    fi
    if [ ${1} == "2024central" ]
    then
    release=CMSSW_15_0_2
    fi
    if [ ${1} == "2025" ]
    then
    release=CMSSW_15_0_15_patch1
    fi
fi

exists=1
if [ ! -d ${release} ]
then
    cmsrel $release
    exists=0
fi

cd $release/src
cmsenv
if [ ${exists} -lt 1 ]
then
    if [ ! -L Scouting ]
    then
	ln -s ../../Scouting/
    fi
    scram b -j12
fi
cd -
