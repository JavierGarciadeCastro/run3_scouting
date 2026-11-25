thisdir=$(pwd)
release=CMSSW_12_6_0_patch1
arch=el8_amd64_gcc10
if [ $# -gt 0 ]
then
    if [ ${1} == "2022central" ]
    then
	release=CMSSW_12_4_16
    fi
    if [ ${1} == "2023" ]
    then
	release=CMSSW_13_0_10
	arch=el8_amd64_gcc11
    fi
    if [ ${1} == "2023central" ]
    then
	release=CMSSW_13_0_14
	arch=el8_amd64_gcc11
    fi
    if [ ${1} == "2024" ]
    then
    release=CMSSW_14_0_18
    arch=el8_amd64_gcc12
    fi
    if [ ${1} == "2024central" ]
    then
    release=CMSSW_15_0_2
    arch=el8_amd64_gcc12
    fi
    if [ ${1} == "2024private" ]
    then
    release=CMSSW_14_0_21
    arch=el8_amd64_gcc12
    fi
    if [ ${1} == "2025" ]
    then
    release=CMSSW_15_0_15_patch1
    arch=el8_amd64_gcc12
    fi
fi

export SCRAM_ARCH=${arch}
source /cvmfs/cms.cern.ch/cmsset_default.sh

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
