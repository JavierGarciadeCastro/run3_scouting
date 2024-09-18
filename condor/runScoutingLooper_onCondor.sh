#!/bin/bash

export X509_USER_PROXY=$(voms-proxy-info -path)

usage()
{
    echo "Usage:"
    echo ""
    echo "  sh condor/runScoutingLooper_onCondor.sh ['notar'] ['2023'] input_dir"
    echo ""
    echo "The output directory will be created in /ceph/cms/store/user/$USER/Run3ScoutingOutput/"
    echo "Control the jobs to be run by editing runScoutingLooper_onCondor.sub or the corresponding 2023 file"
    echo ""
    exit
}

if [ -z $1 ]; then usage; fi

notar=0
year23=0
indir=""

if [ $# -gt 1 ]
then
   if [ "$1" == "notar" ]
   then
      if [ "$2" == "2023" ]
      then
         notar=1
         year23=1
         indir=$3
      else
         notar=1
         indir=$2
      fi
   else
      if [ "$1" == "2023" ]
      then
         year23=1
         indir=$2
      else
         indir=$1
      fi
   fi
fi

export SCOUTINGOUTPUTDIR=${indir}
export STARTDIR=$PWD

mkdir -p condor/plotting_logs
mkdir -p /ceph/cms/store/user/$USER/Run3ScoutingOutput/$SCOUTINGOUTPUTDIR

if [ ${notar} -gt 0 ]
then
   if [ ${year23} -gt 0 ]
   then
      condor_submit condor/runScoutingLooper_onCondor2023.sub
   else
      condor_submit /home/users/garciaja/run3Scouting/run3_scouting/condor/runScoutingLooper_onCondor_BToPhi.sub
   fi
else
   sh condor/create_package.sh
   if [ ${year23} -gt 0 ]
   then
      condor_submit condor/runScoutingLooper_onCondor2023.sub
   else
      condor_submit /home/users/garciaja/run3Scouting/run3_scouting/condor/runScoutingLooper_onCondor_BToPhi.sub
   fi
fi
