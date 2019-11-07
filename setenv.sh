#!/bin/bash

alias cmake='/cvmfs/sft.cern.ch/lcg/contrib/CMake/3.13.4/Linux-x86_64/bin/cmake'

# setup LCG env
#export lcgenv=/cvmfs/sft.cern.ch/lcg/releases/lcgenv/1.3.6-e63db/x86_64-centos7-gcc8-opt/lcgenv

# setup VDT
. /cvmfs/sft.cern.ch/lcg/releases/vdt/0.4.3-992df/x86_64-centos7-gcc8-opt/vdt-env.sh

# setup GCC
#. /cvmfs/sft.cern.ch/lcg/contrib/gcc/6.2/x86_64-centos7/setup.sh
#. /cvmfs/sft.cern.ch/lcg/contrib/gcc/7/x86_64-centos7/setup.sh
. /cvmfs/sft.cern.ch/lcg/contrib/gcc/8/x86_64-centos7/setup.sh

# setup PCRE
#. /cvmfs/sft.cern.ch/lcg/releases/pcre/8.43-511cb/x86_64-centos7-gcc8-opt/pcre-env.sh

#setup TBB
. /cvmfs/sft.cern.ch/lcg/releases/tbb/2019_U7-ba3eb/x86_64-centos7-gcc8-opt/tbb-env.sh

# setup Davix
. /cvmfs/sft.cern.ch/lcg/releases/Davix/0.7.3-d94fa/x86_64-centos7-gcc8-opt/Davix-env.sh

# setup png
. /cvmfs/sft.cern.ch/lcg/releases/png/1.6.37-9c2fe/x86_64-centos7-gcc8-opt/png-env.sh

#setup GSL
. /cvmfs/sft.cern.ch/lcg/releases/GSL/2.5-32fc5/x86_64-centos7-gcc8-opt/GSL-env.sh

# setup ROOT
#. /cvmfs/sft.cern.ch/lcg/releases/LCG_89/ROOT/6.10.02/x86_64-centos7-gcc7-opt/bin/thisroot.sh
#. /cvmfs/sft.cern.ch/lcg/releases/LCG_91/ROOT/6.10.06/x86_64-centos7-gcc7-opt/bin/thisroot.sh
. /cvmfs/sft.cern.ch/lcg/releases/LCG_96/ROOT/6.18.00/x86_64-centos7-gcc8-opt/bin/thisroot.sh
