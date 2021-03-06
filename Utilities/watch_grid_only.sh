#!/bin/bash

if [[ "$1" == "" ]]; then echo "Please pass a Submission directory!"; exit 1; fi;

cd ${1}

for f in $(ls .)
  do

  MASS=$f
  cd $MASS || exit 1
  ls crab_* > /dev/null 2>&1 || exit 0;
  
  for d in crab_*; do 
      crab -c $d -status >& /dev/null
      crab -c $d -status -get >& /dev/null
      crab -c $d -status 2>&1 | tee status.$d;
      
  done
  
  cd ../
  
done