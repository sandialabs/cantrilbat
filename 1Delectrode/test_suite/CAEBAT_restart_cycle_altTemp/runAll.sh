#!/bin/sh
#
function cleanup {
  if test -e runInProgress.txt
  then
    /bin/rm -f runInProgress.txt
    echo "  Run Interrupted: " `date` > runBad
    /bin/rm -f runGood
  fi
  exit 0
}

#
#  trap interrupt and signals. 
#
trap cleanup INT TERM TSTP QUIT

if test ! -d workdir1
then
   mkdir workdir1
fi
if test ! -d workdir2
then
   mkdir workdir2
fi
if test ! -d workdir3
then
   mkdir workdir3
fi

run
cp solutionStartEnd.xml  solutionStartEnd_begin.xml

cp  solutionStartEnd_begin.xml solutionRestart.xml

runOnce.sh   workdir1  0.1


cp  solutionStartEnd1.xml solutionRestart.xml

runOnce.sh   workdir2  0.2


cp  solutionStartEnd2.xml solutionRestart.xml
runOnce.sh   workdir3  0.3

cp  solutionStartEnd3.xml solutionRestart.xml


