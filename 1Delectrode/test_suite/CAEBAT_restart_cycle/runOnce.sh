#!/bin/sh


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


runDir=$1
endTime=$2

DATA="ECsoln_ion.xml  LiCoO2Cathode_electrode.xml  LiCoO2_Margules_1.xml  Li_Metal.xml \
MCMBAnode_electrode.xml  MCMB_RedlichKister.xml  metal_Li_LiIon_electrons.xml  MgO.xml run \
anode.inp cathode.inp"

restartFile=solutionRestart.xml

sss=`echo $runDir | sed s/workdir// `
# echo 'suffix = ' $sss

cd $runDir
/bin/rm -rf  $DATA
/bin/rm -rf LiIon_PorousBat


for i in $DATA 
do
   # echo "ln -s ../$i  $i"
   ln -s ../$i .
done

if test ! -x LiIon_PorousBat   ; then
  if test -x ../../../src_LiIon_PorousBat/LiIon_PorousBat ; then
    ln -s ../../../src_LiIon_PorousBat/LiIon_PorousBat .
  else
    echo 'ERROR:  LiIon_PorousBat  executable can not be found'
    exit -1
  fi
fi

cd ..

cp LiIon_PorousBat.inp.tmpl $runDir
cp $restartFile $runDir

cd $runDir

cat LiIon_PorousBat.inp.tmpl | sed s/XXXENDTIMEXXX/$endTime/ > LiIon_PorousBat.inp


../run

cp accumulatedHeatSource.txt ../accumulatedHeatSource${sss}.txt  
cp solutionStartEnd.xml  ../solutionStartEnd${sss}.xml 
cp accumulatedElectricalOutput.txt ../accumulatedElectricalOutput${sss}.txt  

cd ..


