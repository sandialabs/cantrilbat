#!/bin/sh
#
#  Example of how to 
#
#
#
if [ $# = 0 ] 
then
  exit 0
fi
list=$*
cdir=`pwd`
# echo my name  $0
ddd=`dirname $0`
myabs=$0
if test "x$ddd" = "x."
then
  myabs=$cdir/$0
fi
# echo my abs name = $myabs
# echo 'working directory for process, ' $$ ', is ' `pwd`
for ilist in $list
do
  if test $ilist = anode.inp
  then
    echo `pwd`"/"$ilist
    sed "s/Electrode Diameter/Electrode Gross Diameter/ 
         s/Electrode Thickness/Electrode Gross Thickness/
         s/Electrode Area/Electrode Gross Area/" anode.inp > ased.txt
    mv ased.txt anode.inp
  fi
  if test $ilist = cathode.inp
  then
    echo `pwd`"/"$ilist
    sed "s/Electrode Diameter/Electrode Gross Diameter/ 
         s/Electrode Thickness/Electrode Gross Thickness/
         s/Electrode Area/Electrode Gross Area/" cathode.inp > csed.txt
    mv csed.txt cathode.inp

  fi
  
  if [ -d $ilist ]
  then
     cd $ilist
     $myabs `/bin/ls -A`
     cd ..
  fi
done
exit 0

