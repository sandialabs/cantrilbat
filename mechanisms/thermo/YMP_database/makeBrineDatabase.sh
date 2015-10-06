#!/bin/sh
#
#  Copywrite 2004 Sandia Corporation. Under the terms of Contract
#  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
#  retains certain rights in this software.
#  See file License.txt for licensing information.
#
# 
./makeActCoeff.sh
./makeSpeciesDatabase.sh

echo >HMW_BrineDatabase.xml

cat >HMW_BrineDatabase.xml  <<+ 
<?xml version="1.0"?>
<!--
     Copywrite 2004 Sandia Corporation. Under the terms of Contract
     DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
     retains certain rights in this software.
     See file License.txt for licensing information.
  -->
<ctml>

  <phase id="NaCl_electrolyte" dim="3">
    <speciesArray datasrc="#species_waterSolution">
               H2O(L) Cl- H+ Na+ OH- K+
    </speciesArray>
    <state>
      <temperature units="K"> 298.15 </temperature>
      <pressure units="Pa"> 101325.0 </pressure>
      <soluteMolalities>
             Na+:6.0954
             Cl-:6.0954
             H+:2.1628E-9
             OH-:1.3977E-6
      </soluteMolalities>
    </state> 
    <elementArray datasrc="elements.xml"> O H C K Si N Na Cl E </elementArray>
    <kinetics model="none" >
    </kinetics>
    <thermo model="HMW">
+

cat actCoeffSection.xml >> HMW_BrineDatabase.xml


cat >>HMW_BrineDatabase.xml <<+

      <solvent> H2O(L) </solvent>
    </thermo>
  </phase>

+

cat speciesDatabase.xml >> HMW_BrineDatabase.xml

cat >>HMW_BrineDatabase.xml <<+

</ctml>
+

