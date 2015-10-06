#!/bin/sh
#
# Copywrite 2004 Sandia Corporation. Under the terms of Contract
# DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
# retains certain rights in this software.
# See file License.txt for licensing information.
#
#
#
#if test ! -f cantera_SPEQ06_HKFT_params.xml
#then
#  ln -s cantera_SPEQ06_HKFT_params_rev1.xml cantera_SPEQ06_HKFT_params.xml
#fi


echo >MineralsDatabase.xml


cat  >MineralsDatabase.xml <<+
<ctml>
  <speciesData id="species_MineralsEQ3">
    <!--
       Copywrite 2004 Sandia Corporation. Under the terms of Contract
       DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
       retains certain rights in this software.
       See file License.txt for licensing information.
     -->
+

cat SPEQ06_mins_no-phase_trans_rev2.xml | cat >>MineralsDatabase.xml

cat salt_mins_rev0.xml | cat >> MineralsDatabase.xml

cat  >> MineralsDatabase.xml <<+

  </speciesData>
</ctml>
+


