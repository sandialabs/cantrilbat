#!/bin/sh
#
# Copywrite 2004 Sandia Corporation. Under the terms of Contract
# DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
# retains certain rights in this software.
# See file License.txt for licensing information.
#
if test ! -e cantera_SPEQ06_HKFT_params.xml -o -h cantera_SPEQ06_HKFT_params.xml
then
  /bin/rm -f cantera_SPEQ06_HKFT_params.xml
  ln -s cantera_SPEQ06_HKFT_params_rev6.xml cantera_SPEQ06_HKFT_params.xml
fi


echo >speciesDatabase.xml


cat  >speciesDatabase.xml <<+
<speciesData id="species_waterSolution">
    <!--
       Copywrite 2004 Sandia Corporation. Under the terms of Contract
       DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
       retains certain rights in this software.
       See file License.txt for licensing information.
     -->
    <species name="H2O(L)">
      <!-- H2O(L) liquid standard state -> pure H2O
           The origin of the NASA polynomial is a bit murky. It does
           fit the vapor pressure curve at 298K adequately.
        -->
      <atomArray>H:2 O:1 </atomArray>
      <thermo>
        <NASA Tmax="600.0" Tmin="273.14999999999998" P0="100000.0">
           <floatArray name="coeffs" size="7">
             7.255750050E+01,  -6.624454020E-01,   2.561987460E-03,  -4.365919230E-06,
             2.781789810E-09,  -4.188654990E+04,  -2.882801370E+02
           </floatArray>
        </NASA>
      </thermo>
      <standardState model="waterPDSS">
         <!--
              Molar volume in m3 kmol-1.
              (this is from Pitzer, Peiper, and Busey. However,
               the result can be easily derived from ~ 1gm/cm**3)
           -->
         <molarVolume> 0.018068 </molarVolume>
      </standardState>
    </species>

+

cat cantera_SPEQ06_HKFT_params.xml | cat >>speciesDatabase.xml


cat  >>speciesDatabase.xml <<+

</speciesData>
+


