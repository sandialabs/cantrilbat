from Cantera import *

###################################################################
#  vcs_Cantera   Command File
#
#   Master input file describing all of the options in the file
##################################################################
gas = importPhase('gas.xml')
C_Solid = importPhase('C_Solid.xml')
Fe_Solid = importPhase('Fe_Solid.xml')
CaO_Solid = importPhase('CaO_Solid.xml')
Fe3O4_Solid = importPhase('Fe3O4_Solid.xml')
CaCO3_Solid = importPhase('CaCO3_Solid.xml')
FeO_Solid = importPhase('FeO_Solid.xml')

Temperature = 1050.0
Pressure = 1.0 * OneAtm

mix = Mixture([gas, C_Solid, Fe_Solid, CaO_Solid,
               Fe3O4_Solid, CaCO3_Solid, FeO_Solid])

mix.setMoles('CaO(S):7.562E-4, N2:187.1E-3, Fe_Solid:4.2827E-2, C(d):8.8294E-2, H:1.3766E-2, O2:4.76974E-2')

mix.set(T = Temperature, P = Pressure)
mix.equilibrate('TP')
print mix

#-------------------------------------------------------------------------
#
#-------------------------------------------------------------------------
#START BLOCK KMolar Element Abundances 
#   O = 9.6151e-2    
#   H = 1.3766e-2
#   C = 8.8294e-2
#   Fe = 4.2827e-2
#   Ca = 7.5620e-4
#   N = 374.2E-3
#END BLOCK KMolar Element Abundances
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
#START BLOCK Species Initial KMoles 
#   CaO(S) = 7.562E-4
#   N2    = 187.1E-3
#   Fe_Solid = 4.2827E-2
#   C(d)  = 8.8294E-2 
#   H = 1.3766E-2
#   O2 = 4.76974E-2
#END BLOCK Species Initial KMoles
#-------------------------------------------------------------------------
#
#
#-------------------------------------------------------------------------
#START BLOCK Species Initial KMoles Guess 
#  SiH4  = 0.01
#  H2    = 0.99
#END BLOCK Species Initial KMoles Guess
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------



