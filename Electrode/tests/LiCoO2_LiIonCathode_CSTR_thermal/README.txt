LiCoO2_LiIonCathode_CSTR_thermal

Date of Test            4/6/2014

Owner                   Harry Moffat


Description of the test
---------------------------------------------------------------------

     Electrode Object "CSTR_LiCoO2Cathode"

     The basic idea is to calculate the thermal terms in two ways
  and make sure that the output is the same.
   1-    electrodeC->thermalEnergySourceTerm_overpotential(0)
   2-    electrodeC->thermalEnergySourceTerm_Overpotential_SingleStep();

  Then make sure the source terms are internally consistent.
        DeltaH = DeltaG + deltaS

   Integrate at constant current for a couple of global steps.


Features that are turned on
----------------------------------------------------------------------------------------------

     Thermal Source terms

     There is detailed printing.

     There is an exchange current density specification for the reaction rate


Keywords:

     Thermal Sources
     CSTR_LiCoO2Cathode


History:


Status:
    Active and needed
    valgrind checked

