
Directory tests the OCV override capability for the MCMB anode.
It tests a specific temperature model.

In particular it doesn't test the Electrode object. It creates a reacting
surface domain and then tests the output from that.


Simulation of temperature rise in Li-ion cells at very high currents
Mao, Tiedmann, Newman, JPS (2014).


MCMB executable is a standalone program that prints out the experimental
OCV curve written up in Mao et al.

The MCMBThermo program prints out the thermo functions for the OCV curve
calculated from the straight thermodynamics functions from the XML files 


The MCMBThermo_RSD_Override program prints out the thermo functions from
the ReactingSurDomain object when the OCVoverride capability is used 
to caculate the functions.

In particular this executable prints out DeltaS deltaH and deltaG from the
reaction, and makes sure that it is valid.

MCMBThermo_RSD_Override then prints out the functions at 398 K as well as the 
base 298 K. In the past we have had problems when the temperature deviates
from 298 K. 

MCMBThermo_RSD_Override prints out the exchange current density formulation at
one point. It makes sure that the results are valid and agrees with the ROP
calcualtion carried out by ReactingSurDomain object.
It uses a ButlerVolmer_NoActivityCoeffs calculation for the ROP with an
Exchange_Current_Density formulation for the rate constant.


MCMBThermo_RSD_Override prints out a table of H S and G values for each of the
reactants in the reaction. It prints out the necessary delta values needed to 
change the XML derived delta rxn thermo values into the observed values. In
general the deltas needed are large. This should be cause us to be concerned
for the overall process of fitting and believing results.

MCMBThermo_RSD_Override uses an extended reaction description. The extended
reaction description includes an extra reaction representing the reaction 
mechanism's full cell reaction with respect to the reference electrode.
This involves adding the lithium metal phase to the phase description.
Then adding the production of Li(metal) reaction with a zero reaction rate.
We use the thermo functions from this reaction to calculate the OCV override
values, especially the temperature derivative of the OCV.
In contrast, we have attempted to do this in the past by ignoring the ln (x_i
gamma_i) term for the Li+ species.  However, I think this won't work
completely when we go to a full entropy formulation for the electrolyte phase,
which is needed. This means adding a new reaction to all mechanisms only for
its thermo values.

















HKM   
