

Test
       This routine has added phases in all three blocks.








------------------------ Solver Techniques  -------------------------------------------------------

Turned on Heat Source tracking and Resistance Tracking.
Root finder capability is turned off



----------------------- Geometries -----------------------------------------------------------------

Separator material = MgO.xml -> This is not true, but I'm banking on it being unimportant.
Cross sectional areas are set to 1 cm2

Anode particle diameter = 16 microns
Anode porosity = 0.40
Anode Thickness = 96 microns
Anode initial stoichiometric parameter set to 0.8 = Li_C6-bulk
Anode mole volume changed to 0.03431714 m3/kmole within the MCMB_Redlickkister.xml file. It was worked out 
   in the dualfoil comparison notes that this created a 1800 kg/m3 density for MCMB.


Separator Thickness = 25 microns

Electroyte concentrations are set to mole fractions determined to produce 1000 gmol m-3 salt concentration.
The analysis comes from the dualfoil comparision paper.


Initial mole fraction of LiCoO2 changed to 0.6




-------------------------------- Current and Rxn specification ---------------------------------------------


Replaced MCMBAnode_RedlichKister.xml file with the new one created in cantrilbat/batteryDevelop/MCMBAnode/MCMB_Mao.
This new one has a proper order of magnitude treatment of deltaH for the anode reaction. I noted that I have not
as of yet implemented the exact Mao model for the OCV and deltaH within the Electrode object. I just didn't follow
through on this.

ANODE
-----------------------

Anode Override block changed to be turned off. This is not current with the Mao paper, until I add it into the
Electrode object.

Anode Diffusion coefficient = 7.0 x 10-14 m2 s-1
Set the diffusion model to 1, which will then leave out the activity coefficient in the calculation of the diffusion.

Set the exchange current density formulation to a BV_ButlerVolmer_NoActivityCoeffs formulation. with
an exchange current density of

              i_o = 1.90343E7 exp [ -33.25 kJ/gmol /RT ] 

This i_o corresponds to 28.28 amps / m2 at 298K. The calcualtions behind this number are described in the dualfoil
comparison notes. However, noticing that Mao et al. has changed the i_o, I get

    i_o = (96485.3) (3.0E-9) (1.282E4)**0.5 (2.498E4) = 818.37 amps /m2

So the following results:

           i_o = 5.5082E8 exp [ -33.25 kJ/gmol /RT ]

CATHODE
-----------------------

           i_0 = 8.7E8 exp [ -33.25 kJ/gmol /RT ]

yields a value of 200 A / m2 for a 0.1 X solution at starting conditions of y = 0.6.





---------------------- Heat transfer environment ----------------------------------------------------

Bath temperature = 298.15

Heat transfer coefficients on both sides = 0.184 W m-2 K-1 (0.368 W m-2 K-1 in total)







