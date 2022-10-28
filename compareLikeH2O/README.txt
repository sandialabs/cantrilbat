
Directory to create comparison programs for aqueous species 
which only differ by the addition of H2O to the speciation.


Examples to use


BO2-  + 2H2O(L)  = B(OH)4-

H+ + H2O(L)   =  H3O+

ZnO2--  + 2 H2O(L)  = Zn(OH)4--   

ZnHO2-  +   H2O(L)  = Zn(OH)3-

ZnO(aq) +   H2O(L)  = Zn(OH)2(aq)

The principle is that these molecules are similar. We would like to create
a capability to create equivalent species with Zuzax that at least have
the same deltaG at 25 C and one bar for the above reactions.

We will have to modify the HKFT parameters so that they have a good chance of 
being close to having deltaG being 0 as a function of temperature as well.

The resulting volume change of the reaction is up in the air. There probably
does have a significant volume change, whethere there are covalent bonding
in the molecular or whether there is an inner shell hydrogen bonding interaction.
I don't know atm how to modify this. However, since BO2- for example was
calculated assuming the above formula and fitted to volume data using that formula,
I will go with the deltaVolume equaling zero for all of the reactions as well,
with an input to change the deltaVolume amount by a given value for each reaction,
that is a use input.

Plan
--------------------
The first thing to do is to stude the H+ reaction, making sure that deltaG=0 can be 
closely modeled with the HKFT parameters available.
Then, use those parameters to create the B(OH)4- HKFT parameters, and generalize until
the utility program is done.

Then, I should also invest some time in reviewing the literature to come up with an
literature review on what has been said about this issue.




