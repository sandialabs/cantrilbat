

4/15/15

LiCoO2.cpp  captures the thermodynamics of Mao Tiedemann and Newmann (2014).

It shows DeltaH, DeltaS, and DeltaG for the LiCoO2 reaction. 

The basic result is that DeltaH is almost equal to DeltaG.

  DeltaH = 395 kJoules/ gmol

The blessed file is located in LiCoO2_blessed.out



With this basic result, I've modified the  LiCoO2_RedlichKister.xml file to qualitatively get the same order of magnitude result. for Delta H. It was qualitatively getting the same result for DeltaG. However, combined with 
an error in the RedlichKister object within Cantera, it was getting a significantly different result. 


The real issue still remains that the OCV and DeltaH measurements don't coincide with the thermo results.
The thermo results should be paired and fit to experimental data.







