cfjovec[08-18-2012]:

     Copywrite 2004 Sandia Corporation. Under the terms of Contract
     DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
     retains certain rights in this software.
     See file License.txt for licensing information.

Brief Description of YMP Thermodynamic Databases:

The thermodynamic parameter databases for electrolytes, solids, and aqueous species attached to this directory
are compilations from the scientific literature and handbook sources developed for the
Yucca Mountain Project (YMP) high-level nuclear waste repository program. Three main thermodynamic databases resulted 
from the YMP effort: Pitzer (data0.ypf.R2), dilute systems (data0.ymp.R5), and HKFT (SPEQ06).  The Pitzer 
thermodynamic database was developed to evaluate electrolyte solutions resulting from the compositional evolution 
of waters in evaporative environments. The database for dilute systems (data0.ymp.R5) comprises large sets of 
of aqueous species, solids, and gases mainly derived from the 'SPEQ06' database for the SUPCRT92 code plus 
other additions. HKFT is mainly a compilation and revisions to data developed for use with the SUPCRT92 code which 
also feeds into the aforementioned data0.ymp.R5 database.  The SPEQ06 database is built upon the work and thermodynamic 
data developed by late Prof. Helgeson and students at UC Berkeley.  These YMP datasets were developed as data inputs 
for calculations performed using the computer code EQ3/6 (Wolery and Jarek 2003). However, these have been 
reformatted for use with certain classes in the Cantera code suite.  The HKFT (or SPEQ06) database used in the 
YMP are for aqueous species and minerals to be used with the PDSS_HKFT and MineralEQ3 classes in Cantera. The HKFT
data blocks are based on the HKFT EoS for aqueous species consistent with the parameters developed for the
SUPCRT92 code. Similarly, thermodynamic data for gases and mineral solids with stoichiometric compositions are
consistent with SUPCRT92 input types for parameter extrapolation. It should be noted that Cantera accepts other types
of thermodynamic data inputs based on other database developments such as NASA and NIST.  
The data0.ypf.R2 files contain data blocks with Pitzer parameters to be used with the HMWSoln class.
It should be noted that additional databases have branched from the three main YMP databases and are not
included here.  For example, a modification of the Pitzer database resulted in the addition of Pitzer parameters
for ammonium species that are not included in the current compilation.  These and other additions / modifications
will eventually be incorporated in the current Cantera YMP databases. Some solids in the data0.ymp.R5 are not
included in the current Cantera YMP database compilation due mainly to incomplete thermodynamic parameter sets
that are required as Cantera inputs. This missing data is mainly restricted to actinide solids. 
The 'Reference_Sources_Codes.txt' file has all the reference sources for thermodynamic parameters in the SPEQ06_mins
file. Details of these database developments are documented in the following YMP reports:

Mariner, P. E., 2007, In-Drift Precipitates/Salts Model (ANL-EBS-MD-000045 REV 03).  Las Vegas, Nevada, 
Sandia National Laboratories; OCRWM Lead Laboratory for Repository Systems, Appendix I, 358 pp. 

Wolery, T. J. and Jové Colón, C. F., 2007, Qualification of Thermodynamic Data for Geochemical Modeling 
of Mineral–Water Interactions in Dilute Systems (ANL-WIS-GS-000003 REV 01): Las Vegas, Nevada, 
Sandia National Laboratories; OCRWM Lead Laboratory for Repository Systems, 412 pp.

EQ3/6 User Manual: Wolery TJ, Jarek RL. 2003. EQ3/6, Version 8.0: Software User's Manual. Albuquerque: 
Sandia National Laboratories, 376 pp.

These documents are available through through the U.S. Nuclear Regulatory Commission (NRC) web-based
ADAMS Public Documents (see http://www.nrc.gov). 

The following data blocks files are formatted for use in Cantera input data files:

cantera_SPEQ06_HKFT_params_rev5.xml
data0.ypf.R2_betas_ca_rev2.xml 
data0.ypf.R2_extract_Pitzer_params_lambdas_rev3.xml
data0.ypf.R2_extract_Pitzer_params_psi_rev3.xml
data0.ypf.R2_extract_Pitzer_params_zetas_rev2.xml
data0.ypf.R2_Thetas_cc_aa_rev2.xml
SPEQ06_mins_no-phase_trans_rev1.xml
Reference_Sources_Codes.txt
