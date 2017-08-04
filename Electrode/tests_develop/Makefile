#/bin/sh
###############################################################

all:
	cd LiSi_1_Newman ;                                          $(MAKE)
	cd MCMB_OCVOverride_ConstTemp ;                             $(MAKE)
	cd LiCoO2_SimpleDiff_Dualfoil;				    $(MAKE)
	cd LiCoO2_SimpleDiff_Dualfoil_cc ;                          $(MAKE)
	cd MCMBAnode_LiIon_SimpleDiff_Dualfoil; 		    $(MAKE)
	cd MCMBAnode_LiIon_SimpleDiff_Dualfoil_cc; 		    $(MAKE)
	cd CSTR_ReadCompareSoln ;				    $(MAKE)
	cd RadialDistrib_ReadCompareSoln ;			    $(MAKE)
	cd LiCoO2_SimpleDiff_Jac ;				    $(MAKE) 
#	cd GFCEO_CSTR_LiCoO2_LiIonCathode_cV ;                      $(MAKE) all

test:
#	cd LiSi_1_Newman ;     (not ready yet)                      $(MAKE) -s test
	cd MCMB_OCVOverride_ConstTemp ;                             $(MAKE) -s test
	cd LiCoO2_SimpleDiff_Dualfoil ;                             $(MAKE) -s test
	cd LiCoO2_SimpleDiff_Dualfoil_cc ;                          $(MAKE) -s test
	cd MCMBAnode_LiIon_SimpleDiff_Dualfoil; 		    $(MAKE) -s test
	cd MCMBAnode_LiIon_SimpleDiff_Dualfoil_cc; 		    $(MAKE) -s test
	cd CSTR_ReadCompareSoln ;				    $(MAKE) -s test
	cd RadialDistrib_ReadCompareSoln ;			    $(MAKE) -s test
	cd LiCoO2_SimpleDiff_Jac ;				    $(MAKE) -s test
#	cd GFCEO_CSTR_LiCoO2_LiIonCathode_cV ;                      $(MAKE) -s test

rebaseline:
	cd MCMB_OCVOverride_ConstTemp ;                             $(MAKE) rebaseline
	cd LiCoO2_SimpleDiff_Dualfoil ;                             $(MAKE) rebaseline
	cd LiCoO2_SimpleDiff_Dualfoil_cc ;                          $(MAKE) rebaseline
	cd MCMBAnode_LiIon_SimpleDiff_Dualfoil; 		    $(MAKE) rebaseline
	cd MCMBAnode_LiIon_SimpleDiff_Dualfoil_cc; 		    $(MAKE) rebaseline
	cd CSTR_ReadCompareSoln ;				    $(MAKE) rebaseline
	cd RadialDistrib_ReadCompareSoln ;			    $(MAKE) rebaseline
	cd LiSi_1_Newman ;                                          $(MAKE) rebaseline
	cd LiCoO2_SimpleDiff_Jac ;				    $(MAKE) rebaseline
	cd GFCEO_CSTR_LiCoO2_LiIonCathode_cV ;                      $(MAKE) rebaseline

clean:
	cd LiSi_1_Newman ;                                          /bin/rm -f .depends; $(MAKE) clean
	cd MCMB_OCVOverride_ConstTemp ;                             $(MAKE) clean
	cd LiCoO2_SimpleDiff_Dualfoil;				    $(MAKE) clean
	cd LiCoO2_SimpleDiff_Dualfoil_cc ;                          $(MAKE) clean
	cd MCMBAnode_LiIon_SimpleDiff_Dualfoil; 		    $(MAKE) clean
	cd MCMBAnode_LiIon_SimpleDiff_Dualfoil_cc; 		    $(MAKE) clean
	cd CSTR_ReadCompareSoln ;				    $(MAKE) clean
	cd RadialDistrib_ReadCompareSoln ;			    $(MAKE) clean
	cd LiCoO2_SimpleDiff_Jac ;				    $(MAKE) clean
	cd GFCEO_CSTR_LiCoO2_LiIonCathode_cV ;                      $(MAKE) clean

clean_dist:
	cd LiSi_1_Newman ;                                          /bin/rm -f .depends; $(MAKE) clean
	cd MCMB_OCVOverride_ConstTemp ;                             $(MAKE) clean
	cd LiCoO2_SimpleDiff_Dualfoil;                              $(MAKE) clean
	cd LiCoO2_SimpleDiff_Dualfoil_cc ;                          $(MAKE) clean
	cd MCMBAnode_LiIon_SimpleDiff_Dualfoil;                     $(MAKE) clean
	cd MCMBAnode_LiIon_SimpleDiff_Dualfoil_cc; 		    $(MAKE) clean
	cd CSTR_ReadCompareSoln ;				    $(MAKE) clean
	cd RadialDistrib_ReadCompareSoln ;			    $(MAKE) clean
	cd LiCoO2_SimpleDiff_Jac ;				    $(MAKE) clean
	cd GFCEO_CSTR_LiCoO2_LiIonCathode_cV ;                      $(MAKE) clean_dist

depends:
	cd LiSi_1_Newman ;                                          /bin/rm -rf .depends; $(MAKE) depends
	cd MCMB_OCVOverride_ConstTemp ;                             $(MAKE) depends
	cd LiCoO2_SimpleDiff_Dualfoil;				    $(MAKE) depends
	cd LiCoO2_SimpleDiff_Dualfoil_cc ;                          $(MAKE) depends
	cd  MCMBAnode_LiIon_SimpleDiff_Dualfoil; 		    $(MAKE) depends
	cd MCMBAnode_LiIon_SimpleDiff_Dualfoil_cc; 		    $(MAKE) depends
	cd CSTR_ReadCompareSoln ;				    $(MAKE) depends
	cd RadialDistrib_ReadCompareSoln ;			    $(MAKE) depends
	cd LiCoO2_SimpleDiff_Jac ;				    $(MAKE) depends
	cd GFCEO_CSTR_LiCoO2_LiIonCathode_cV ;                      $(MAKE) depends


install:


