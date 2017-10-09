#/bin/sh
###############################################################

.PHONY : all
all:
	$(MAKE) alllink

.PHONY : alllink test testd clean clean_test clean_dist clean_test_flags rebaseline install depends


LiSi_1_Newman_test: 
	cd LiSi_1_Newman &&                                         $(MAKE) $(TASK)

MCMB_OCVOverride_ConstTemp_test: 
	cd MCMB_OCVOverride_ConstTemp &&                       $(MAKE) $(TASK)

LiCoO2_SimpleDiff_Dualfoil_test: 
	cd LiCoO2_SimpleDiff_Dualfoil &&                       $(MAKE) $(TASK)

LiCoO2_SimpleDiff_Dualfoil_cc_test:
	cd LiCoO2_SimpleDiff_Dualfoil_cc &&                       $(MAKE) $(TASK)

MCMBAnode_LiIon_SimpleDiff_Dualfoil_test:
	cd MCMBAnode_LiIon_SimpleDiff_Dualfoil  &&                 $(MAKE) $(TASK)

MCMBAnode_LiIon_SimpleDiff_Dualfoil_cc_test:
	cd MCMBAnode_LiIon_SimpleDiff_Dualfoil_cc    &&          $(MAKE) $(TASK)

CSTR_ReadCompareSoln_test:
	cd CSTR_ReadCompareSoln                       &&          $(MAKE) $(TASK)

RadialDistrib_ReadCompareSoln_test:
	cd RadialDistrib_ReadCompareSoln             &&          $(MAKE) $(TASK)

LiCoO2_SimpleDiff_Jac_test:
	cd LiCoO2_SimpleDiff_Jac                     &&          $(MAKE) $(TASK)

TEST_LIST =                       \
         MCMB_OCVOverride_ConstTemp_test              \
         LiCoO2_SimpleDiff_Dualfoil_test              \
         LiCoO2_SimpleDiff_Dualfoil_cc_test           \
         MCMBAnode_LiIon_SimpleDiff_Dualfoil_test     \
         MCMBAnode_LiIon_SimpleDiff_Dualfoil_cc_test  \
         CSTR_ReadCompareSoln_test                    \
         RadialDistrib_ReadCompareSoln_test           \
         LiCoO2_SimpleDiff_Jac_test                   \


ALL_LIST = LiSi_1_Newman_test \
           $(TEST_LIST)

#
# TASK Rules
#

alllink: TASK = all
alllink: $(ALL_LIST)

test:  TASK = test
test:  $(TEST_LIST)

testd: TASK = testd
testd: $(TEST_LIST)

clean_test: TASK = clean_test
clean_test: $(ALL_LIST)

clean_dist: TASK = clean_dist
clean_dist: $(ALL_LIST)

clean_test_flags: TASK = clean_test_flags
clean_test_flags: $(TEST_LIST)

clean: TASK = clean
clean: $(ALL_LIST)

test_clean: TASK = clean_test
test_clean: $(TEST_LIST)

depends: TASK = depends
depends: $(ALL_LIST)

rebaseline: TASK = rebaseline
rebaseline: $(TEST_LIST)

