#/bin/sh
#
# Multitask format for Makefiles
#

.PHONY : all
all:
	$(MAKE) alllink

.PHONY : multiPhase_test
multiPhase_test:
	cd multiPhase  &&                          $(MAKE) $(TASK)

.PHONY : she_electrode_test
she_electrode_test:
	cd she_electrode &&                       $(MAKE) $(TASK)

.PHONY : sce_electrode_test
sce_electrode_test:
	cd she_electrode  &&                      $(MAKE) $(TASK)

.PHONY : ex_gri_test
ex_gri_test:
	cd ex_gri  &&                              $(MAKE) $(TASK)

.PHONY : CO2_solubility_test_t
CO2_solubility_test_t:
	cd CO2_solubility_test &&                 $(MAKE) $(TASK)

TEST_LIST = multiPhase_test she_electrode_test sce_electrode_test ex_gri_test CO2_solubility_test_t

#
# TASK Rules
#
	
alllink: TASK = all
alllink: $(TEST_LIST)

test:  TASK = test
test:  $(TEST_LIST)

testd: TASK = testd
testd: $(TEST_LIST)

clean: TASK = clean
clean: $(TEST_LIST)

rebaseline: TASK = rebaseline
rebaseline: $(TEST_LIST)

depends: TASK = depends
depends: $(TEST_LIST)

clean_test: TASK = clean_test
clean_test: $(TEST_LIST)

clean_dist: TASK = clean_dist
clean_dist: $(TEST_LIST)

clean_test_flags: TASK = clean_test_flags
clean_test_flags: $(TEST_LIST)


