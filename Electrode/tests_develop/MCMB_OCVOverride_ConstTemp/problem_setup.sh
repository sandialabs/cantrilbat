#!/bin
#
#     Problem setup for runtest_Electrode
#
TEST_NAME=MCMB_OCVOverride_ConstTemp
#
#  Name of the program to run
#
PROGRAM=MCMBThermo
#
#  Any other program options
#
PROGRAM_OPTS='  '
PROGRAM_STD_OUTPUT="MCMBThermo.out"
#
#  
#
#   Text files to compare against
#
BLESSED_DATA_FILES=" MCMB_blessed.out  MCMBThermo_blessed.out "
DATA_FILES="         MCMB.out          MCMBThermo.out "
DIFF_NAMES="         diff_out.txt      diff_tout.txt"
DIFF_REQ="           True              True "
#
#  CSV file to compare against
#
BLESSED_CSV_FILES=" thermo_blessed.csv   "
CSV_FILES="         thermo.csv "
DIFF_CSV_NAMES="    diff_csv.txt    "
DIFF_CSV_REQ="      True              "
#
#  Extra files to be removed before the test starts
#
EXTRA_WHACKED_FILES=''
