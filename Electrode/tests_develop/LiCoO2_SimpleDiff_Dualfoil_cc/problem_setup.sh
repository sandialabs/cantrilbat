#!/bin
#
#     Problem setup for runtest_Electrode
#
TEST_NAME=LiCoO2_SimpleDiff_Dualfoil_cc
#
#  Name of the program to run
#
PROGRAM=LiCoO2_Cathode_Radial_cc
#
#  Any other program options
#
PROGRAM_OPTS=' cathodeCC.inp  '
PROGRAM_STD_OUTPUT="ccout.txt"
#
#  
#
#   Text files to compare against
#
BLESSED_DATA_FILES="   good_ccout.txt"
DATA_FILES="                   ccout.txt "
DIFF_NAMES="               diff_ccout.txt"
DIFF_REQ="                        True "
#
#  CSV file to compare against
#
BLESSED_CSV_FILES=" LiCoO2_CC_intResults_blessed_0_0.csv   "
CSV_FILES="          LiCoO2_CC_intResults_0_0.csv  "
DIFF_CSV_NAMES="    diff_csv.txt    "
DIFF_CSV_REQ="      True              "
#
#  Extra files to be removed before the test starts
#
EXTRA_WHACKED_FILES=''
