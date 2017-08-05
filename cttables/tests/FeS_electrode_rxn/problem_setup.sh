#!/bin
#
#     Problem setup for runtest_Electrode
#
TEST_NAME=cttables_FeS_electrode_rxn
#
#  Name of the program to run
#
PROGRAM=cttables
#
#  Any other program options
#
PROGRAM_OPTS=' '
#
#  
#
#   Text files to compare against
#
BLESSED_DATA_FILES=" good_out.txt "
DATA_FILES="         out.txt     "
DIFF_NAMES="         diff_out.txt  "
DIFF_REQ="           True          "
#
#  CSV file to compare against
#
BLESSED_CSV_FILES="  polarizationData_blessed.csv "
CSV_FILES="          polarizationData.csv     "
DIFF_CSV_NAMES="     diff_csv.txt "
DIFF_CSV_REQ="       True        "
#
#  Extra files to be removed before the test starts
#
EXTRA_WHACKED_FILES=''
