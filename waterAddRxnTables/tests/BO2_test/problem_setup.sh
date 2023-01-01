#!/bin
#
#
TEST_NAME=UO2_pH_Speciation
#
#  Name of the program to run
#
PROGRAM=UO2_pH_Speciation
#
#  Any other program options
#
PROGRAM_OPTS=' '

#
#  
#
#   Text files to compare against
#
BLESSED_DATA_FILES=" output_0_blessed.txt  "
DATA_FILES="         out.txt       "
DIFF_NAMES="         diff_out.txt  "
DIFF_REQ="           True          "
#
#  CSV file to compare against
#
BLESSED_CSV_FILES=" molality_vs_pH_blessed.csv   "
CSV_FILES="         molality_vs_pH.csv      "
DIFF_CSV_NAMES="    diff_csv.txt   "
DIFF_CSV_REQ="      True "
#
#  Extra files to be removed before the test starts
#
EXTRA_WHACKED_FILES=''
