#!/bin
#    runtest_Electrode input file:
#
#     Problem setup for runtest_Electrode
#
TEST_NAME=FeS2_4_MPNewman_cc
#
#  Name of the program to run
#
PROGRAM=FeS2_4_MPNewman_cc
#
#  Any other program options
#
PROGRAM_OPTS=''
#
#  
SOLUTION_BLESSED_XML="   soln_0_0_blessed.xml "
SOLUTION_XML="           soln_0_0.xml         "
SOLUTION_XML_DIFFNAME="  diff_xml.txt         "
#
#   Text files to compare against
#
BLESSED_DATA_FILES=" good_out.txt  "
DATA_FILES="         out.txt       "
DIFF_NAMES="         diff_out.txt  "
DIFF_REQ="           True          "
#
#  CSV file to compare against
#
BLESSED_CSV_FILES="  FeS2_intResults_blessed.csv  outputTable_blessed.csv"
CSV_FILES="          FeS2_intResults_0_0.csv      outputTable.csv    "
DIFF_CSV_NAMES="     diff_csv.txt                 diff_ocsv.txt "
DIFF_CSV_REQ="       True                         True "
#
#  Extra files to be removed before the test starts
#
EXTRA_WHACKED_FILES=''
