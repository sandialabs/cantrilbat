#!/bin
#
#     Problem setup for runtest_Electrode
#
TEST_NAME=MCMBAnode_SimpleDiff_cc_restart
#
#  Name of the program to run
#
PROGRAM=MCMBAnode_SimpleDiff_cc
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
BLESSED_CSV_FILES="  MCMB_globalResults_blessed_0_0.csv  outputTable_blessed.csv  "
CSV_FILES="          MCMB_globalResults_0_0.csv          outputTable.csv          "
DIFF_CSV_NAMES="     diff_gcsv.txt                       diff_ocsv.txt            "
DIFF_CSV_REQ="       True                                True                     "
#
#  Extra files to be removed before the test starts
#
EXTRA_WHACKED_FILES=''
