#!/bin
#
#     Problem setup for runtest_Electrode
#
TEST_NAME=MCMBAnode_LiIon_SimpleDiff_cV
#
#  Name of the program to run
#
PROGRAM=MCMBAnode_SimpleDiff
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
BLESSED_DATA_FILES=" output_test_blessed.txt  "
DATA_FILES="         out.txt       "
DIFF_NAMES="         diff_out.txt  "
DIFF_REQ="           True          "
#
#  CSV file to compare against
#
BLESSED_CSV_FILES=" MCMB_intResults_blessed_0_0.csv  MCMB_globalResults_blessed_0_0.csv   "
CSV_FILES="         MCMB_intResults_0_0.csv          MCMB_globalResults_0_0.csv "
DIFF_CSV_NAMES="     diff_csv.txt                    diff_gcsv.txt "
DIFF_CSV_REQ="       True                            True"
#
#  Extra files to be removed before the test starts
#
EXTRA_WHACKED_FILES=''
