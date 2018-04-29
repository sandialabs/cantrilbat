#!/bin
#
#     Problem setup for runtest_Electrode
#
WD=`pwd`
TEST_NAME=`basename $WD`
#
#  Name of the program to run
#
PROGRAM=LiCoO2_Cathode_3
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
DIFF_REQ="           False          "
#
#  CSV file to compare against
#

BLESSED_CSV_FILES="  LiCoO2_globalResults_0_0_blessed.csv   results_blessed.csv  "
CSV_FILES="          LiCoO2_globalResults_0_0.csv           results.csv  "
DIFF_CSV_NAMES="     diff_gcsv.txt                          diff_res.txt "
DIFF_CSV_REQ="       True                                   True "

#
#  Extra files to be removed before the test starts
#
EXTRA_WHACKED_FILES=''
