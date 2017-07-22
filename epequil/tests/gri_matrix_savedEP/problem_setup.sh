#!/bin
#
#     Problem setup for runtest_Electrode
#
TEST_NAME=gri_matrix_savedEP
#
#  Name of the program to run
#
PROGRAM=gri_matrix_savedEP
#
#  Any other program options
#
PROGRAM_OPTS=' '
PROGRAM_STD_OUTPUT=out.txt
#
#  
#
#   Text files to compare against
#
BLESSED_DATA_FILES=" output_blessed.txt  "
DATA_FILES="         out.txt       "
DIFF_NAMES="         diff_out.txt  "
DIFF_REQ="           True          "
#
#  CSV file to compare against
#
BLESSED_CSV_FILES="   "
CSV_FILES="               "
DIFF_CSV_NAMES="      "
DIFF_CSV_REQ="               "
#
#  Extra files to be removed before the test starts
#
EXTRA_WHACKED_FILES=''
