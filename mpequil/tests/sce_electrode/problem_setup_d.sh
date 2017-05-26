#!/bin
#
#     Problem setup for runtest_Electrode
#
TEST_NAME=mpequil_sce_electrode
#
#  Name of the program to run
#
PROGRAM=mpequil
#
#  Any other program options
#
PROGRAM_OPTS=' mpequil.inp -d 8'
#
PROGRAM_STD_OUTPUT='dout.txt'
#  
#
#   Text files to compare against
#
BLESSED_DATA_FILES=" good_dout.txt  "
DATA_FILES="         dout.txt       "
DIFF_NAMES="         diff_dout.txt  "
DIFF_REQ="           True          "
#
#  CSV file to compare against
#
BLESSED_CSV_FILES="  vcs_equilibrate_blessed.csv "
CSV_FILES="          vcs_equilibrate_res.csv        "
DIFF_CSV_NAMES="     diff_csv.txt    "
DIFF_CSV_REQ="       True            "
#
#  Extra files to be removed before the test starts
#
EXTRA_WHACKED_FILES=''
