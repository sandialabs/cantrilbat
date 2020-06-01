#!/bin
#
#     Problem setup for runtest_Electrode
#
TEST_NAME=CO2_solubility
TESTNAME=CO2_solubility
#
#  Name of the program to run
#
PROGRAM=CO2_solubility
#
#  Any other program options
#
PROGRAM_OPTS=' '
#
#  
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
BLESSED_CSV_FILES=" vcs_equilibrate_blessed.csv  vcs_equilibrate_blessed_1.csv vcs_equilibrate_blessed_2.csv vcs_equilibrate_blessed_3.csv vcs_equilibrate_blessed_4.csv  "
CSV_FILES="         vcs_equilibrate_res.csv      vcs_equilibrate_res_1.csv     vcs_equilibrate_res_2.csv     vcs_equilibrate_res_3.csv     vcs_equilibrate_res_4.csv "
DIFF_CSV_NAMES="    diff_csv.txt                 diff_csv_1.txt                diff_csv_2.txt                diff_csv_3.txt                diff_csv_4.txt                 "
DIFF_CSV_REQ="      True                         True                          True                          True                          True        "
#
#  Extra files to be removed before the test starts
#
EXTRA_WHACKED_FILES=''
