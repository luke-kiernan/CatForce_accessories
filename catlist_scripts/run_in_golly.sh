#!/bin/bash
# Goal: run a golly script with arguments, or read them in from golly interface
# APG's idea: define a variable pseudo_argv at the top of the script
#             then in the script, we test if said variables are defined
#             via try-except.
#                   NameError => we're running in Golly, so we ask user for input
#                   no NameError => use pseudo_argv values, put there by this script
# so now ./run_in_golly.sh my_golly_script.py arg1 arg2 ... argN does the desired thing.
GOLLYPATH="/Applications/golly-4.1-mac/Golly.app/contents/MacOS/Golly"
echo "pseudo_argv = r'''" > tmpfile.py
for var in "$@"; do
    echo "$var" >> tmpfile.py
done
echo "'''.split(chr(10))[1:-1]" >> tmpfile.py
cat "$1" >> tmpfile.py
$GOLLYPATH tmpfile.py