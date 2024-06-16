#!/bin/sh

function help ()
{

    echo "Tests qi with a set of examples taken from \"quadrics.txt\""
    echo "in the range specified in arguments."
    echo
    echo "Syntax is: $0 [from] [to] [optional qi flags]"
    echo
    echo "Ex: To test examples from 10 to 40, do: $0 10 40"
    echo
    exit


}

if test $# -lt 2; then help; fi
if test ! -f ./quadrics.txt; then echo "File: \"quadrics.txt\" not found."; exit; fi
if test ! -f ./pick.sed; then echo "Sed file: \"pick.sed\" not found."; exit; fi
if test ! -f ./pick.sh; then echo "Script file: \"pick.sh\" not found."; exit; fi
if test ! -f ./qi; then echo "qi program not found."; exit; fi

from=$1
to=$2
test_id=$from

# Allows to pass every parameters from the third one
# to the qi program
shift
shift

while test $test_id -le $to
do
#echo "Test number $test_id" > /dev/stderr
  cat quadrics.txt | sed -f pick.sed | ./pick.sh $test_id | ./qi $*
  test_id=$(($test_id+1))
done

# eof

