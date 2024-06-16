#!/bin/bash

function help () {
    
    echo "Syntax is: $0 [example_number]  <  [example_file]"
    echo
    echo "Expected format of \"example_file\" is a list of peers of"
    echo "matrices, such as:"
    echo "[1 0 0 -1 -1 0 2 32 -50 2]"
    echo "[2 6 3 0 0 0 -1 -3 -2 0]"
    echo "etc."
    echo
    echo "Produced output is a script file, understandable by the standalone"
    echo "\"qi\" program."
    echo
    
    if test -f "./pick.sed"
	then
	echo "Use the \"pick.sed\" sed file to convert the original list of quadrics into a"
	echo "simplified list of matrices, understandable by $0."
	echo
    fi

    # Exit error
    exit 1

}

if test $# -eq 0
then
    help
fi

# $1 gives the example number
# We can deduce the position of
# the pair of quadrics
second=$((2*$1))
first=$(($second-1))

sed -n "${first},${second}p" > /tmp/.pick_sh_matrices
echo "quadric_1=$(head -1 < /tmp/.pick_sh_matrices)"
echo "quadric_2=$(tail -1 < /tmp/.pick_sh_matrices)"
echo "intersect quadric_1 quadric_2"

# Exit success
exit 0

# eof
