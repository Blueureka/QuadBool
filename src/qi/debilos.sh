#!/bin/bash
cat quadrics.txt | sed -f pick.sed | pick.sh $1 | ./qi   
#src/qi/qi -l 3 < src/qi/$1
