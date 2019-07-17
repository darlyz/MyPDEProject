#!/bin/sh
echo "test matrcalc.c ..."
if [ -f "./testm.out" ]
then
rm -f testm.out
fi
gcc -I../ test_matr.c ../matrcalc.c ../elemcalc.c -o testm.out
./testm.out
