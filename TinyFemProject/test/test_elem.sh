#!/bin/sh
echo "test elemcalc.c ..."
if [ -f "./teste.out" ]
then
rm -f teste.out
fi
gcc -I../ test_elem.c ../elemcalc.c -o teste.out
./teste.out
