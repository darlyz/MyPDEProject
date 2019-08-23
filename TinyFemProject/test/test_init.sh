#!/bin/sh
echo "test eleminit.c ..."
if [ -f "./testi.out" ]
then
rm -f testi.out
fi
gcc -I../ test_init.c ../ReadMesh.c ../initial.c -o testi.out
./testi.out
