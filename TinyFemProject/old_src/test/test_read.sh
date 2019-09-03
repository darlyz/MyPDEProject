#!/bin/sh
echo "test ReadMesh.c ..."
if [ -f "./testr.out" ]
then
rm -f testr.out
fi
gcc -I../ test_read.c ../ReadMesh.c -o testr.out
./testr.out ../mesh/exam.gid/exam.dat -mesh
