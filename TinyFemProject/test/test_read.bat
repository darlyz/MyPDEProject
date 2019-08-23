@echo off
echo "test ReadMesh.c ..."
set exec=testr.exe
if exist %exec% (
    del %exec%
)
gcc -I..\ test_read.c ..\ReadMesh.c -o testr.exe
testr.exe ..\mesh\exam.gid\exam.dat -init