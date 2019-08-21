@echo off
echo "test ReadMesh.c ..."
set meshexec=testr.exe
if exist %meshexec% (
    del %meshexec%
)
gcc -I..\ test_read.c ..\ReadMesh.c -o testr.exe
testr.exe ..\mesh\exam.gid\exam.dat -init