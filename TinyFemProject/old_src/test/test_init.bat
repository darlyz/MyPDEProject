@echo off
echo "test eleminit.c ..."
set exec=testi.exe
if exist %exec% (
    del %exec%
)
gcc -I..\ test_init.c ..\ReadMesh.c ..\initial.c -o testi.exe
testi.exe