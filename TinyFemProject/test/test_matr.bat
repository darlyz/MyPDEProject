@echo off
echo "test matrcalc.c ..."
set matrexec=testm.exe
if exist %matrexec% (
    del %matrexec%
)
gcc -I..\ test_matr.c ..\matrcalc.c ..\elemcalc.c -o testm.exe
testm.exe