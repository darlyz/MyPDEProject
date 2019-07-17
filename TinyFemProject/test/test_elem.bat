@echo off
echo "test elemcalc.c ..."
set elemexec=teste.exe
if exist %elemexec% (
    del %elemexec%
)
gcc -I..\ test_elem.c ..\elemcalc.c -o teste.exe
teste.exe