@echo off
echo "test matrcalc.c ..."
set exec=testm.exe
if exist %exec% (
    del %exec%
)
gcc -I..\ test_matr.c ..\matrcalc.c ..\elemcalc.c ..\gaussain.c ..\shap_func.c ..\calc_shap.c -o testm.exe
testm.exe