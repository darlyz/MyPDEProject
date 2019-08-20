@echo off
echo "test matrcalc.c ..."
set matrexec=testm.exe
if exist %matrexec% (
    del %matrexec%
)
gcc -I..\ test_matr.c ..\matrcalc.c ..\elemcalc.c ..\gaussain.c ..\shap_func.c ..\calc_shap.c -o testm.exe
testm.exe