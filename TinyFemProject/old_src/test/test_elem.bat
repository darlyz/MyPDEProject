@echo off
echo "test elemcalc.c ..."
set exec=teste.exe
if exist %exec% (
    del %exec%
)
gcc -I..\ test_elem.c ..\elemcalc.c ..\gaussain.c ..\shap_func.c ..\calc_shap.c -o teste.exe
teste.exe