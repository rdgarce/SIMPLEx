@echo off

set incl_test=-I ../src -I ../src/lina
set src_test=dynamicTestPL.c SIMPLEx_test.c ../src/SIMPLEx.c ../src/lina/lina.c ../src/lina/qr.c

set incl_lplib=-I lpsolve -I lpsolve/bfp -I lpsolve/bfp/bfp_LUSOL -I lpsolve/bfp/bfp_LUSOL/LUSOL -I lpsolve/colamd -I lpsolve/shared
set src_lplib=lpsolve/lp_MDO.c lpsolve/shared/commonlib.c lpsolve/colamd/colamd.c lpsolve/shared/mmio.c lpsolve/shared/myblas.c lpsolve/ini.c lpsolve/lp_rlp.c lpsolve/lp_crash.c lpsolve/bfp/bfp_LUSOL/lp_LUSOL.c lpsolve/bfp/bfp_LUSOL/LUSOL/lusol.c lpsolve/lp_Hash.c lpsolve/lp_lib.c lpsolve/lp_wlp.c lpsolve/lp_matrix.c lpsolve/lp_mipbb.c lpsolve/lp_MPS.c lpsolve/lp_params.c lpsolve/lp_presolve.c lpsolve/lp_price.c lpsolve/lp_pricePSE.c lpsolve/lp_report.c lpsolve/lp_scale.c lpsolve/lp_simplex.c lpsolve/lp_SOS.c lpsolve/lp_utils.c lpsolve/yacc_read.c
set lplibpath=lpsolve/lib/liblpsolve55.a

set cc_flags=-DBFP_CALLMODEL=__stdcall -DYY_NEVER_INTERACTIVE -DPARSER_LP -DINVERSE_ACTIVE=INVERSE_LUSOL -DRoleIsExternalInvEngine

rem gcc -Wall -Wextra -O0 -ggdb %incl_test% %incl_lplib% %cc_flags% %src_test% %src_lplib% -lm -o dynamicTestPL

gcc -O0 -Wno-unused-function %incl_test% %incl_lplib% %cc_flags% %src_test% -lm -o dynamicTestPL %lplibpath%

gcc -O0 -Wno-unused-function -I ../src -I ../src/lina staticTestPL.c ../src/SIMPLEx.c ../src/lina/lina.c ../src/lina/qr.c -o staticTestPL.exe -lm

staticTestPL.exe
dynamicTestPL.exe