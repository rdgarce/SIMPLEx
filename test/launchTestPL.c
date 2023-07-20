#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>
#include "testsPL.h"
#include "SIMPLEx.h"

int main()
{

    uint16_t total_num_tests = 0, passed_tests = 0;
    int res;
    bool test_success;
    for (uint16_t i = 0; i < N_TESTS; i++)
    {
        res = simplex_revisedSolve( PLtest_array[i].m,
                                    PLtest_array[i].n,
                                    (double *)PLtest_array[i].AT,
                                    (double *)PLtest_array[i].b,
                                    PLtest_array[i].nbIDx,
                                    PLtest_array[i].bIDx,
                                    (double *)PLtest_array[i].Zcnb,
                                    (double *)PLtest_array[i].Zcb,
                                    (double *)PLtest_array[i].Wcnb,
                                    (double *)PLtest_array[i].Wcb,
                                    PLtest_array[i].Wb,
                                    &PLtest_array[i].Zopt
                                );
        test_success = true;
        total_num_tests++;

        if (res != PLtest_array[i].ret)
        {
            fprintf(stderr,"** Error in test n° %d: expected ret. val. %d while it was %d\n",i+1, PLtest_array[i].ret,res);
        }
        

        for (uint16_t j = 0; j < PLtest_array[i].m; j++)
        {
            if (fabs(PLtest_array[i].b[j] - PLtest_array[i].b_opt[j]) >= FLT_EPSILON)
            {
                fprintf(stderr,"** Error in test n° %d: b_%d = %.17f while bopt_%d = %.17f\n",i+1,j,PLtest_array[i].b[j],j,PLtest_array[i].b_opt[j]);
                test_success = false;
            }
        }
        if (fabs(PLtest_array[i].Zopt - PLtest_array[i].Zopt_opt) >= FLT_EPSILON)
        {
            fprintf(stderr,"** Error in test n° %d: Zopt = %.17f while Zopt_opt = %.17f\n",i+1, PLtest_array[i].Zopt, PLtest_array[i].Zopt_opt);
            for (uint16_t j = 0; j < PLtest_array[i].m; j++)
                fprintf(stderr, "\tx_%d = %.17f while x_opt_%d = %.17f\n", PLtest_array[i].bIDx[j], PLtest_array[i].b[j],
                                                                    PLtest_array[i].bIDx_opt[j], PLtest_array[i].b_opt[j]);
            test_success = false;
        }
        
        if (test_success)
            passed_tests++;
    }
    fprintf(stderr,"PL Tests executed: %d/%d passed\n",passed_tests,total_num_tests);
}

