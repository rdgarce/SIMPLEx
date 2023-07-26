#include <stdint.h>
#include <stdlib.h>
#include "SIMPLEx_test.h"

PLtest_t *create_PLtest(uint16_t m, uint16_t n)
{
    PLtest_t *test = (PLtest_t *)malloc(sizeof(PLtest_t));
    if (!test)
        return NULL;
    
    test->AT = (double *)malloc(sizeof(double) * n * m);
    test->b = (double *)malloc(sizeof(double) * m);
    test->b_opt = (double *)malloc(sizeof(double) * m);
    test->bIDx = (uint16_t *)malloc(sizeof(uint16_t) * m);
    test->bIDx_opt = (uint16_t *)malloc(sizeof(uint16_t) * m);
    test->nbIDx = (uint16_t *)malloc(sizeof(uint16_t) * (n-m));
    test->Wcb = (double *)malloc(sizeof(double) * m);
    test->Wcnb = (double *)malloc(sizeof(double) * (n-m));
    test->Zcb = (double *)malloc(sizeof(double) * m);
    test->Zcnb = (double *)malloc(sizeof(double) * (n-m));

    if (!test->AT || !test->b || !test->b_opt || !test->bIDx || !test->bIDx_opt || 
            !test->nbIDx || !test->Wcb || !test->Wcnb || !test->Zcb || !test->Zcnb)
    {
        free(test->AT);
        free(test->b);
        free(test->b_opt);
        free(test->bIDx);
        free(test->bIDx_opt);
        free(test->nbIDx);
        free(test->Wcb);
        free(test->Wcnb);
        free(test->Zcb);
        free(test->Zcnb);
        return NULL;
    }
    test->m = m;
    test->n = n;
    return test;
}

void free_PLtest(PLtest_t *test)
{
    if (!test)
        return;
    
    free(test->AT);
    free(test->b);
    free(test->b_opt);
    free(test->bIDx);
    free(test->bIDx_opt);
    free(test->nbIDx);
    free(test->Wcb);
    free(test->Wcnb);
    free(test->Zcb);
    free(test->Zcnb);

    free(test);
}