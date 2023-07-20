#include <stdint.h>

typedef struct PLtest
{
    uint16_t m;
    uint16_t n;
    uint16_t *nbIDx;
    uint16_t *bIDx;
    uint16_t *bIDx_opt;
    double *AT;
    double *b;
    double *b_opt;
    double *Zcnb;
    double *Zcb;
    double Zopt;
    double Zopt_opt;
    double *Wcnb;
    double *Wcb;
    double Wb;
    int ret;
}PLtest_t;

PLtest_t *create_PLtest(uint16_t m, uint16_t n);
void free_PLtest(PLtest_t *test);