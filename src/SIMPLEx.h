#include <stdint.h>

#ifndef SIMPLEx_H
#define SIMPLEx_H

// For debug purpose
//#define SIMPLEX_DEBUG

int simplex_revisedSolve(uint16_t m, uint16_t n, double *AT, double *b, uint16_t *nbIDx, uint16_t *bIDx,
                            double *Zcnb, double *Zcb, double *Wcnb, double *Wcb, double Wb,double *Zopt);

#endif