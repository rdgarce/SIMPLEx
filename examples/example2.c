#include "SIMPLEx.h"
#include <stdio.h>

#define M 3
#define N 7

int main()
{
    double AT[N][M] = {
                        {1.0/4.0, 1.0/2.0, 0},
                        {-8, -12, 0},
                        {-1, -1.0/2.0, 1},
                        {9, 3, 0},
                        {1, 0, 0},
                        {0, 1, 0},
                        {0, 0, 1},
                    };

    double Zcnb[N-M] = {-3.0/4.0, 20, -1.0/2.0, 6};
    double Zcb[M] = {0, 0, 0};

    double *Wcnb = NULL;
    double *Wcb = NULL;
    double Wb = 0;

    uint16_t nbIDx[N-M] = {0, 1, 2, 3};
    uint16_t bIDx[M] = {4, 5, 6};

    double b[M] = {0, 0, 1};
    double Zopt;

    if(!simplex_revisedSolve(M,N, (double *)AT, (double *)b, (uint16_t *)nbIDx, (uint16_t *)bIDx,
                            (double *)Zcnb, (double *)Zcb, (double *)Wcnb, (double *)Wcb, Wb,&Zopt))
    {
        printf("Optimal solution found:\n");
        for (uint16_t i = 0; i < M; i++)
            printf("x%d = %.3f\n",bIDx[i],b[i]);
        
        printf("Value: %.3f\n",Zopt);
    }
}