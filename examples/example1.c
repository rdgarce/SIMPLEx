#include "SIMPLEx.h"
#include <stdio.h>

#define M 4
#define N 8

int main()
{
    double AT[N][M] = {
                        {1,0,2,4},
                        {0,1,5,2},
                        {1,0,0,0},
                        {0,1,0,0},
                        {0,0,-1,0},
                        {0,0,0,-1},
                        {0,0,1,0},
                        {0,0,0,1},
                    };

    double Zcnb[N-M] = {-1,-1,0,0};
    double Zcb[M] = {0,0,0,0};

    double Wcnb[N-M] = {-6,-7,1,1};
    double Wcb[M] = {0,0,0,0};
    double Wb = -18;

    uint16_t nbIDx[N-M] ={0,1,4,5};
    uint16_t bIDx[M] ={2,3,6,7};

    double b[M] = {6,4,10,8};
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