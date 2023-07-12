#include "SIMPLEx.h"
#include "lina/lina.h"

#define M 3
#define N 9

int main()
{

    double NBT[M][N-M] = {
                        {1,1,1,1,1,1},
                        {2,-1,-2,1,0,0},
                        {0,0,1,1,2,1}
                    };

    lina_transpose((double *)NBT,(double *)NBT,M,N-M);

    double B[M][M] = {
                        {1,0,0},
                        {0,1,0},
                        {0,0,1}
                    };
    
    double Cnb[N-M] = {-1,-2,1,-1,-4,2};
    double Cb[M] = {0,0,0};

    uint16_t nbIDx[N-M] ={1,2,3,4,5,6};
    uint16_t bIDx[M] ={7,8,9};

    double b[M] = {6,4,4};

    simplex_saveToStream(stdout, M, N, (double *)NBT, (double *)B, Cnb, nbIDx, Cb, bIDx, b);

    getc(stdin);

    simplex_solve(M, N, (double *)NBT, (double *)B, Cnb, nbIDx, Cb, bIDx, b);

    


}