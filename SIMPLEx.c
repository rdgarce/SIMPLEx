#include "SIMPLEx.h"
#include "lina/lina.h"
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include <float.h>

int simplex_solve(uint16_t m, uint16_t n, double *NBT, double *B, double *Cnb, uint16_t *nbIDx, 
                                                        double *Cb, uint16_t *bIDx, double *b)
{
    // n has to be bigger than m because n = vars + slacks + artificials while m <= slacks + artificial
    assert(n > m);
    // all pointers must be non NULL
    assert(NBT && B && Cnb && nbIDx && Cb && bIDx && b);
    // all b elements must be non negative and all bIDx must be non negative and less than n
    for (uint16_t i = 0; i < m; i++)
        assert(b[i] >= 0 && bIDx[i] >= 0 && bIDx[i] <= n);
    // all nbIDx must be non negative and less than n
    for (uint16_t i = 0; i < n-m; i++)
        assert(nbIDx[i] >= 0 && nbIDx[i] <= n);

    // ----- END ASSERT ZONE ----- //

    // generic purpose m sized buffer for operations
    double *m_size_buffer = (double *)malloc(sizeof(double)*m);
    double *n_minus_m_size_buffer = (double *)malloc(sizeof(double)*(n-m));
    if (!m_size_buffer || !n_minus_m_size_buffer)
        return -2;

    double *m_size_buffer_old_ptr = m_size_buffer;
    double *b_old_ptr = b;
    double *temp_ptr;
    double temp_val;
    uint16_t temp_uint_val;

    int end_cause = 1;
    while (end_cause > 0)
    {
        // find the non-basic variable that has to enter the base according to the Bland rule:
        // the first negative Cnb value found determines the entering variable
        uint16_t CnbIDx = 0;

        while (Cnb[CnbIDx] >= 0 && CnbIDx < n-m)
            CnbIDx++;
        
        // no entering variable found, means we reach the optimal solution.
        // exiting with end_cause 0
        if (CnbIDx == n-m)
        {
            end_cause = 0;
            break;
        }

        // calculate the new b vector with m_size_buffer = B * b
        lina_dot(B,b,m_size_buffer,m,m,1);
        // swap b and m_size_buffer pointer so that b always point to new calculated b
        temp_ptr = b;
        b = m_size_buffer;
        m_size_buffer = temp_ptr;

        // update the CnbIDx row in the NBT according to NBT[CnbIDx] = B * NBT[CnbIDx]
        // m_size_buffer used to temporary store the result that is the new NBT[CnbIDx]
        lina_dot(B,&NBT[CnbIDx*m],m_size_buffer,m,m,1);

        // find the basic variable that has to exit the base according to the Bland rule:
        // find the basic variable with the minimum b_i/a_{ij} with a_{ij} > 0
        // and if multiple variable have the same b_i/a_{ij} than select the one with the smallest row index
        uint16_t NBColIDx;
        double min_value = DBL_MAX;
        bool foundExitVar = false;

        for (uint16_t i = 0; i < m; i++)
        {
            // a_{ij} <= 0, then continue
            if (m_size_buffer[i] <= 0)
                continue;
            
            temp_val = b[i]/m_size_buffer[i];
            foundExitVar = true;

            if (temp_val < min_value)
            {
                NBColIDx = i;
                min_value = temp_val;
            }
        }
        // if no exit var is found we exit because the problem is unlimited
        if (!foundExitVar)
        {
            end_cause = -1;
            break;
        }

        // now we pivot the NBT[CnbIDx] row and so we scale the corresponding rows in B matrix

        // if this is the column of the exiting variable we have to just scale the corresponding B row 
        // by the value of NBT[CnbIDx][NBColIDx]
        lina_scale(&B[NBColIDx*m],&B[NBColIDx*m],1.0/NBT[CnbIDx*m + NBColIDx],1,m);
        NBT[CnbIDx*m + NBColIDx] = 1;

        for (uint16_t i = 0; i < NBColIDx; i++)
        {
            if (NBT[CnbIDx*m + i] != 0)
            {
                lina_scale(&B[i*m],&B[i*m],-1.0/NBT[CnbIDx*m + i],1,m);
                lina_add(&B[NBColIDx*m],&B[i*m], &B[i*m],1,m);
                NBT[CnbIDx*m + i] = 0;
            }
        }
        for (uint16_t i = NBColIDx+1; i < m; i++)
        {
            if (NBT[CnbIDx*m + i] != 0)
            {
                lina_scale(&B[i*m],&B[i*m],-1.0/NBT[CnbIDx*m + i],1,m);
                lina_add(&B[NBColIDx*m],&B[i*m], &B[i*m],1,m);
                NBT[CnbIDx*m + i] = 0;
            }
        }
        
        // we have an updated NBT and B (that is B^(-1)), we now have to swap coefficients and indices of
        // entering and exiting variables
        temp_val = Cnb[CnbIDx];
        Cnb[CnbIDx] = Cb[NBColIDx];
        Cb[NBColIDx] = temp_val;
        temp_uint_val = nbIDx[CnbIDx];
        nbIDx[CnbIDx] = bIDx[NBColIDx];
        bIDx[NBColIDx] = temp_uint_val;

        // calculate Cb * B^(-1) and store result in m_size_buffer
        lina_dot(Cb,B,m_size_buffer,1,m,m);

        // we have to calculate m_size_buffer * NB.
        // we have NBT (NB^T) so (m_size_buffer * NB)^T = NBT * m_size_buffer^T
        // since m_size_buffer is a row vector, we have to consider it now a column vector
        // lina_dot would consider it memory contiguous anyway, so it's free!
        // the result is stored in n_minus_m_size_buffer
        lina_dot(NBT,m_size_buffer,n_minus_m_size_buffer,n-m,m,1);

        // now we calculate the updated Cnb = Cnb - n_minus_m_size_buffer
        lina_scale(n_minus_m_size_buffer,n_minus_m_size_buffer,-1.0,1,n-m);
        lina_add(Cnb,n_minus_m_size_buffer,Cnb,1,n-m);
    }
    
    // at the end we must copy the last b vector to the caller location
    if (b != b_old_ptr)
    {
        for (uint16_t i = 0; i < m; i++)
            b_old_ptr[i] = b[i];
    }
    
    free(m_size_buffer_old_ptr);
    free(n_minus_m_size_buffer);
    return end_cause;
}

void simplex_saveToStream(FILE *st, uint16_t m, uint16_t n, double *NBT, double *B, double *Cnb, 
                                        uint16_t *nbIDx, double *Cb, uint16_t *bIDx, double *b)
{
    fprintf(st, "       ");
    for (uint16_t i = 1; i <= n; i++)
    {
        fprintf(st, "x%d      ",i);
    }
    fprintf(st, "-z      b\n");

    for (uint16_t i = 0; i < m; i++)
    {
        fprintf(st, "x%d | ",bIDx[i]);
        fprintf(st, "%.3f",NBT[i]);
        for (uint16_t j = 1; j < n-m; j++)
        {
            NBT[j*m+i] >= 0 ? fprintf(st, "   %.3f",NBT[j*m+i]) : fprintf(st, "  %.3f",NBT[j*m+i]);
        }
        for (uint16_t k = 0; k < m; k++)
        {
            B[i*m + k] >= 0 ? fprintf(st, "   %.3f",B[i*m + k]) : fprintf(st, "  %.3f",B[i*m + k]);
        }
        fprintf(st, "   %.3f   %.3f\n",0.0, b[i]);        
    }
    
    fprintf(st, "-z | ");
    bool found = false;
    uint16_t j,k;
    j = k = 0;
    while (!found && j < n-m)
    {
        if (nbIDx[j] == 1)
        {
            fprintf(st, "%.3f",Cnb[j]);
            found = true;
        }
        j++;
    }
    while (!found && k < m)
    {
        if (bIDx[k] == 1)
        {
            fprintf(st, "%.3f",Cb[k]);
            found = true;
        }
        k++;
    }

    for (uint16_t i = 2; i <= n; i++)
    {
        found = false;
        j = k = 0;
        while (!found && j < n-m)
        {
            if (nbIDx[j] == i)
            {
                Cnb[j] >= 0 ? fprintf(st, "   %.3f",Cnb[j]) : fprintf(st, "  %.3f",Cnb[j]);
                found = true;
            }
            j++;
        }
        while (!found && k < m)
        {
            if (bIDx[k] == i)
            {
                Cb[k] >= 0 ? fprintf(st, "   %.3f",Cb[k]) : fprintf(st, "  %.3f",Cb[k]);
                found = true;
            }
            k++;
        }
    }

    fprintf(st, "   %.3f   ", 1.0); 
    
    // -z = -Cb*B*b
    double *buff = (double *)malloc(sizeof(double) * m);
    if (!buff)
    {
        fprintf(st, "MALLOC ERROR!");
        return;
    }
    
    double z;
    lina_dot(Cb,B,buff,1,m,m);
    lina_dot(buff,b,&z,1,m,1);
    if (z != 0)
        z = -z;

    fprintf(st, "%.3f\n",z);
    free(buff);
}