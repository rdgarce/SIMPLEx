#include "SIMPLEx.h"
#include "lina/lina.h"
#include <stdint.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#ifdef SIMPLEX_DEBUG
    #include <assert.h>
#endif

/*
*   Pivot the rows of the B matrix around the index-th value of the p vector.
*   B is mxm and p is 1xm.
*/
static void pivot_matrix(double *p, double *B, uint16_t m, uint16_t index)
{
    for (uint16_t i = 0; i < index; i++)
    {
        if (p[i] != 0)
        {
            for (uint16_t j = 0; j < m; j++)
                B[i*m + j] -= B[index*m + j]*p[i]/p[index];            
        }
    }
    for (uint16_t i = index+1; i < m; i++)
    {
        if (p[i] != 0)
        {
            for (uint16_t j = 0; j < m; j++)
                B[i*m + j] -= B[index*m + j]*p[i]/p[index];
        }
    }
    for (uint16_t j = 0; j < m; j++)
        B[index*m + j] *= 1.0/p[index];   
}

/*
*  Given the following linear programming problem in the canonical form:
*   
*   min z = Zc * x
*   subject to: 
*   A * x = b
*   x >= 0
*
*  Execute the two-phase revised simplex algorithm.
*  
*   INPUT PARAMETERS:
*       - m:     the number of constraints of the problem
*       - n:     the number of variables of the system of equations in 
*                   canonical form (n also counts slack and artificial variables)
*       - AT:    the nxm matrix corresponding to the transposed of the A matrix 
*                   of the system
*       - b:     the 1xm vector of the values of the basic variables in the 
*                   initial basic feasible solution
*       - nbIDx: the 1x(n-m) vector of the indices of the A columns (AT rows)
*                   representing the non-basic variables in the initial basic
*                   feasible solution
*       - bIDx:  the 1xm vector of the indices of the A columns (AT rows)
*                   representing the basic variables in the initial basic
*                   feasible solution
*       - Zcnb:  the 1x(n-m) vector of the coefficients of the non-basic variables
*                   of the objective function in the initial basic feasible solution
*       - Zcb:   the 1xm vector of the coefficients of the basic variables
*                   of the objective function in the initial basic feasible solution
*       - Wcnb:  the 1x(n-m) vector of the coefficients of the non-basic variables
*                   of the artificial objective function in the initial basic 
*                   feasible solution. If this parameter is set to a NULL pointer,
*                   only phase 2 of the algorithm is executed
*       - Wcb:   the 1xm vector of the coefficients of the basic variables
*                   of the artificial objective function in the initial basic 
*                   feasible solution. If this parameter is set to a NULL pointer,
*                   only phase 2 of the algorithm is executed
*       - Wb:    the scalar containing the numerical value of the artificial objective
*                   function in the initial basic feasible solution
*   
*   OUTPUT PARAMETERS:
*       - b:     the 1xm vector of the values of the basic variables in the 
*                   last basic feasible solution
*       - nbIDx: the 1x(n-m) vector of the indices of the A columns (AT rows)
*                   representing the non-basic variables in the last basic
*                   feasible solution
*       - bIDx:  the 1xm vector of the indices of the A columns (AT rows)
*                   representing the basic variables in the last basic
*                   feasible solution
*       - Zcnb:  the 1x(n-m) vector of the coefficients of the non-basic variables
*                   of the objective function in the last basic feasible solution
*       - Zcb:   the 1xm vector of the coefficients of the basic variables
*                   of the objective function in the last basic feasible solution
*       - Zopt:  the scalar containing the numerical value of the optimal solution
*                   of the objective function. This parameter should be considered only
*                   if the return value is 0
*
*   RETURNS:
*       - (0):   the input problem is successfully minimized and the optimal solution
*                   is found. The output parameters contain the structure of the optimal
*                   solution.
*       - (-1):  the problem is unlimited
*       - (-2):  the problem is infeasible
*       - (-3):  a memory allocation problem occurred
*
*   IMPORTANT NOTES:
*       - All the input and output parameters memory allocation must be done by the caller
*       - All the input and output matrices and vectors must be allocated in memory-contiguous 
*            locations and values stored in row-major order
*/
int simplex_revisedSolve(uint16_t m, uint16_t n, double *AT, double *b, uint16_t *nbIDx, uint16_t *bIDx,
                            double *Zcnb, double *Zcb, double *Wcnb, double *Wcb, double Wb,double *Zopt)
{
    #ifdef SIMPLEX_DEBUG
    // n has to be bigger than m because n = vars + slacks + artificials while m <= slacks + artificial
    assert(n > m);
    // all pointers must be non NULL (except Wcnb and Wcb)
    assert(AT && b && nbIDx && bIDx && Zcnb && Zcb && Zopt);
    // all b elements must be non negative and all bIDx must be non negative and less than n
    for (uint16_t i = 0; i < m; i++)
        assert(b[i] >= 0 && bIDx[i] >= 0 && bIDx[i] < n);
    // all nbIDx must be non negative and less than n
    for (uint16_t i = 0; i < n-m; i++)
        assert(nbIDx[i] >= 0 && nbIDx[i] <= n);
    #endif

    double *B_1 = (double *)calloc(m*m,sizeof(double));
    double *memory_buffer = (double *)malloc(sizeof(double)*(m + n));
    // Bytemask tale che se lo i-esimo elemento = 1, allora
    // la var con IDx = i non può entrare in base
    bool *no_entering_varIDx_mask = (bool *)calloc(n,sizeof(bool));
    
    // returning condition initialized at -3
    int exit_condition = -3;

    // Allocation error
    // returns -3
    if (!B_1 || !memory_buffer || !no_entering_varIDx_mask)
        goto free_and_return;
    
    double *m_size_buffer = memory_buffer;
    double *m_size_buffer_2 = memory_buffer + m;
    double *n_minus_m_size_buffer = memory_buffer + (2*m);
    double simple_buffer;

    // B_1 matrix is set to be an identity matrix for later uses
    for (uint16_t i = 0; i < m; i++)
        B_1[i*m + i] = 1;
    

    // phase 1
    if (Wcnb && Wcb)
    {
        // Wcnb and Wcb are not NULL, so we execute the phase 1: 
        // min w = [ Wcnb | Wcb ] * [ xnb | xb]^T + Wb

        // first iteration of the phase 1 is not in a while cicle because on first iteration:
        // - We have not to calculate b' = B_1 * b
        // - We have not to calculate p_x' = B_1 * p_x (the entering variable new AT's row)

        // search for the first entering variable according to the Bland rule:
        // the first negative Wcnb value found determines the entering variable
        uint16_t cnb_index = 0;
        while ((Wcnb[cnb_index] >= 0 || no_entering_varIDx_mask[nbIDx[cnb_index]]) && cnb_index < n-m)
            cnb_index++;
        
        if (cnb_index == n-m)
        {
            // no entering variable found.
            // check for the obj. function value
            lina_dot(Wcb,b,&simple_buffer,1,m,1);
            simple_buffer = simple_buffer - Wb;
            
            if (fabs(simple_buffer) < FLT_EPSILON)
                goto phase2;
            else
            {
                exit_condition = -2;
                goto free_and_return;
            }
        }

        uint16_t enteringVarIDx = nbIDx[cnb_index];
        
        // find the basic variable that has to exit the base according to the Bland rule:
        // find the basic variable with the minimum b_i/a_{ij} with a_{ij} > 0
        // and if multiple variable have the same b_i/a_{ij} than select the one with the smallest row index
        uint16_t at_index;
        double min_value = DBL_MAX;
        bool found = false;

        for (uint16_t i = 0; i < m; i++)
        {
            if(AT[enteringVarIDx*m + i] <= 0)
                continue;
            
            simple_buffer = b[i]/AT[enteringVarIDx*m + i];
            found = true;

            if (simple_buffer < min_value)
            {
                at_index = i;
                min_value = simple_buffer;
            }
        }
        // exit variable not found: problem is unlimited
        if (!found)
        {
            exit_condition = -1;
            goto free_and_return;
        }

        uint16_t exitingVarIDx = bIDx[at_index];
        
        // Update nbIDx and bIDx
        nbIDx[cnb_index] = exitingVarIDx;
        bIDx[at_index] = enteringVarIDx;
        // Update Wcnb and Wcb
        simple_buffer = Wcnb[cnb_index];
        Wcnb[cnb_index] = Wcb[at_index];
        Wcb[at_index] = simple_buffer;
        // Update Zcnb and Zcb
        simple_buffer = Zcnb[cnb_index];
        Zcnb[cnb_index] = Zcb[at_index];
        Zcb[at_index] = simple_buffer;
        // Set the no entering variable flag
        no_entering_varIDx_mask[exitingVarIDx] = true;

        // now we have the updated Wcb and Wcnb
        // we compute the new B^(-1) matrix pivoting around the p vector.
        // the first time the p vector is the one from the AT matrix while
        // on the next iteration we must calculate it before
        pivot_matrix(&AT[enteringVarIDx*m],B_1,m,at_index);

        // Compute m_size_buffer = Wcb * B_1 (Cb * B^(-1))
        lina_dot(Wcb,B_1,m_size_buffer,1,m,m);

        // we should calculate m_size_buffer * NB but we have NB^T, a submatrix of AT
        // (CbB^(-1)*NB)^T = NB^T*CbB^(-1) = NBT * m_size_buffer
        // Compute n_minus_m_size_buffer = NBT * m_size_buffer^T
        // NBT is the matrix composed of the non-basic variable rows of the AT
        // so we execute the dot product row by row to select the particular rows of AT
        // that refers to non-basic variables
        for (uint16_t i = 0; i < n-m; i++)
            lina_dot(&AT[nbIDx[i]*m],m_size_buffer,&n_minus_m_size_buffer[i],1,m,1);

        // Compute n_minus_m_size_buffer = Wcnb - n_minus_m_size_buffer (Wcnb' = Wcnb - Wcb * B^(-1))
        //lina_scale(n_minus_m_size_buffer,n_minus_m_size_buffer,-1.0,1,n-m);
        lina_add(Wcnb,n_minus_m_size_buffer,n_minus_m_size_buffer,-1.0,1,n-m);

        // now n_minus_m_size_buffer contains the Wcnb' vector

        // the first iteration is completed, now we enter the loop that does all the steps
        while (1)
        {
            // search for the first entering variable according to the Bland rule:
            // the first negative Wcnb value found determines the entering variable
            uint16_t cnb_index = 0;
            while ((n_minus_m_size_buffer[cnb_index] >= 0 || no_entering_varIDx_mask[nbIDx[cnb_index]]) && cnb_index < n-m)
                cnb_index++;
            
            if (cnb_index == n-m)
            {
                // no entering variable found.
                // check for the obj. function value
                // w = CbB^(-1)b - Wb
                
                // Compute m_size_buffer = Wcb * B_1 (Cb * B^(-1))
                lina_dot(Wcb,B_1,m_size_buffer,1,m,m);
                // Compute simple_buffer = m_size_buffer * b
                lina_dot(m_size_buffer,b,&simple_buffer,1,m,1);
                simple_buffer = simple_buffer - Wb;
                
                if (fabs(simple_buffer) < FLT_EPSILON)
                    goto phase2;
                else
                {
                    exit_condition = -2;
                    goto free_and_return;
                }
            }

            uint16_t enteringVarIDx = nbIDx[cnb_index];

            // calculate b' = B_1 * b
            lina_dot(B_1,b,m_size_buffer,m,m,1);
            // calculate p_x' = B_1 * p_x (the entering variable new AT's row)
            lina_dot(B_1,&AT[enteringVarIDx*m],m_size_buffer_2,m,m,1);

            // find the basic variable that has to exit the base according to the Bland rule:
            // find the basic variable with the minimum b_i/a_{ij} with a_{ij} > 0
            // and if multiple variable have the same b_i/a_{ij} than select the one with the smallest row index
            uint16_t at_index;
            double min_value = DBL_MAX;
            bool found = false;

            for (uint16_t i = 0; i < m; i++)
            {
                if(m_size_buffer_2[i] <= 0)
                    continue;
                
                simple_buffer = m_size_buffer[i]/m_size_buffer_2[i];
                found = true;

                if (simple_buffer < min_value)
                {
                    at_index = i;
                    min_value = simple_buffer;
                }
            }
            // exit variable not found: problem is unlimited
            if (!found)
            {
                exit_condition = -1;
                goto free_and_return;
            }

            uint16_t exitingVarIDx = bIDx[at_index];

            // Update nbIDx and bIDx
            nbIDx[cnb_index] = exitingVarIDx;
            bIDx[at_index] = enteringVarIDx;
            // Update Wcnb and Wcb
            simple_buffer = Wcnb[cnb_index];
            Wcnb[cnb_index] = Wcb[at_index];
            Wcb[at_index] = simple_buffer;
            // Update Zcnb and Zcb
            simple_buffer = Zcnb[cnb_index];
            Zcnb[cnb_index] = Zcb[at_index];
            Zcb[at_index] = simple_buffer;
            // Set the no entering variable flag
            no_entering_varIDx_mask[exitingVarIDx] = true;
            
            // now we have the updated Wcb and Wcnb
            // we compute the new B^(-1) matrix pivoting around the p vector.
            // the first time the p vector is the one from the AT matrix while
            // on the next iteration we must calculate it before
            pivot_matrix(m_size_buffer_2,B_1,m,at_index);

            // Compute m_size_buffer = Wcb * B_1 (Cb * B^(-1))
            lina_dot(Wcb,B_1,m_size_buffer,1,m,m);

            // we should calculate m_size_buffer * NB but we have NB^T, a submatrix of AT
            // (CbB^(-1)*NB)^T = NB^T*CbB^(-1) = NBT * m_size_buffer
            // Compute n_minus_m_size_buffer = NBT * m_size_buffer^T
            // NBT is the matrix composed of the non-basic variable rows of the AT
            // so we execute the dot product row by row to select the particular rows of AT
            // that refers to non-basic variables
            for (uint16_t i = 0; i < n-m; i++)
                lina_dot(&AT[nbIDx[i]*m],m_size_buffer,&n_minus_m_size_buffer[i],1,m,1);

            // Compute n_minus_m_size_buffer = Wcnb - n_minus_m_size_buffer (Wcnb' = Wcnb - Wcb * B^(-1))
            //lina_scale(n_minus_m_size_buffer,n_minus_m_size_buffer,-1.0,1,n-m);
            lina_add(Wcnb,n_minus_m_size_buffer,n_minus_m_size_buffer,-1.0,1,n-m);
        }
    }

    // phase 2
    phase2:

    //before starting the while cicle, we must calculate the Zcnb'
    // Compute m_size_buffer = Zcb * B_1 (Cb * B^(-1))
    lina_dot(Zcb,B_1,m_size_buffer,1,m,m);

    // we should calculate m_size_buffer * NB but we have NB^T, a submatrix of AT
    // (CbB^(-1)*NB)^T = NB^T*CbB^(-1) = NBT * m_size_buffer
    // Compute n_minus_m_size_buffer = NBT * m_size_buffer^T
    // NBT is the matrix composed of the non-basic variable rows of the AT
    // so we execute the dot product row by row to select the particular rows of AT
    // that refers to non-basic variables
    for (uint16_t i = 0; i < n-m; i++)
        lina_dot(&AT[nbIDx[i]*m],m_size_buffer,&n_minus_m_size_buffer[i],1,m,1);

    // Compute n_minus_m_size_buffer = Zcnb - n_minus_m_size_buffer (Zcnb' = Zcnb - Zcb * B^(-1))
    //lina_scale(n_minus_m_size_buffer,n_minus_m_size_buffer,-1.0,1,n-m);
    lina_add(Zcnb,n_minus_m_size_buffer,n_minus_m_size_buffer,-1.0,1,n-m);

    while (1)
    {
        // search for the first entering variable according to the Bland rule:
        // the first negative Zcnb value found determines the entering variable
        uint16_t cnb_index = 0;
        while ((n_minus_m_size_buffer[cnb_index] >= 0 || no_entering_varIDx_mask[nbIDx[cnb_index]]) && cnb_index < n-m)
            cnb_index++;
        
        // if no entering variable is found then the optimal solution is found
        if (cnb_index == n-m)
        {
            exit_condition = 0;
            // compute b_opt and z_opt
            lina_dot(B_1,b,m_size_buffer,m,m,1);
            for (uint16_t i = 0; i < m; i++)
                b[i] = m_size_buffer[i];
            
            lina_dot(Zcb,b,Zopt,1,m,1);

            goto free_and_return;
        }

        uint16_t enteringVarIDx = nbIDx[cnb_index];

        // calculate b' = B_1 * b
        lina_dot(B_1,b,m_size_buffer,m,m,1);
        // calculate p_x' = B_1 * p_x (the entering variable new AT's row)
        lina_dot(B_1,&AT[enteringVarIDx*m],m_size_buffer_2,m,m,1);

        // find the basic variable that has to exit the base according to the Bland rule:
        // find the basic variable with the minimum b_i/a_{ij} with a_{ij} > 0
        // and if multiple variable have the same b_i/a_{ij} than select the one with the smallest row index
        uint16_t at_index;
        double min_value = DBL_MAX;
        bool found = false;

        for (uint16_t i = 0; i < m; i++)
        {
            if(m_size_buffer_2[i] <= 0)
                continue;
            
            simple_buffer = m_size_buffer[i]/m_size_buffer_2[i];
            found = true;

            if (simple_buffer < min_value)
            {
                at_index = i;
                min_value = simple_buffer;
            }
        }
        // exit variable not found: problem is unlimited
        if (!found)
        {
            exit_condition = -1;
            goto free_and_return;
        }

        uint16_t exitingVarIDx = bIDx[at_index];

        // Update nbIDx and bIDx
        nbIDx[cnb_index] = exitingVarIDx;
        bIDx[at_index] = enteringVarIDx;
        // Update Zcnb and Zcb
        simple_buffer = Zcnb[cnb_index];
        Zcnb[cnb_index] = Zcb[at_index];
        Zcb[at_index] = simple_buffer;

        // now we have the updated Wcb and Wcnb
        // we compute the new B^(-1) matrix pivoting around the p vector.
        // the first time the p vector is the one from the AT matrix while
        // on the next iteration we must calculate it before
        pivot_matrix(m_size_buffer_2,B_1,m,at_index);

        // Compute m_size_buffer = Wcb * B_1 (Cb * B^(-1))
        lina_dot(Zcb,B_1,m_size_buffer,1,m,m);

        // we should calculate m_size_buffer * NB but we have NB^T, a submatrix of AT
        // (CbB^(-1)*NB)^T = NB^T*CbB^(-1) = NBT * m_size_buffer
        // Compute n_minus_m_size_buffer = NBT * m_size_buffer^T
        // NBT is the matrix composed of the non-basic variable rows of the AT
        // so we execute the dot product row by row to select the particular rows of AT
        // that refers to non-basic variables
        for (uint16_t i = 0; i < n-m; i++)
            lina_dot(&AT[nbIDx[i]*m],m_size_buffer,&n_minus_m_size_buffer[i],1,m,1);

        // Compute n_minus_m_size_buffer = Wcnb - n_minus_m_size_buffer (Wcnb' = Wcnb - Wcb * B^(-1))
        //lina_scale(n_minus_m_size_buffer,n_minus_m_size_buffer,-1.0,1,n-m);
        lina_add(Zcnb,n_minus_m_size_buffer,n_minus_m_size_buffer,-1.0,1,n-m);
    }

    free_and_return:
    free(B_1);
    free(memory_buffer);
    free(no_entering_varIDx_mask);
    return exit_condition;
}