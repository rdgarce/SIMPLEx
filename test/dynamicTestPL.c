#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <stdbool.h>
#include <float.h>
#include "lp_lib.h"
#include "SIMPLEx_test.h"
#include "SIMPLEx.h"
#include "assert.h"

#define SIMPLEX_DEBUG

// The number of tests to execute
#define NUM_TESTS 5000
// MAX num of constraints
#define MAX_CONSTR_NUM 10
// MAX num of variables
#define MAX_VAR_NUM 10

// Min value of a variable coefficients and b values
#define MIN_VAR_COEFF 1
// Max value of a variable coefficients and b values
#define MAX_VAR_COEFF 10

// The average percentage of <= constraints in the lp problems
#define LE_CONSTR_PCT 100
// The average percentage of >= constraints in the lp problems
#define GE_CONSTR_PCT 0
// The average percentage of = constraints in the lp problems
#define EQ_CONSTR_PCT 0

static void fpritf_dbl_matrix(FILE *fp, double *mat, unsigned int rows, unsigned int cols);
static void fpritf_uint_matrix(FILE *fp, uint16_t *mat, unsigned int rows, unsigned int cols);
static bool check_SIMPLEx_lplib_solution(lprec *lplib_test, int lplib_ret_val, PLtest_t *SIMPLEx_test, int SIMPLEx_ret_val);

// Do not edit this!
#define RAND_MIN_MAX_ ((rand() % (MAX_VAR_COEFF+1)) + MIN_VAR_COEFF)

int main()
{
    // Just to be sure you know what you are doing ;)
    assert(LE_CONSTR_PCT + GE_CONSTR_PCT + EQ_CONSTR_PCT == 100);
    assert(MIN_VAR_COEFF < MAX_VAR_COEFF);

    PLtest_t *SIMPLEx_test;
    lprec *lplib_test;

    // the for cicle that executes all the tests
    for (unsigned int i = 0; i < NUM_TESTS; i++)
    {
        uint16_t num_constr = rand() % (MAX_CONSTR_NUM+1);
        if (num_constr == 0)
            num_constr+=1;
        
        uint16_t num_vars = (rand() % (MAX_VAR_NUM+1) % (UINT16_MAX - 2*num_constr)); // this has to be added with slacks and artifical to become the "n" of a SIMPLEx lp problem
        if (num_vars == 0)
            num_vars+=1;
        uint16_t n,m;

        // determines the number of LE, GE, EQ constraints
        uint16_t LE_num = 0, GE_num = 0, EQ_num = 0;
        for (uint16_t j = 0; j < num_constr; j++)
        {
            uint16_t random_num = rand() % 100;
            if (random_num < LE_CONSTR_PCT)
                LE_num++;
            else if (random_num >= LE_CONSTR_PCT && random_num < LE_CONSTR_PCT + GE_CONSTR_PCT)
                GE_num++;
            else
                EQ_num++;
        }

        // initialize the tests structures
        lplib_test = make_lp(0,num_vars);
        if (!lplib_test)
        {
            fprintf(stderr, "Error while allocating lplib testcase. Exiting...\n");
            return 1;
        }

        m = num_constr;
        n = num_vars + LE_num + 2*GE_num + EQ_num;

        //fprintf(stderr, "num_constr = %d\nnum_vars = %d\nn = %d\nm = %d\n",num_constr,num_vars,n,m);
        SIMPLEx_test = create_PLtest(m,n);
        if (!SIMPLEx_test)
        {
            fprintf(stderr, "Error while allocating SIMPLEx testcase. Exiting...\n");
            return 1;
        }

        double Wb = 0;
        // non basic Z obj func. values
        double *Zcnb = (double *)calloc(n-m+1, sizeof(double));
        // Zcb is all zeros
        double *Zcb = (double *)calloc(m,sizeof(double));
        // non basic W obj func. values
        double *Wcnb = (double *)calloc(n-m, sizeof(double));
        // Wcb is all zeros
        double *Wcb = (double *)calloc(m,sizeof(double));
        if (!Zcnb || !Zcb || !Wcnb || !Wcb)
        {
            fprintf(stderr, "Error while allocating the obj functions rows. Exiting...\n");
            return 1;
        }

        // Now we set the first num_vars nbIDx elements for SIMPLEx problem
        // The firsts non basic variable are the user variables
        for (unsigned int j = 0; j < num_vars; j++)
            SIMPLEx_test->nbIDx[j] = j;

        // We have to add all the constraints to both the problems
        double b_val;
        double *constr_row = (double *)malloc(sizeof(double) * (num_vars+1));
        if (!constr_row)
        {
            fprintf(stderr, "Error while allocating the constraint row values. Exiting...\n");
            return 1;
        }
        // this is done because lplib wants a 0 in first position
        constr_row[0] = 0;
        // First add the LE constraints
        for (unsigned int j = 0; j < LE_num; j++)
        {
            b_val = RAND_MIN_MAX_;
            for (unsigned int k = 0; k < num_vars; k++)
                constr_row[k+1] = RAND_MIN_MAX_;
            
            if (!add_constraint(lplib_test,constr_row,LE,b_val))
            {
                fprintf(stderr, "Error while adding a LE constraint to a lplib problem. Exiting...\n");
                return 1;
            }
            
            // Adding the first coefficients to the AT matrix
            for (unsigned int k = 0; k < num_vars; k++)
                SIMPLEx_test->AT[k*m + j] = constr_row[k+1];

            // Adding 0 from num_vars to num_vars+j
            for (unsigned int k = num_vars; k < num_vars + j; k++)
                SIMPLEx_test->AT[k*m + j] = 0;
            
            SIMPLEx_test->AT[(num_vars + j)*m + j] = 1;

            // Adding 0 from num_vars+j+1 to n
            for (unsigned int k = num_vars+j+1; k < n; k++)
                SIMPLEx_test->AT[k*m + j] = 0;

            // setting the b val
            SIMPLEx_test->b[j] = b_val;
            // setting the bIDx val
            SIMPLEx_test->bIDx[j] = num_vars+j;
        }

        // Now add the EQ constraints
        for (uint16_t j = LE_num; j < LE_num + EQ_num; j++)
        {
            b_val = RAND_MIN_MAX_;
            for (unsigned int k = 0; k < num_vars; k++)
                constr_row[k+1] = RAND_MIN_MAX_;
            
            if (!add_constraint(lplib_test,constr_row,EQ,b_val))
            {
                fprintf(stderr, "Error while adding a LE constraint to a lplib problem. Exiting...\n");
                return 1;
            }

            // Adding the first coefficients to the AT matrix
            for (unsigned int k = 0; k < num_vars; k++)
                SIMPLEx_test->AT[k*m + j] = constr_row[k+1];
            
            // Adding 0 from num_vars to num_vars+j
            for (unsigned int k = num_vars; k < num_vars + j; k++)
                SIMPLEx_test->AT[k*m + j] = 0;
            
            SIMPLEx_test->AT[(num_vars + j)*m + j] = 1;

            // Adding 0 from num_vars+j+1 to n
            for (unsigned int k = num_vars+j+1; k < n; k++)
                SIMPLEx_test->AT[k*m + j] = 0;

            // setting the b val
            SIMPLEx_test->b[j] = b_val;
            // setting the bIDx val
            SIMPLEx_test->bIDx[j] = num_vars+j;

            // incrementing the Wcnb function coefficients
            for (unsigned int k = 0; k < num_vars; k++)
            {
                Wcnb[k] -= constr_row[k+1];
            }
            // incrementing the Wb val
            Wb -= b_val;
        }

        // Now add the GE constraints
        for (uint16_t j = LE_num + EQ_num, z = num_vars; j < LE_num + EQ_num + GE_num; j++, z++)
        {
            unsigned int x_ = (j > LE_num + EQ_num);

            b_val = RAND_MIN_MAX_;
            for (unsigned int k = 0; k < num_vars; k++)
                constr_row[k+1] = RAND_MIN_MAX_;
            
            if (!add_constraint(lplib_test,constr_row,GE,b_val))
            {
                fprintf(stderr, "Error while adding a GE constraint to a lplib problem. Exiting...\n");
                return 1;
            }

            // Adding the first coefficients to the AT matrix
            for (unsigned int k = 0; k < num_vars; k++)
                SIMPLEx_test->AT[k*m + j] = constr_row[k+1];

            // Adding 0 from num_vars to num_vars+j
            for (unsigned int k = num_vars; k < num_vars + j + x_; k++)
                SIMPLEx_test->AT[k*m + j] = 0;
            
            SIMPLEx_test->AT[(num_vars + j + x_)*m + j] = -1;
            SIMPLEx_test->AT[(num_vars + j+1 + x_)*m + j] = 1;

            // Adding 0 from num_vars+j+2 to n
            for (unsigned int k = num_vars+j+2 + x_; k < n; k++)
                SIMPLEx_test->AT[k*m + j] = 0;

            // setting the b val
            SIMPLEx_test->b[j] = b_val;
            // setting the bIDx val
            SIMPLEx_test->bIDx[j] = num_vars+j+1 + x_;
            // setting the nbIDx val
            SIMPLEx_test->nbIDx[z] = num_vars+j + x_;

            // incrementing the Wcnb function coefficients
            for (unsigned int k = 0; k < num_vars; k++)
            {
                Wcnb[k] -= constr_row[k+1];
            }
            Wcnb[num_vars + j] = 1;
            // incrementing the Wb val
            Wb -= b_val;
        }
        
        // Create the objective functions
        Zcnb[0] = 0;
        for (unsigned int j = 0; j < num_vars; j++)
        {
            Zcnb[j+1] = RAND_MIN_MAX_;
        }
        
        if (!set_obj_fn(lplib_test, Zcnb))
        {
            fprintf(stderr, "Error while adding the objective function to a lplib problem. Exiting...\n");
            return 1;
        }

        // Choose if the problem is a MIN or a MAX
        // 0 = MIN, 1 = MAX
        uint8_t min_max;
        min_max = rand() % 2;
        if (min_max == 1)
        {
            set_maxim(lplib_test);
            SIMPLEx_test->min_max = true;
            for (unsigned int j = 0; j < n-m; j++)
                SIMPLEx_test->Zcnb[j] = -Zcnb[j+1];            
        }
        else
        {
            set_minim(lplib_test);
            SIMPLEx_test->min_max = false;
            for (unsigned int j = 0; j < n-m; j++)
                SIMPLEx_test->Zcnb[j] = Zcnb[j+1];
        }
        for (unsigned int j = 0; j < m; j++)
            SIMPLEx_test->Zcb[j] = Zcb[j];
        
        // If there are GE or EQ constraints, than we need to set the artificial obj func
        if (GE_num > 0 || EQ_num > 0)
        {
            for (unsigned int j = 0; j < n-m; j++)
                SIMPLEx_test->Wcnb[j] = Wcnb[j];
            
            for (unsigned int j = 0; j < m; j++)
                SIMPLEx_test->Wcb[j] = Wcb[j];
            
            SIMPLEx_test->Wb = Wb;
        }
        else
        {
            SIMPLEx_test->Wcnb = NULL;
            SIMPLEx_test->Wcb = NULL;
        }
        
        
        // Now we optimize both the problems
        
        set_verbose(lplib_test,0);
        
        //fprintf(stderr,"\nPrinting AT:\n");
        //fpritf_dbl_matrix(stderr,SIMPLEx_test->AT,n,m);

        /*fprintf(stderr,"\nPrinting b:\n");
        fpritf_dbl_matrix(stderr,SIMPLEx_test->b,1,m);

        fprintf(stderr,"\nPrinting nbIDx:\n");
        fpritf_uint_matrix(stderr,SIMPLEx_test->nbIDx,1,n-m);

        fprintf(stderr,"\nPrinting bIDx:\n");
        fpritf_uint_matrix(stderr,SIMPLEx_test->bIDx,1,m);

        fprintf(stderr,"\nPrinting Zcnb:\n");
        fpritf_dbl_matrix(stderr,SIMPLEx_test->Zcnb,1,n-m);

        fprintf(stderr,"\nPrinting Zcb:\n");
        fpritf_dbl_matrix(stderr,SIMPLEx_test->Zcb,1,m); 
        
        fprintf(stderr,"\nPrinting Wcb:\n");
        fpritf_dbl_matrix(stderr,Wcb,1,m);

        fprintf(stderr,"\nPrinting Wcnb:\n");
        fpritf_dbl_matrix(stderr,Wcnb,1,n-m);

        fprintf(stderr,"\nPrinting Wb:\n");
        fprintf(stderr,"Wb: %f\n",Wb); */
        

        int res_lplib = solve(lplib_test);
        SIMPLEx_test->ret = simplex_revisedSolve(m,n, (double *)SIMPLEx_test->AT, (double *)SIMPLEx_test->b, 
                                            (uint16_t *)SIMPLEx_test->nbIDx, (uint16_t *)SIMPLEx_test->bIDx,
                                            (double *)SIMPLEx_test->Zcnb, (double *)SIMPLEx_test->Zcb, 
                                            (double *)SIMPLEx_test->Wcnb, (double *)SIMPLEx_test->Wcb, 
                                            SIMPLEx_test->Wb,&SIMPLEx_test->Zopt);
        
        if (check_SIMPLEx_lplib_solution(lplib_test,res_lplib,SIMPLEx_test,SIMPLEx_test->ret))
            fprintf(stderr,"Test match\n");
        else
        {
            fprintf(stderr,"Tests DONT match\n");
            fprintf(stderr,"\nSIMPLEx returned %d while lplib %d\n",SIMPLEx_test->ret,res_lplib);

            fprintf(stderr,"\nPrinting AT:\n");
            fpritf_dbl_matrix(stderr,SIMPLEx_test->AT,n,m);
            fprintf(stderr,"\nPrinting b:\n");
            fpritf_dbl_matrix(stderr,SIMPLEx_test->b,1,m);

            fprintf(stderr,"\nPrinting nbIDx:\n");
            fpritf_uint_matrix(stderr,SIMPLEx_test->nbIDx,1,n-m);

            fprintf(stderr,"\nPrinting bIDx:\n");
            fpritf_uint_matrix(stderr,SIMPLEx_test->bIDx,1,m);

            fprintf(stderr,"\nPrinting Zcnb:\n");
            fpritf_dbl_matrix(stderr,SIMPLEx_test->Zcnb,1,n-m);

            fprintf(stderr,"\nPrinting Zcb:\n");
            fpritf_dbl_matrix(stderr,SIMPLEx_test->Zcb,1,m);

            fprintf(stderr,"\nPrinting Zopt:\n");
            fprintf(stderr,"%f\n",SIMPLEx_test->Zopt);

            print_lp(lplib_test);
            print_solution(lplib_test,1);
            print_objective(lplib_test);
            getc(stdin);
        }
                    

        // free the memory
        free_PLtest(SIMPLEx_test);
        delete_lp(lplib_test);
        free(constr_row);
        free(Zcnb);
        free(Zcb);
        free(Wcnb);
        free(Wcb);
    }
}

static void fpritf_dbl_matrix(FILE *fp, double *mat, unsigned int rows, unsigned int cols)
{
    for (int j = 0; j < rows; j++)
    {
        fprintf(fp,"| ");
        for (int i = 0; i < cols; i++)
        {
            fprintf(fp,"%f ", mat[j*cols + i]);
        }
        fprintf(fp,"|\n");
    }
}

static void fpritf_uint_matrix(FILE *fp, uint16_t *mat, unsigned int rows, unsigned int cols)
{
    for (int j = 0; j < rows; j++)
    {
        fprintf(fp,"| ");
        for (int i = 0; i < cols; i++)
        {
            fprintf(fp,"%d ", mat[j*cols + i]);
        }
        fprintf(fp,"|\n");
    }
}

static bool check_SIMPLEx_lplib_solution(lprec *lplib_test, int lplib_ret_val, PLtest_t *SIMPLEx_test, int SIMPLEx_ret_val)
{
    // First check for ret vals that must be coherent
    if ((SIMPLEx_ret_val == 0 && lplib_ret_val != 0) || 
        (SIMPLEx_ret_val == -1 && lplib_ret_val != 3) || 
        (SIMPLEx_ret_val == -2 && lplib_ret_val != 2))
        return false;
    
    // If the problem is unlimited or infeasible and return values match than return true
    if (SIMPLEx_ret_val == -1 || SIMPLEx_ret_val == -2)
        return true;
    
    // Check for obj function value
    double lplib_test_obj_val = get_objective(lplib_test);

    if (SIMPLEx_test->min_max == false)
    {
        // min problem
        if (fabs(lplib_test_obj_val - SIMPLEx_test->Zopt) > FLT_EPSILON)
            return false;
    } 
    else
    {
        // max problem
        if (fabs(lplib_test_obj_val + SIMPLEx_test->Zopt) > FLT_EPSILON)
            return false;
    }

    // The following lines perform a check that sees if all the variable values
    // in the both the solutions match between each others.
    // It is commented since different algorithms could find different solutions
    // but they are still optimal.
    // e.g. when obj function level curve is aligned with one constraint:
    // in this case there are INFINITE different optimal solutions

    /* double *lplib_test_var_values;
    if(get_ptr_variables(lplib_test, &lplib_test_var_values) != 1)
        return false;

    uint16_t n_vars = get_Ncolumns(lplib_test);

    // if obj functions have the same value we test for optimal basic variables values
    for (uint16_t i = 0; i < SIMPLEx_test->m; i++)
    {
        if (SIMPLEx_test->bIDx[i] > n_vars-1)
            continue;
        
        if (fabs(SIMPLEx_test->b[i] - lplib_test_var_values[SIMPLEx_test->bIDx[i]]) > FLT_EPSILON)
            return false;
    } */
        
    return true;
}