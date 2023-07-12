#include <stdint.h>
#include <stdio.h>

/*
* Execute the revised simplex algorithm minimizing the input linear programming problem in the canonical form:
* 
*   INPUT PARAMETERS:
*       - m:     the number of constraints of the problem
*       - n:     the number of variables of the system of equations in canonical form 
*                   (n also counts slack and artificial variables)
*       - NBT:   the (n-m) X m matrix of the coefficients of the non-basic variables of
*                   the system of equations as columns. It has to be considered as the 
*                   transposed of the classic NB
*       - B:     the m X m matrix of the coefficients of the basic variables of the system 
*                   of equations
*       - Cnb:   the 1 X (n-m) vector of the coefficients of the non-basic variables in the 
*                   objective function
*       - nbIDx: the 1 X (n-m) vector of the indices of the non-basic variables in the 
*                   current solution
*       - Cb:    the 1 X m vector of the coefficients of the basic variables in the 
*                   objective function
*       - bIDx:  the 1 X m vector of the indices of the basic variables in the current solution
*       - b:     the 1 X m vector of the values of the constant terms of the system of equations
*    
*   RETURNS:
*       - (0) if the algorithm terminated with success: in this case bIDx contains the indices 
*           of the basic variables of the solution and b contains the values of these variables.
*           The value of the objective function is found by calculating z = Cb*B*b
*       - (-1) if the problem is unlimited
*       - (-2) if memory allocation problems occurs
* 
*   IMPORTANT NOTES: 
*       - All the input matrices and vectors MUST be allocated in memory-contiguous locations 
*           and values stored in row-major order
*       - All b elements have to be non-negative because the algorithm assumes the system of 
*           equations is in the canonical form
*       - All values x in nbIDx and bIDx must satisfy x>=0 and x<n
*/
int simplex_solve(uint16_t m, uint16_t n, double *NBT, double *B, double *Cnb, uint16_t *nbIDx, 
                                                        double *Cb, uint16_t *bIDx, double *b);

void simplex_saveToStream(FILE *st, uint16_t m, uint16_t n, double *NBT, double *B, double *Cnb, 
                                        uint16_t *nbIDx, double *Cb, uint16_t *bIDx, double *b);