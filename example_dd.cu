#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dbldbl.h"

// Macro to catch CUDA errors in CUDA runtime calls
#define CUDA_SAFE_CALL(call)                                          \
do {                                                                  \
    cudaError_t err = call;                                           \
    if (cudaSuccess != err) {                                         \
        fprintf (stderr, "Cuda error in file '%s' in line %i : %s.\n",\
                 __FILE__, __LINE__, cudaGetErrorString(err) );       \
        exit(EXIT_FAILURE);                                           \
    }                                                                 \
} while (0)

__global__ void solve_quadratic_eq (double a, double b, double c, double *res)
                              
{
    /* Compute solutions in double precision using standard quadratic formula*/
    res[0] = (-b + sqrt (b*b - 4.0*a*c)) / (2.0 * a);
    res[1] = (-b - sqrt (b*b - 4.0*a*c)) / (2.0 * a);

    /* Compute solutions in double-double using standard quadratic formula */
    dbldbl aa = make_dbldbl (a, 0.0);
    dbldbl bb = make_dbldbl (b, 0.0);
    dbldbl cc = make_dbldbl (c, 0.0);
    dbldbl four = make_dbldbl (4.0, 0.0);

    dbldbl zz = neg_dbldbl (bb);           // -b
    dbldbl yy = mul_double_to_dbldbl (b,b);// b*b
    dbldbl ww = mul_dbldbl (four, aa);     // 4*a
    dbldbl vv = mul_dbldbl (ww, cc);       // 4*a*c
    dbldbl uu = sub_dbldbl (yy, vv);       // b*b - 4*a*c
    dbldbl tt = sqrt_dbldbl (uu);          // sqrt (b*b - 4*a*c)
    dbldbl rr = add_double_to_dbldbl (a,a);// 2*a
    dbldbl qq = add_dbldbl (zz, tt);       // -b + sqrt (b*b - 4*a*c)
    dbldbl pp = sub_dbldbl (zz, tt);       // -b - sqrt (b*b - 4*a*c)
    dbldbl xx1 = div_dbldbl (qq, rr);      // (-b + sqrt (b*b - 4*a*c)) / (2*a)
    dbldbl xx2 = div_dbldbl (pp, rr);      // (-b - sqrt (b*b - 4*a*c)) / (2*a)

    res[2] = get_dbldbl_head(xx1) + get_dbldbl_tail(xx1);
    res[3] = get_dbldbl_head(xx2) + get_dbldbl_tail(xx2);

    /* Compute solutions in double precision using more robust formula */
    double q = -0.5 * (b + copysign (sqrt (b*b - 4.0*a*c), b));
    res[4] = q / a;
    res[5] = c / q;
}

int main (void)
{
    /* Naive computation of the solution of a quadratic equation using both 
       double precision and double-double computation, using an example from
       George E. Forsythe, How Do You Solve a Quadratic Equation, Technical
       Report No. CS40, Computer Science Department, Stanford University, 
       June 1966.
    */
    double a = 1.0;
    double b = -100000.0;
    double c = 1.0;
    double *res = 0;
    double x1d, x2d, x1dd, x2dd, x1r, x2r;

    printf ("\nSolving quadratic equation with a = %g  b = %g  c = %g\n", 
            a, b, c);

    CUDA_SAFE_CALL (cudaMalloc ((void**)&res, 6*sizeof(double)));
    solve_quadratic_eq<<<1,1>>>(a, b, c, res);
    CUDA_SAFE_CALL (cudaThreadSynchronize());
    CUDA_SAFE_CALL (cudaMemcpy (&x1d, &res[0], sizeof(x1d),
                                cudaMemcpyDeviceToHost));
    CUDA_SAFE_CALL (cudaMemcpy (&x2d, &res[1], sizeof(x2d),
                                cudaMemcpyDeviceToHost));   
    CUDA_SAFE_CALL (cudaMemcpy (&x1dd, &res[2], sizeof(x1dd),
                                cudaMemcpyDeviceToHost));   
    CUDA_SAFE_CALL (cudaMemcpy (&x2dd, &res[3], sizeof(x2dd),
                                cudaMemcpyDeviceToHost));
    CUDA_SAFE_CALL (cudaMemcpy (&x1r, &res[4], sizeof(x1r),
                                cudaMemcpyDeviceToHost));   
    CUDA_SAFE_CALL (cudaMemcpy (&x2r, &res[5], sizeof(x2r),
                                cudaMemcpyDeviceToHost));

    printf ("\nUsing double precision (std. quadratic formula):\n");
    printf ("x1 =% 18.11e   a*x1**2+b*x1+c =% 18.11e\n",
            x1d, a*x1d*x1d+b*x1d+c);
    printf ("x2 =% 18.11e   a*x2**2+b*x2+c =% 18.11e\n", 
            x2d, a*x2d*x2d+b*x2d+c);

    printf ("\nUsing double-double (std. quadratic formula):\n");
    printf ("x1 =% 18.11e   a*x1**2+b*x1+c =% 18.11e\n", 
            x1dd, a*x1dd*x1dd+b*x1dd+c);
    printf ("x2 =% 18.11e   a*x2**2+b*x2+c =% 18.11e\n",
            x2dd, a*x2dd*x2dd+b*x2dd+c);

    printf ("\nUsing double precision (more robust formula):\n");
    printf ("x1 =% 18.11e   a*x1**2+b*x1+c =% 18.11e\n",
            x1r, a*x1r*x1r+b*x1r+c);
    printf ("x2 =% 18.11e   a*x2**2+b*x2+c =% 18.11e\n", 
            x2r, a*x2r*x2r+b*x2r+c);

    CUDA_SAFE_CALL (cudaFree(res));
    
    return EXIT_SUCCESS;
}
