/****************************************************************

file festkomma.h

Elkies package for the three cubes problem,
Version 1.0

*****************************************************************/


#include<gmp.h>

typedef unsigned long int UDItype;
typedef long int           DItype;


#define MAX(h, i) ((h) > (i) ? (h) : (i))
#define MIN(h, i) ((h) < (i) ? (h) : (i))

/* The g ++ compiler relies on the normal macros add_ssaaaa and subf_ddmmss
   not because there is a strange cast to be made while writing. */
#define addf_ssaaaa(sh, sl, ah, al, bh, bl)                       \
  __asm__ ("addq %5,%1\n\tadcq %3,%0"                             \
           : "=r" (sh),           "=&r" (sl)                      \
           : "0"  ((UDItype)(ah)), "g" ((UDItype)(bh)),           \
             "%1" ((UDItype)(al)), "g" ((UDItype)(bl)))

#define subf_ddmmss(sh, sl, ah, al, bh, bl)                       \
  __asm__ ("subq %5,%1\n\tsbbq %3,%0"                             \
           : "=r" (sh), "=&r" (sl)                                \
           : "0" ((UDItype)(ah)), "g" ((UDItype)(bh)),            \
             "1" ((UDItype)(al)), "g" ((UDItype)(bl)))

#define umul_ppmm(w1, w0, u, v)                                   \
  __asm__ ("mulq %3"                                              \
	   : "=a" (w0), "=d" (w1)                                 \
	   : "%0" ((UDItype)(u)), "rm" ((UDItype)(v)))

#define smul_ppmm(w1, w0, u, v)                                   \
  __asm__ ("imulq %3"                                             \
	   : "=a" (w0), "=d" (w1)                                 \
	   : "%0" ((DItype)(u)), "rm" ((DItype)(v)))


/* Can encode two types of fixed-point numbers.
   1. A number between 0 and 1 with precision 2 * (- 128) according to the rule
          fixed [1] * 2 ** (- 64) + fixed [0] * 2 ** (- 128).
      Data type "mpx_t".
   2. A signed number between (-2 ** 63) and 2 ** 63 with accuracy
      2 * (- 64) according to the rule
          (signed long) fixed [1] + fixed [0] * 2 ** (- 64).
      Data type "mpxg_t".*/
typedef unsigned long int mpx_t [2];


/* The naive conversion of mpxg_t into double.*/
#define mpxg_get_d_simple(fixed)                                  \
  (ldexp (((double) (fixed)[0]), -64)                             \
+ ((double) ((signed long) (fixed)[1])))


/* Conversion of mpxg_t into double.
   For values ​​close to 0, you want very high accuracy.
   Even at those just under 0, as
         (-1) + fixed [0] * 2 ** (- 64)
   are. */
inline double mpxg_get_d (mpx_t fixed) {
 /* fixed[0] >= 2**63 <==> gebrochener Teil >= 0.5. */
 if ((signed long) fixed[0] < 0) {
  /* Artfully complicated conversion into double.
     The first addend is [erg] - 1, which may be small
     and can be converted into double exactly. */
  return (- ldexp (((double) -fixed[0]), -64)
                 + ((double) (1 + (signed long) fixed[1])));
 }
 /* The naive conversion into double.
  fixed [1] should be understood as signed. */
 return (   ldexp (((double) fixed[0]), -64) 
                 + ((double) ((signed long) fixed[1])));
}


/* The naive conversion of mpx_t into double. */
#define mpx_get_d(fixed)                                          \
   (ldexp (((double) (fixed)[1]), -64)                            \
 + ldexp (((double) (fixed)[0]), -128))


/* requirement: 0 < floa < 1. */
inline void mpx_set_d (mpx_t erg, double floa) {
 floa = ldexp (floa, 64);
 erg[1] = (ulong) floor (floa);
 floa -= erg[1];
 floa = ldexp (floa, 64);
 erg[0] = (ulong) floor (floa);
 /* Hier darf kein lround stehen. Gibt Probleme bei Zahlen >2**63! */
}


/* requirement: 0 < floa < 1. Precision von fl >=128 Bit. 
   Attention: extremely inefficient! */
inline void mpx_set_mpf (mpx_t erg, mpf_t fl) {
 /* Efficient code looks like this:
 mp_ptr  ptr;
 long    pos;
 ptr = PTR (fl);
 pos = ABSIZ (fl) - 1;
 erg[1] = ptr[pos];
 erg[0] = ptr[pos - 1]; */

 mpf_t  tmp;

 mpf_init (tmp);
 mpf_mul_2exp (tmp, fl, 64);
 erg[1] = mpf_get_ui (tmp);
 mpf_sub_ui (tmp, tmp, erg[1]);
 mpf_mul_2exp (tmp, tmp, 64);
 erg[0] = mpf_get_ui (tmp);
 mpf_clear (tmp);
}


/* Intended for output with gmp_printf. Has certainly optimization potential. */
inline void mpf_set_mpx (mpf_t erg, mpx_t fixed) {
 mpf_set_ui (erg, fixed[0]);
 mpf_div_2exp (erg, erg, 64);
 mpf_add_ui (erg, erg, fixed[1]);
 mpf_div_2exp (erg, erg, 64);
}


/* Product of two mpx_t numbers (128-bit fixed-point numbers between 0 and 1).
   Not for data type mpxg_t. */
inline void mpx_mul (mpx_t erg, mpx_t fak1, mpx_t fak2) {
 ulong  argh, argl;

 /* Main part of the product. */
 umul_ppmm (erg[1], erg[0], fak1[1], fak2[1]);
 /* Multiplication over cross. Ignore the third limb. */
 umul_ppmm (argh, argl, fak1[1], fak2[0]);
 addf_ssaaaa (erg[1], erg[0], erg[1], erg[0], 0L, argh);
 /* Multiplication over cross the other way round. Ignore again the third limb. */
 umul_ppmm (argh, argl, fak1[0], fak2[1]);
 addf_ssaaaa (erg[1], erg[0], erg[1], erg[0], 0L, argh);
 /* Leave the product of low-end limbs completely gone. */
}


/* Unsigned fixed point * long.
   fix is ​​unsigned fixed-point, that is, of the data type mpx_t.
   The result erg is a signed fixed point of the data type mpxg_t.*/
inline void mpx_mul_si (mpx_t erg, mpx_t fix, long ganz) {
 ulong  argh, argl;

 /* Treat as unsigned. */
 umul_ppmm (erg[1], erg[0], fix[1], ganz);
 umul_ppmm (argh, argl, fix[0], ganz);
 /* Ignore argl. */
 addf_ssaaaa (erg[1], erg[0], erg[1], erg[0], 0, argh);

 /* Correction if completely <0.
     We counted on quite an error of exactly 2 ** 64. */
 if (ganz < 0)
  subf_ddmmss (erg[1], erg[0], erg[1], erg[0], fix[1], fix[0]);
}


/* Works with mpx_t as mpxg_t. */
#define mpx_add(summe, summand1, summand2)                        \
 addf_ssaaaa (summe[1], summe[0], summand1[1], summand1[0],       \
                                  summand2[1], summand2[0])


/* Works with mpx_t as mpxg_t. */
#define mpx_sub(differenz, minuend, subtrahend)                   \
 subf_ddmmss (differenz[1], differenz[0], minuend[1], minuend[0], \
                                    subtrahend[1], subtrahend[0])


/* Works on summand1 of the data type mpx_t or mpxg_t. summand2 is a long one.
   Overflow is also here, as always with add and sub, in the application
   to pay attention. */
#define mpxb_add_si(summe, summand1, summand2)                    \
 addf_ssaaaa (summe[1], summe[0], summand1[1], summand1[0],       \
                                                    summand2, 0)
